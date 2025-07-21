import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# --- Mathematical Models for Fitting ---
def binding_curve(log_conc, bottom, top, log_kd, hill_slope):
    """
    Hill binding model for log-transformed concentrations.
    Used for direct binding affinity (Kd).
    """
    return bottom + (top - bottom) / (1 + 10 ** ((log_kd - log_conc) * hill_slope))


def kd_model(c_ligand, unbound, bound, kd, c_target):
    """
    Quadratic binding model (Kd model) accounting for target depletion.
    c_ligand and kd should be in Molar.
    """
    c_ligand = np.asarray(c_ligand)
    c_target = np.asarray(c_target)

    # Ensure c_target is not zero to avoid division by zero
    if np.any(c_target == 0):
        logging.warning("Target concentration (c_target) cannot be zero for Kd model.")
        # Return an array of NaN or a suitable default to prevent errors
        return np.full_like(c_ligand, np.nan)

    # Calculate the discriminant, ensuring it's non-negative for the square root
    discriminant = (c_ligand + c_target + kd) ** 2 - 4 * c_ligand * c_target
    numerator = (c_ligand + c_target + kd) - np.sqrt(np.maximum(0, discriminant))

    # Calculate fraction bound, ensuring it's within [0, 1]
    fraction_bound = numerator / (2 * c_target)
    fraction_bound = np.clip(fraction_bound, 0, 1)
    return unbound + (bound - unbound) * fraction_bound


def dose_response_curve(log_conc, bottom, top, log_ec50, hill_slope):
    """
    4-Parameter Logistic (4PL) model for dose-response curves.
    log_conc is log10(concentration).
    log_ec50 is log10(EC50).
    """
    return bottom + (top - bottom) / (1 + 10 ** ((log_ec50 - log_conc) * hill_slope))


# --- GUI Functions ---
class MSTAnalyzerGUI:
    def __init__(self, master):
        self.master = master
        master.title("MST Batch Analysis (Interactive)")
        master.geometry("1400x900")  # Increased size for plot

        self.input_dir = ""
        self.ligand_info = {}
        self.f0_start = -1.5
        self.f0_end = 0.0
        self.f1_start = 2.5
        self.f1_end = 4.0
        self.c_target = 1.0e-6  # Default 1 µM
        self.model_choice = tk.StringVar(value="Dose-Response")
        self.display_unit = tk.StringVar(value="nM")  # Default display unit for Kd/EC50

        # Stores all processed data points for interactive use:
        # {ligand_group_key:
        #   {'file_name':
        #       {'concentrations': [...], 'fnorms': [...], 'log_concs': [...],
        #        'excluded_indices': set(), 'F0': ..., 'F1': ..., 'Times': [...], 'Ratios': [...]}
        #   }
        # }
        self.all_data_points = {}

        # **MODIFIED**: Stores which files are selected, nested by ligand group for persistence.
        # {ligand_group_key: {file_name: tk.BooleanVar}}
        self.replicate_selection_vars = {}

        self.create_widgets()
        self.create_plot_area()

    def create_widgets(self):
        # Left Panel for Controls
        # Increased width to 5 to accommodate all labels/entries better
        control_frame = ttk.Frame(self.master, width=550)
        control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        control_frame.grid_propagate(False)  # Prevent frame from resizing to content

        # Frame for folder path
        frame_folder = ttk.LabelFrame(control_frame, text="Input Folder")
        frame_folder.grid(row=0, column=0, columnspan=2, padx=5, pady=5, sticky="ew")
        frame_folder.columnconfigure(1, weight=1)  # Allow entry field to expand

        ttk.Label(frame_folder, text="Folder Path:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.folder_path_entry = ttk.Entry(frame_folder, width=40)
        self.folder_path_entry.grid(row=0, column=1, padx=5, pady=5, sticky="ew")  # Sticky for expansion
        ttk.Button(frame_folder, text="Browse", command=self.browse_folder).grid(row=0, column=2, padx=5, pady=5)

        # Frame for ligand info and exclusion (dynamic part)
        frame_ligand_data = ttk.LabelFrame(control_frame, text="Ligand/Replicate Groups")
        frame_ligand_data.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="ew")
        # Configure columns to expand for input fields
        frame_ligand_data.columnconfigure(0, weight=1)
        frame_ligand_data.columnconfigure(1, weight=1)
        frame_ligand_data.columnconfigure(2, weight=1)

        # Using a Canvas with a scrollbar for ligand rows if many are added
        ligand_canvas = tk.Canvas(frame_ligand_data, height=150)  # Fixed height for the canvas
        ligand_scrollbar = ttk.Scrollbar(frame_ligand_data, orient="vertical", command=ligand_canvas.yview)
        self.ligand_rows_frame = ttk.Frame(ligand_canvas)

        self.ligand_rows_frame.bind(
            "<Configure>",
            lambda e: ligand_canvas.configure(
                scrollregion=ligand_canvas.bbox("all")
            )
        )
        ligand_canvas.create_window((0, 0), window=self.ligand_rows_frame, anchor="nw")
        ligand_canvas.configure(yscrollcommand=ligand_scrollbar.set)

        ligand_canvas.grid(row=0, column=0, columnspan=3, sticky="nsew", padx=5, pady=5)
        ligand_scrollbar.grid(row=0, column=3, sticky="ns", padx=0, pady=5)
        frame_ligand_data.rowconfigure(0, weight=1)  # Allow the canvas row to expand vertically
        frame_ligand_data.columnconfigure(0, weight=1)  # Allow canvas column to expand horizontally

        self.ligand_entries = []  # (key_entry, max_conc_entry, df_entry)

        ligand_header_row = 0
        # Changed label to reflect the stricter matching
        ttk.Label(self.ligand_rows_frame, text="Filename Identifier").grid(row=ligand_header_row, column=0, padx=2,
                                                                           pady=2, sticky="w")
        ttk.Label(self.ligand_rows_frame, text="Max Conc (mM)").grid(row=ligand_header_row, column=1, padx=2, pady=2,
                                                                     sticky="w")
        ttk.Label(self.ligand_rows_frame, text="Dilution Factor").grid(row=ligand_header_row, column=2, padx=2, pady=2,
                                                                       sticky="w")

        # Make header columns expandable too
        self.ligand_rows_frame.columnconfigure(0, weight=1)
        self.ligand_rows_frame.columnconfigure(1, weight=1)
        self.ligand_rows_frame.columnconfigure(2, weight=1)

        self.add_ligand_button = ttk.Button(frame_ligand_data, text="Add Ligand Group", command=self.add_ligand_row)
        self.add_ligand_button.grid(row=1, column=0, columnspan=4, pady=5)  # Placed below canvas

        self.add_ligand_row()  # Add one row by default

        # Frame for Analysis Parameters
        frame_params = ttk.LabelFrame(control_frame, text="Analysis Parameters")
        frame_params.grid(row=2, column=0, columnspan=2, padx=5, pady=5, sticky="ew")
        frame_params.columnconfigure(1, weight=1)  # Allow parameter entries to expand
        frame_params.columnconfigure(3, weight=1)

        ttk.Label(frame_params, text="Model Choice:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        ttk.Radiobutton(frame_params, text="Hill Binding (logKd)", variable=self.model_choice, value="Hill",
                        command=self.update_params).grid(row=0, column=1, sticky="w", padx=2, pady=5)
        ttk.Radiobutton(frame_params, text="Kd (Quadratic)", variable=self.model_choice, value="Kd",
                        command=self.update_params).grid(row=0, column=2, sticky="w", padx=2, pady=5)
        ttk.Radiobutton(frame_params, text="Dose-Response (4PL)", variable=self.model_choice, value="Dose-Response",
                        command=self.update_params).grid(row=0, column=3, sticky="w", padx=2, pady=5)

        ttk.Label(frame_params, text="F0 Window (s):").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.f0_start_entry = ttk.Entry(frame_params, width=10)
        self.f0_start_entry.grid(row=1, column=1, padx=2, pady=5, sticky="ew")
        self.f0_start_entry.insert(0, str(self.f0_start))
        self.f0_start_entry.bind("<Return>", self.update_params)
        self.f0_start_entry.bind("<FocusOut>", self.update_params)

        ttk.Label(frame_params, text="to").grid(row=1, column=2, padx=0, pady=5)
        self.f0_end_entry = ttk.Entry(frame_params, width=10)
        self.f0_end_entry.grid(row=1, column=3, padx=2, pady=5, sticky="ew")
        self.f0_end_entry.insert(0, str(self.f0_end))
        self.f0_end_entry.bind("<Return>", self.update_params)
        self.f0_end_entry.bind("<FocusOut>", self.update_params)

        ttk.Label(frame_params, text="F1 Window (s):").grid(row=2, column=0, sticky="w", padx=5, pady=5)
        self.f1_start_entry = ttk.Entry(frame_params, width=10)
        self.f1_start_entry.grid(row=2, column=1, padx=2, pady=5, sticky="ew")
        self.f1_start_entry.insert(0, str(self.f1_start))
        self.f1_start_entry.bind("<Return>", self.update_params)
        self.f1_start_entry.bind("<FocusOut>", self.update_params)

        ttk.Label(frame_params, text="to").grid(row=2, column=2, padx=0, pady=5)
        self.f1_end_entry = ttk.Entry(frame_params, width=10)
        self.f1_end_entry.grid(row=2, column=3, padx=2, pady=5, sticky="ew")
        self.f1_end_entry.insert(0, str(self.f1_end))
        self.f1_end_entry.bind("<Return>", self.update_params)
        self.f1_end_entry.bind("<FocusOut>", self.update_params)

        self.c_target_label = ttk.Label(frame_params, text="Target Conc (mM):")
        self.c_target_entry = ttk.Entry(frame_params, width=15)
        self.c_target_entry.insert(0, str(self.c_target * 1e3))  # Display in mM
        self.c_target_entry.bind("<Return>", self.update_params)
        self.c_target_entry.bind("<FocusOut>", self.update_params)

        self.toggle_c_target_visibility()  # Set initial visibility

        # Unit Selection for Kd/EC50
        ttk.Label(frame_params, text="Result Unit:").grid(row=4, column=0, sticky="w", padx=5, pady=5)
        self.unit_selector = ttk.Combobox(frame_params, textvariable=self.display_unit,
                                          values=["M", "mM", "µM", "nM", "pM"], state="readonly")
        self.unit_selector.grid(row=4, column=1, columnspan=3, padx=5, pady=5, sticky="ew")
        self.unit_selector.bind("<<ComboboxSelected>>", self.update_plot)

        # Frame for Replicate Selection
        self.frame_replicate_select = ttk.LabelFrame(control_frame, text="Select Replicates & Exclude Points")
        self.frame_replicate_select.grid(row=3, column=0, columnspan=2, padx=5, pady=5, sticky="ew")
        self.frame_replicate_select.columnconfigure(1, weight=1)  # Make combobox expand

        ttk.Label(self.frame_replicate_select, text="Select Ligand Group:").grid(row=0, column=0, sticky="w", padx=5,
                                                                                 pady=5)
        self.ligand_group_selector = ttk.Combobox(self.frame_replicate_select, state="readonly")
        self.ligand_group_selector.grid(row=0, column=1, sticky="ew", padx=5, pady=5)
        self.ligand_group_selector.bind("<<ComboboxSelected>>", self.load_replicate_checkboxes)

        # Using a Canvas with a scrollbar for replicate checkboxes
        replicate_canvas = tk.Canvas(self.frame_replicate_select, height=150)  # Fixed height
        replicate_scrollbar = ttk.Scrollbar(self.frame_replicate_select, orient="vertical",
                                            command=replicate_canvas.yview)
        self.replicate_checkbox_frame = ttk.Frame(replicate_canvas)

        self.replicate_checkbox_frame.bind(
            "<Configure>",
            lambda e: replicate_canvas.configure(
                scrollregion=replicate_canvas.bbox("all")
            )
        )
        replicate_canvas.create_window((0, 0), window=self.replicate_checkbox_frame, anchor="nw")
        replicate_canvas.configure(yscrollcommand=replicate_scrollbar.set)

        replicate_canvas.grid(row=1, column=0, columnspan=2, sticky="nsew", padx=5, pady=5)
        replicate_scrollbar.grid(row=1, column=2, sticky="ns", padx=0, pady=5)
        self.frame_replicate_select.rowconfigure(1, weight=1)  # Allow the canvas row to expand vertically

        self.results_label = ttk.Label(control_frame, text="Results:", font=("Arial", 12, "bold"))
        self.results_label.grid(row=4, column=0, columnspan=2, padx=5, pady=10, sticky="w")

        # Export Main Plot Data Button
        self.export_main_button = ttk.Button(control_frame, text="Export Current Plot Data",
                                             command=self.export_main_plot_data)
        self.export_main_button.grid(row=5, column=0, columnspan=2, pady=5)

        # --- NEW: Overlay Plot Section ---
        self.frame_overlay_plot = ttk.LabelFrame(control_frame, text="Overlay Two Ligand Groups")
        self.frame_overlay_plot.grid(row=6, column=0, columnspan=2, padx=5, pady=5, sticky="ew")  # Adjusted row
        self.frame_overlay_plot.columnconfigure(1, weight=1)

        ttk.Label(self.frame_overlay_plot, text="Group 1:").grid(row=0, column=0, sticky="w", padx=5, pady=5)
        self.overlay_group1_selector = ttk.Combobox(self.frame_overlay_plot, state="readonly")
        self.overlay_group1_selector.grid(row=0, column=1, sticky="ew", padx=5, pady=5)

        ttk.Label(self.frame_overlay_plot, text="Group 2:").grid(row=1, column=0, sticky="w", padx=5, pady=5)
        self.overlay_group2_selector = ttk.Combobox(self.frame_overlay_plot, state="readonly")
        self.overlay_group2_selector.grid(row=1, column=1, sticky="ew", padx=5, pady=5)

        self.generate_overlay_button = ttk.Button(self.frame_overlay_plot, text="Generate Overlay Plot",
                                                  command=self.generate_overlay_plot)
        self.generate_overlay_button.grid(row=2, column=0, columnspan=2, pady=5)

    def create_plot_area(self):
        # Right Panel for Plotting
        plot_frame = ttk.Frame(self.master)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

        self.fig, self.ax = plt.subplots(figsize=(8, 6))  # Adjust figsize as needed
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        self.toolbar.update()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.ax.set_xscale('log')
        self.ax.set_xlabel('Ligand Concentration [M]')
        self.ax.set_ylabel('Fnorm')
        self.ax.set_title('Interactive Binding Curve')
        self.ax.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.7)

        # Connect the pick event
        self.canvas.mpl_connect('pick_event', self.on_pick)

    def on_button_press(self, event):
        # This is for general clicks, can be used for things like right-click menu or
        # clearing selected points if nothing specific is picked.
        if event.button == 3:  # Right-click
            # Example: Could clear exclusions or open a context menu
            pass

    def on_pick(self, event):
        # The pick event can be triggered by various artists. We are interested in scatter points.
        # Check if the picked artist is a PathCollection (from scatter plot)
        if not hasattr(event.artist, 'get_offsets'):
            return  # Not a scatter plot point

        # Get the index of the clicked point within the artist's data
        ind = event.ind[0]  # event.ind is an array of indices if multiple points are picked, take the first one

        # Get the coordinates of the picked point
        offsets = event.artist.get_offsets()
        clicked_conc = offsets[ind, 0]  # X-coordinate (concentration)
        clicked_fnorm = offsets[ind, 1]  # Y-coordinate (Fnorm)

        logging.info(f"Picked point details: X={clicked_conc:.2e}, Y={clicked_fnorm:.2f}, Index in Artist: {ind}")

        # Determine which replicate and which data point it is
        current_ligand_group = self.ligand_group_selector.get()
        if not current_ligand_group:
            logging.warning("No ligand group selected in main plot for picking.")
            return

        # **MODIFIED**: Use the correct nested dictionary for replicate selections
        group_selections = self.replicate_selection_vars.get(current_ligand_group, {})

        found = False
        # Iterate through the files (replicates) in the current ligand group
        for file_name, file_data in self.all_data_points.get(current_ligand_group, {}).items():
            # Only consider selected replicates for the main plot's interactive exclusion
            if not group_selections.get(file_name, tk.BooleanVar(value=False)).get():
                continue

            # Find the exact point in the original data structure for this file
            # We use np.isclose for floating point comparison accuracy
            # It's important to iterate over the *original* indices (0 to len-1)
            # as these are what we store in excluded_indices
            for i, (conc, fnorm) in enumerate(zip(file_data['concentrations'], file_data['fnorms'])):
                if np.isclose(conc, clicked_conc, rtol=1e-9) and np.isclose(fnorm, clicked_fnorm, rtol=1e-9):
                    logging.info(f"Found match in file: {os.path.basename(file_name)} at original index {i}")
                    logging.info(f"Excluded indices BEFORE: {file_data['excluded_indices']}")
                    if i in file_data['excluded_indices']:
                        # Point is currently excluded, so include it
                        file_data['excluded_indices'].remove(i)
                        logging.info(f"Included point {i} ({conc:.2e} M) from {os.path.basename(file_name)}")
                    else:
                        # Point is currently included, so exclude it
                        file_data['excluded_indices'].add(i)
                        logging.info(f"Excluded point {i} ({conc:.2e} M) from {os.path.basename(file_name)}")
                    logging.info(f"Excluded indices AFTER: {file_data['excluded_indices']}")
                    found = True
                    break  # Found the point in this file, no need to check other points in this file
            if found:
                break  # Found the point in this file, no need to check other files

        if not found:
            logging.warning(
                f"Could not find exact matching data point for picked coordinates ({clicked_conc:.2e}, {clicked_fnorm:.2f}). This might happen if axes limits or data processing subtly changes values.")

        # Re-draw the plot and re-fit the curve with the new exclusion
        self.update_plot()

    def browse_folder(self):
        folder_selected = filedialog.askdirectory(title="Select folder containing MST .csv files")
        if folder_selected:
            self.folder_path_entry.delete(0, tk.END)
            self.folder_path_entry.insert(0, folder_selected)
            self.input_dir = folder_selected
            self.load_all_data()  # Load all data once after folder selection

    def add_ligand_row(self):
        current_row = len(self.ligand_entries)

        key_entry = ttk.Entry(self.ligand_rows_frame, width=20)
        key_entry.grid(row=current_row + 1, column=0, padx=2, pady=2, sticky="ew")  # +1 for header row
        key_entry.insert(0, f"Ligand_Group_{current_row + 1}")
        key_entry.bind("<Return>", self.update_all_data_and_ui)
        key_entry.bind("<FocusOut>", self.update_all_data_and_ui)

        max_conc_entry = ttk.Entry(self.ligand_rows_frame, width=10)
        max_conc_entry.grid(row=current_row + 1, column=1, padx=2, pady=2, sticky="ew")
        max_conc_entry.insert(0, "100.0")
        max_conc_entry.bind("<Return>", self.update_all_data_and_ui)
        max_conc_entry.bind("<FocusOut>", self.update_all_data_and_ui)

        df_entry = ttk.Entry(self.ligand_rows_frame, width=10)
        df_entry.grid(row=current_row + 1, column=2, padx=2, pady=2, sticky="ew")
        df_entry.insert(0, "2.0")
        df_entry.bind("<Return>", self.update_all_data_and_ui)
        df_entry.bind("<FocusOut>", self.update_all_data_and_ui)

        self.ligand_entries.append((key_entry, max_conc_entry, df_entry))

        # Re-pack the button to be below the last added row
        self.add_ligand_button.grid(row=len(self.ligand_entries) + 1, column=0, columnspan=4, pady=5)
        # Update canvas scroll region after adding new widgets
        self.ligand_rows_frame.update_idletasks()
        self.ligand_rows_frame.master.config(scrollregion=self.ligand_rows_frame.master.bbox("all"))

    def update_params(self, event=None):
        """Update analysis parameters and trigger a re-plot/re-fit."""
        try:
            self.f0_start = float(self.f0_start_entry.get())
            self.f0_end = float(self.f0_end_entry.get())
            self.f1_start = float(self.f1_start_entry.get())
            self.f1_end = float(self.f1_end_entry.get())

            # Convert mM input to M for c_target
            self.c_target = float(self.c_target_entry.get()) * 1e-3

            if not (self.f0_start < self.f0_end and self.f1_start < self.f1_end):
                messagebox.showwarning("Input Error", "Window start time must be less than end time.")
                return
        except ValueError:
            messagebox.showwarning("Input Error", "Please enter valid numbers for time windows and concentration.")
            return

        self.toggle_c_target_visibility()

        # Recalculate Fnorms for all raw data based on new windows
        self.recalculate_all_fnorms()

        # Then update the plot with the new Fnorms
        self.update_plot()

    def toggle_c_target_visibility(self):
        if self.model_choice.get() == "Kd":
            self.c_target_label.grid(row=3, column=0, sticky="w", padx=5, pady=5)
            self.c_target_entry.grid(row=3, column=1, columnspan=2, padx=5, pady=5, sticky="ew")
        else:
            self.c_target_label.grid_remove()
            self.c_target_entry.grid_remove()

    def update_all_data_and_ui(self, event=None):
        """Called when ligand info changes, triggers full reload of data."""
        self.load_all_data()

    def load_all_data(self):
        """Loads all CSV data, calculates initial Fnorms, and populates data structures."""
        if not self.input_dir or not os.path.isdir(self.input_dir):
            return

        temp_ligand_info = {}
        for key_entry, max_conc_entry, df_entry in self.ligand_entries:
            ligand_key = key_entry.get().strip()
            if not ligand_key: continue
            try:
                max_conc = float(max_conc_entry.get()) * 1e-3  # Convert mM to M
                df_val = float(df_entry.get())
                if max_conc <= 0 or df_val <= 0: raise ValueError
                temp_ligand_info[ligand_key] = {'max_conc': max_conc, 'df': df_val}
            except ValueError:
                messagebox.showwarning("Input Error", f"Invalid input for {ligand_key} (Max Conc/Dilution Factor).")
                return

        self.ligand_info = temp_ligand_info

        # Preserve excluded indices if possible when reloading, but reset if file/group changes
        new_all_data_points = {}

        csv_files = [f for f in os.listdir(self.input_dir) if f.endswith(".csv")]
        if not csv_files:
            messagebox.showinfo("No Files", f"No .csv files found in the selected directory: {self.input_dir}")
            return

        for file_name in csv_files:
            filepath = os.path.join(self.input_dir, file_name)

            # Find the matching ligand key (group) using startsWith for stricter matching
            ligand_group_key = None
            for key in self.ligand_info:
                # Use os.path.basename to match against just the filename without path
                # Use .lower() for case-insensitive matching
                if os.path.basename(file_name).lower().startswith(key.lower()):
                    ligand_group_key = key
                    break  # Assign to the first matching group (prioritize order in GUI)

            if ligand_group_key is None:
                logging.warning(f"Skipping {file_name} — no matching ligand identifier found in '{self.ligand_info}'.")
                continue

            try:
                df = pd.read_csv(filepath, sep=';', encoding='utf-8-sig')
                # print(f"Columns found in {file_name}: {df.columns.tolist()}") # For debugging purposes

                # Determine the Fnorm/Ratio column name
                ratio_col_name = None
                possible_ratio_cols = ['Ratio 670nm / 650nm', 'Relative Fluorescence 670nm',
                                       'Relative Fluorescence 650nm']
                for col in possible_ratio_cols:
                    if col in df.columns:
                        ratio_col_name = col
                        break
                if ratio_col_name is None:
                    logging.warning(
                        f"Skipping {file_name} — missing expected ratio column (either 'Ratio 670nm / 650nm' , 'Relative Fluorescence 670nm' or 'Relative Fluorescence 650nm').")
                    continue

                # Ensure Capillary Position column exists and is numeric
                capillary_pos_col = 'Capillary Position'
                if capillary_pos_col not in df.columns:
                    logging.warning(f"Skipping {file_name} — missing required column: '{capillary_pos_col}'.")
                    continue

                df[capillary_pos_col] = pd.to_numeric(df[capillary_pos_col], errors='coerce')
                # Drop rows where capillary position could not be converted (e.g., empty or non-numeric)
                df.dropna(subset=[capillary_pos_col], inplace=True)

                required_cols = {'Time [s]', ratio_col_name, capillary_pos_col}
                if not required_cols.issubset(df.columns):
                    logging.warning(f"Skipping {file_name} — missing one or more required columns ({required_cols}).")
                    continue

                if ligand_group_key not in new_all_data_points:
                    new_all_data_points[ligand_group_key] = {}

                # Check if this file was already processed, preserve excluded indices if so
                existing_excluded_indices = set()
                if ligand_group_key in self.all_data_points and file_name in self.all_data_points[ligand_group_key]:
                    existing_excluded_indices = self.all_data_points[ligand_group_key][file_name]['excluded_indices']

                new_all_data_points[ligand_group_key][file_name] = {
                    'Times': df['Time [s]'].values,
                    'Ratios': df[ratio_col_name].values,  # Use the determined ratio column
                    'Capillary_Positions': df[capillary_pos_col].values,  # Use the determined capillary pos column
                    'concentrations': [], 'fnorms': [], 'log_concs': [],
                    'excluded_indices': existing_excluded_indices,  # Preserve here
                    'F0': np.nan, 'F1': np.nan  # Store these for debugging/future plots
                }

            except Exception as e:
                logging.error(f"Error loading {file_name}: {e}")
                continue

        self.all_data_points = new_all_data_points  # Update the main data structure
        self.recalculate_all_fnorms()  # Calculate Fnorms for all loaded data
        self.update_ligand_group_selector()  # Populate dropdowns
        self.load_replicate_checkboxes()  # Load checkboxes for selected ligand group

    def recalculate_all_fnorms(self):
        """Recalculates Fnorm values for all stored raw data based on current window settings."""
        for ligand_group_key, files_data in self.all_data_points.items():
            # Check if the ligand_group_key still exists in self.ligand_info
            if ligand_group_key not in self.ligand_info:
                logging.warning(f"Skipping Fnorm recalculation for {ligand_group_key}: Ligand info not found.")
                continue

            max_conc = self.ligand_info[ligand_group_key]['max_conc']
            dilution = self.ligand_info[ligand_group_key]['df']

            # Generate expected concentrations for this ligand group
            expected_concentrations = [max_conc / (dilution ** i) for i in range(16)]

            for file_name, data in files_data.items():
                temp_df = pd.DataFrame({
                    'Time [s]': data['Times'],
                    'Ratios': data['Ratios'],
                    'Capillary Position': data['Capillary_Positions']
                })

                unique_caps = sorted(temp_df['Capillary Position'].unique())
                cap_to_conc_idx = {cap: i for i, cap in enumerate(unique_caps)}

                current_file_fnorms = []
                current_file_concs = []
                current_file_log_concs = []

                for cap_pos in unique_caps:
                    cap_df = temp_df[temp_df['Capillary Position'] == cap_pos].copy()
                    conc_idx = cap_to_conc_idx.get(cap_pos)
                    if conc_idx is None or conc_idx >= len(expected_concentrations): continue
                    current_concentration = expected_concentrations[conc_idx]

                    F0_values = cap_df[(cap_df['Time [s]'] >= self.f0_start) & (cap_df['Time [s]'] <= self.f0_end)][
                        'Ratios']
                    F1_values = cap_df[(cap_df['Time [s]'] >= self.f1_start) & (cap_df['Time [s]'] <= self.f1_end)][
                        'Ratios']

                    F0 = F0_values.mean()
                    F1 = F1_values.mean()
                    Fnorm = F1 / F0 if F0 not in [np.nan, 0] else np.nan

                    if not np.isnan(Fnorm) and current_concentration > 0:
                        current_file_fnorms.append(Fnorm)
                        current_file_concs.append(current_concentration)
                        current_file_log_concs.append(np.log10(current_concentration))

                # Update the stored data for this file
                data['concentrations'] = np.array(current_file_concs)
                data['fnorms'] = np.array(current_file_fnorms)
                data['log_concs'] = np.array(current_file_log_concs)
                # Note: excluded_indices are NOT reset here, they persist across window changes.
                # If you want them reset, clear the set here: data['excluded_indices'] = set()

    def update_ligand_group_selector(self):
        """Populates the ligand group comboboxes (main and overlay)."""
        current_selection_main = self.ligand_group_selector.get()
        current_selection_o1 = self.overlay_group1_selector.get()
        current_selection_o2 = self.overlay_group2_selector.get()

        ligand_keys = list(self.all_data_points.keys())

        # Update main selector
        self.ligand_group_selector['values'] = ligand_keys
        if ligand_keys:
            if current_selection_main in ligand_keys:
                self.ligand_group_selector.set(current_selection_main)
            else:
                self.ligand_group_selector.set(ligand_keys[0])
            self.load_replicate_checkboxes()  # Load checkboxes for the newly selected main group
        else:
            self.ligand_group_selector.set("")
            self.clear_replicate_checkboxes()
            self.update_plot()  # Clear plot if no data

        # Update overlay selectors
        self.overlay_group1_selector['values'] = ligand_keys
        self.overlay_group2_selector['values'] = ligand_keys

        if ligand_keys:
            if current_selection_o1 in ligand_keys:
                self.overlay_group1_selector.set(current_selection_o1)
            elif len(ligand_keys) >= 1:
                self.overlay_group1_selector.set(ligand_keys[0])
            else:
                self.overlay_group1_selector.set("")  # No groups to select

            if current_selection_o2 in ligand_keys:
                self.overlay_group2_selector.set(current_selection_o2)
            elif len(ligand_keys) >= 2:  # Default to second group if available
                self.overlay_group2_selector.set(ligand_keys[1])
            elif len(ligand_keys) == 1:
                self.overlay_group2_selector.set("")  # Only one group, so no second option
            else:
                self.overlay_group2_selector.set("")  # No groups to select

    def clear_replicate_checkboxes(self):
        """Removes all replicate checkbox widgets from the frame."""
        for widget in self.replicate_checkbox_frame.winfo_children():
            widget.destroy()
        # **MODIFIED**: No longer clears the underlying data dictionary.

    def load_replicate_checkboxes(self, event=None):
        """Loads checkboxes for replicates of the selected ligand group, preserving their state."""
        self.clear_replicate_checkboxes()  # Clear existing checkbox widgets

        selected_group = self.ligand_group_selector.get()
        if not selected_group or selected_group not in self.all_data_points:
            return

        # **MODIFIED**: Ensure the nested dictionary for this group exists.
        self.replicate_selection_vars.setdefault(selected_group, {})
        group_selections = self.replicate_selection_vars[selected_group]

        replicates = self.all_data_points[selected_group]
        sorted_replicate_names = sorted(replicates.keys())

        for i, file_name in enumerate(sorted_replicate_names):
            # **MODIFIED**: Get the existing BooleanVar or create a new one (defaulting to True).
            # This is the core of the state-preservation fix.
            var = group_selections.setdefault(file_name, tk.BooleanVar(value=True))

            # The trace might be added multiple times if not handled carefully,
            # but BooleanVar handles this gracefully (replaces the old command).
            var.trace_add("write", lambda *args, fn=file_name: self.on_replicate_toggle(fn))

            chk = ttk.Checkbutton(self.replicate_checkbox_frame, text=os.path.basename(file_name), variable=var)
            chk.grid(row=i, column=0, sticky="w", padx=2, pady=2)

        # Update canvas scroll region after adding new widgets
        self.replicate_checkbox_frame.update_idletasks()
        self.replicate_checkbox_frame.master.config(scrollregion=self.replicate_checkbox_frame.master.bbox("all"))

        self.update_plot()  # Update plot after loading new checkboxes

    def on_replicate_toggle(self, file_name):
        """Called when a replicate checkbox is toggled."""
        self.update_plot()

    def _get_processed_group_data(self, ligand_group_key, use_selected_replicates=False):
        """
        Helper function to get aggregated and fitted data for a specific ligand group.
        :param ligand_group_key: The identifier for the ligand group.
        :param use_selected_replicates: If True, only use replicates currently selected
                                        in the main GUI's checkboxes for this group.
                                        Otherwise, use all non-excluded points.
        :return: A dictionary containing averaged data, fit results, and success status.
        """
        result = {
            'has_data': False,
            'avg_concs': np.array([]), 'avg_fnorms': np.array([]), 'avg_fnorms_errors': np.array([]),
            'fit_popt': None, 'fit_success': False, 'fit_comment': "No data for fit.",
            'kd_or_ec50_result': np.nan,
            'individual_replicates_data': []  # To store data for individual replicates for export
        }

        if ligand_group_key not in self.all_data_points:
            result['fit_comment'] = f"Ligand group '{ligand_group_key}' not found."
            return result

        all_concs_for_averaging = []
        all_fnorms_for_averaging = []

        files_data = self.all_data_points[ligand_group_key]

        # **MODIFIED**: Get the correct dictionary of selections for this group
        group_selections = self.replicate_selection_vars.get(ligand_group_key, {})

        for file_name, data in files_data.items():
            # **MODIFIED**: Check against the persistent selection state
            if use_selected_replicates and not group_selections.get(file_name, tk.BooleanVar(value=True)).get():
                continue  # Skip if not selected in main GUI's checkboxes

            current_replicate_concs = []
            current_replicate_fnorms = []
            current_replicate_excluded_concs = []
            current_replicate_excluded_fnorms = []

            for idx, (c, f) in enumerate(zip(data['concentrations'], data['fnorms'])):
                if idx not in data['excluded_indices']:  # Always respect individual point exclusions
                    all_concs_for_averaging.append(c)
                    all_fnorms_for_averaging.append(f)
                    current_replicate_concs.append(c)
                    current_replicate_fnorms.append(f)
                else:
                    current_replicate_excluded_concs.append(c)
                    current_replicate_excluded_fnorms.append(f)

            result['individual_replicates_data'].append({
                'file_name': os.path.basename(file_name),
                'Concentration': np.array(current_replicate_concs),
                'Fnorm': np.array(current_replicate_fnorms),
                'Excluded_Concentration': np.array(current_replicate_excluded_concs),
                'Excluded_Fnorm': np.array(current_replicate_excluded_fnorms)
            })

        if not all_concs_for_averaging:
            result['fit_comment'] = "No included data points for fitting in this group."
            return result

        reps_df = pd.DataFrame({
            'Concentration': all_concs_for_averaging,
            'Fnorm': all_fnorms_for_averaging
        })
        reps_df.dropna(subset=['Fnorm'], inplace=True)  # Drop any NaN Fnorm values before grouping

        if reps_df.empty:
            result['fit_comment'] = "No valid Fnorm data after filtering for this group."
            return result

        grouped_data = reps_df.groupby('Concentration')['Fnorm'].agg(['mean', 'std', 'count']).reset_index()
        grouped_data.rename(columns={'mean': 'Fnorm_mean', 'std': 'Fnorm_std', 'count': 'Replicate_Count'},
                            inplace=True)
        # Handle cases where Replicate_Count is 1 (std and sem would be NaN/Inf)
        grouped_data['Fnorm_sem'] = grouped_data.apply(
            lambda row: row['Fnorm_std'] / np.sqrt(row['Replicate_Count']) if row['Replicate_Count'] > 1 else np.nan,
            axis=1
        )

        avg_actual_concs = grouped_data['Concentration'].values
        avg_log_concs = np.log10(avg_actual_concs)
        avg_fnorms = grouped_data['Fnorm_mean'].values
        avg_fnorms_errors = grouped_data['Fnorm_sem'].values

        result.update({
            'has_data': True,
            'avg_concs': avg_actual_concs,
            'avg_fnorms': avg_fnorms,
            'avg_fnorms_errors': avg_fnorms_errors
        })

        if len(avg_fnorms) >= 3:  # Need at least 3 unique points for a meaningful fit
            try:
                min_val_fnorms_avg = np.min(avg_fnorms)
                max_val_fnorms_avg = np.max(avg_fnorms)

                # Set reasonable initial guesses for bottom and top based on current data range
                p0_bottom_avg = min_val_fnorms_avg
                p0_top_avg = max_val_fnorms_avg

                # Define general bounds for Fnorm, ensuring positive values for Fnorm and Kd
                # Adjusted to be more robust, especially for Fnorm values
                # Allow some wiggle below current min and above current max
                bounds_lower_fnorm_avg = [min(0.0, min_val_fnorms_avg - 0.5), min(0.0, min_val_fnorms_avg - 0.5)]
                bounds_upper_fnorm_avg = [max(3.0, max_val_fnorms_avg + 0.5), max(3.0, max_val_fnorms_avg + 0.5)]

                # Adjust p0 if top < bottom (e.g. inverted curve)
                if p0_top_avg < p0_bottom_avg:
                    p0_top_avg, p0_bottom_avg = p0_bottom_avg, p0_top_avg  # Swap if inverted
                if p0_bottom_avg == p0_top_avg:  # Avoid division by zero if they are the same
                    p0_top_avg += 0.01

                # Clamp initial guesses within some overall broad bounds before specific model bounds
                p0_bottom_avg = np.clip(p0_bottom_avg, bounds_lower_fnorm_avg[0], bounds_upper_fnorm_avg[0])
                p0_top_avg = np.clip(p0_top_avg, bounds_lower_fnorm_avg[1], bounds_upper_fnorm_avg[1])

                if self.model_choice.get() == "Hill":
                    p0_hill_avg = [p0_bottom_avg, p0_top_avg, np.median(avg_log_concs), 1.0]
                    bounds_hill_avg = (
                        [bounds_lower_fnorm_avg[0], bounds_lower_fnorm_avg[1], np.min(avg_log_concs) - 2, 0.1],
                        [bounds_upper_fnorm_avg[0], bounds_upper_fnorm_avg[1], np.max(avg_log_concs) + 2, 5.0]
                    )
                    popt, pcov = curve_fit(binding_curve, avg_log_concs, avg_fnorms, p0=p0_hill_avg,
                                           bounds=bounds_hill_avg, maxfev=10000)
                    result['kd_or_ec50_result'] = 10 ** popt[2]
                    result['fit_success'] = True
                    result['fit_comment'] = "Hill Model"
                    result['fit_popt'] = popt
                elif self.model_choice.get() == "Kd":
                    kd_lower_bound = 1e-12
                    kd_upper_bound = 1e-2
                    median_ligand_conc_avg = np.median(avg_actual_concs)
                    # p0_kd_val_avg is a starting guess for Kd itself, should be within realistic range
                    p0_kd_val_avg = max(1e-9, min(median_ligand_conc_avg, self.c_target, 1e-6))
                    p0_kd_avg = [p0_bottom_avg, p0_top_avg, p0_kd_val_avg]
                    bounds_kd_avg = ([bounds_lower_fnorm_avg[0], bounds_lower_fnorm_avg[1], kd_lower_bound],
                                     [bounds_upper_fnorm_avg[0], bounds_upper_fnorm_avg[1], kd_upper_bound])
                    popt, pcov = curve_fit(
                        lambda c, unbound, bound, kd_fit: kd_model(c, unbound, bound, kd_fit, self.c_target),
                        avg_actual_concs, avg_fnorms, p0=p0_kd_avg, bounds=bounds_kd_avg, maxfev=10000)
                    result['kd_or_ec50_result'] = popt[2]
                    result['fit_success'] = True
                    result['fit_comment'] = "Kd Model"
                    result['fit_popt'] = popt
                elif self.model_choice.get() == "Dose-Response":
                    p0_log_ec50_avg = np.median(avg_log_concs)
                    p0_hill_slope_avg = 1.0
                    p0_dr_avg = [p0_bottom_avg, p0_top_avg, p0_log_ec50_avg, p0_hill_slope_avg]
                    bounds_dr_avg = (
                    [bounds_lower_fnorm_avg[0], bounds_lower_fnorm_avg[1], np.min(avg_log_concs) - 2, -5.0],
                    [bounds_upper_fnorm_avg[0], bounds_upper_fnorm_avg[1], np.max(avg_log_concs) + 2, 5.0])
                    popt, pcov = curve_fit(dose_response_curve, avg_log_concs, avg_fnorms, p0=p0_dr_avg,
                                           bounds=bounds_dr_avg, maxfev=10000)
                    result['kd_or_ec50_result'] = 10 ** popt[2]
                    result['fit_success'] = True
                    result['fit_comment'] = "Dose-Response Model"
                    result['fit_popt'] = popt

            except (RuntimeError, ValueError) as e:
                result['fit_comment'] = f"Fit failed: {e}"
                logging.error(f"Fit failed for {ligand_group_key}: {e}")
            except Exception as e:
                result['fit_comment'] = f"An unexpected error occurred during fit for {ligand_group_key}: {e}"
                logging.exception(f"Unexpected error for {ligand_group_key}")
        else:
            result['fit_comment'] = "Not enough unique data points for fit (need >= 3)."

        return result

    def update_plot(self):
        """Re-draws the main plot with current selected data and fit."""
        self.ax.clear()
        self.ax.set_xscale('log')
        self.ax.set_xlabel('Ligand Concentration [M]')
        self.ax.set_ylabel('Fnorm')
        self.ax.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.7)

        current_ligand_group = self.ligand_group_selector.get()
        if not current_ligand_group or current_ligand_group not in self.all_data_points:
            self.ax.set_title('No Data Selected')
            self.canvas.draw()
            self.results_label.config(text="Results: No data to display.")
            return

        replicates_data = self.all_data_points[current_ligand_group]

        # **MODIFIED**: Get the correct selection dictionary for the current group
        group_selections = self.replicate_selection_vars.get(current_ligand_group, {})

        # Use a distinct color for each replicate for better visibility
        # The 'tab10' colormap has 10 distinct colors
        cmap = plt.cm.get_cmap('tab10', len(replicates_data))

        # Separate and plot individual replicate points first (included and excluded)
        for i, (file_name, data) in enumerate(sorted(replicates_data.items())):
            # **MODIFIED**: Check against the persistent selection state
            if not group_selections.get(file_name, tk.BooleanVar(value=True)).get():
                continue  # Skip if replicate is not selected

            # Separate included and excluded points for plotting
            included_concs = []
            included_fnorms = []
            excluded_concs = []
            excluded_fnorms = []

            for idx, (c, f) in enumerate(zip(data['concentrations'], data['fnorms'])):
                if idx in data['excluded_indices']:
                    excluded_concs.append(c)
                    excluded_fnorms.append(f)
                else:
                    included_concs.append(c)
                    included_fnorms.append(f)

            # Plot included points with normal marker, excluded points with 'x'
            self.ax.scatter(included_concs, included_fnorms, color=cmap(i), s=50, zorder=5, picker=True,
                            # picker=True for interactivity
                            label=f'{os.path.basename(file_name).replace(".csv", "")} (Included)')
            self.ax.scatter(excluded_concs, excluded_fnorms, color=cmap(i), marker='x', s=100, zorder=5, picker=True,
                            label=f'{os.path.basename(file_name).replace(".csv", "")} (Excluded)', alpha=0.6)

        self.ax.set_title(f'Binding Curve: {current_ligand_group}')
        self.ax.legend(fontsize=8, loc='best')

        # --- Get data for fitting using the new helper ---
        # Note: For the main plot, we pass use_selected_replicates=True
        # This ensures only currently checked replicates contribute to the fit.
        processed_data = self._get_processed_group_data(current_ligand_group, use_selected_replicates=True)

        if not processed_data['has_data']:
            self.ax.text(0.5, 0.5, processed_data['fit_comment'],
                         horizontalalignment='center', verticalalignment='center',
                         transform=self.ax.transAxes, color='red', fontsize=12)
            self.results_label.config(text=f"Results: {processed_data['fit_comment']}")
            self.canvas.draw()
            return

        avg_actual_concs = processed_data['avg_concs']
        avg_log_concs = np.log10(avg_actual_concs)  # Recalculate log for plotting consistency
        avg_fnorms = processed_data['avg_fnorms']
        avg_fnorms_errors = processed_data['avg_fnorms_errors']

        # Plot averaged data points with error bars
        self.ax.errorbar(avg_actual_concs, avg_fnorms, yerr=avg_fnorms_errors, fmt='o', color='black', capsize=5,
                         zorder=6,
                         label='Averaged Fnorm (SEM)', markersize=8)

        fit_success = processed_data['fit_success']
        popt = processed_data['fit_popt']
        fit_comment = processed_data['fit_comment']
        kd_or_ec50_result = processed_data['kd_or_ec50_result']

        if fit_success and popt is not None:
            # Generate points for the fit curve
            min_log_c_avg = np.min(avg_log_concs)
            max_log_c_avg = np.max(avg_log_concs)
            # Extend range slightly beyond actual data points for smoother curve
            fine_log_concs_avg = np.linspace(min_log_c_avg - (max_log_c_avg - min_log_c_avg) * 0.1,
                                             max_log_c_avg + (max_log_c_avg - min_log_c_avg) * 0.1, 500)
            fine_concs_avg = 10 ** fine_log_concs_avg

            if self.model_choice.get() == "Hill":
                fit_vals_avg = binding_curve(fine_log_concs_avg, *popt)
                kd_text = self.format_molar_conc(kd_or_ec50_result, "Kd")
                self.ax.plot(fine_concs_avg, fit_vals_avg, 'r--', linewidth=2, label=f'Fit ({kd_text})')
                self.results_label.config(text=f"Results ({fit_comment}): {kd_text}")
            elif self.model_choice.get() == "Kd":
                fit_vals_avg = kd_model(fine_concs_avg, popt[0], popt[1], popt[2], self.c_target)
                kd_text = self.format_molar_conc(kd_or_ec50_result, "Kd")
                self.ax.plot(fine_concs_avg, fit_vals_avg, 'r--', linewidth=2, label=f'Fit ({kd_text})')
                self.results_label.config(text=f"Results ({fit_comment}): {kd_text}")
            elif self.model_choice.get() == "Dose-Response":
                fit_vals_avg = dose_response_curve(fine_log_concs_avg, *popt)
                ec50_text = self.format_molar_conc(kd_or_ec50_result, "EC50")
                self.ax.plot(fine_concs_avg, fit_vals_avg, 'r--', linewidth=2, label=f'Fit ({ec50_text})')
                self.results_label.config(text=f"Results ({fit_comment}): {ec50_text}")
            self.ax.legend(fontsize=8, loc='best')
        else:
            self.ax.text(0.5, 0.5, fit_comment, horizontalalignment='center', verticalalignment='center',
                         transform=self.ax.transAxes, color='red', fontsize=12)
            self.results_label.config(text=f"Results: {fit_comment}")

        self.canvas.draw()

    def generate_overlay_plot(self):
        """Generates a new plot window with overlaid binding curves for two selected ligand groups."""
        group1_key = self.overlay_group1_selector.get()
        group2_key = self.overlay_group2_selector.get()

        if not group1_key or not group2_key:
            messagebox.showwarning("Selection Error", "Please select two ligand groups for overlay.")
            return
        if group1_key == group2_key:
            messagebox.showwarning("Selection Error", "Please select two *different* ligand groups for overlay.")
            return
        if group1_key not in self.all_data_points or group2_key not in self.all_data_points:
            messagebox.showwarning("Data Error",
                                   "Selected ligand groups not found in loaded data. Please reload data or check names.")
            return

        # **MODIFIED**: Process data for both groups, now respecting the replicate selections from the UI.
        data1 = self._get_processed_group_data(group1_key, use_selected_replicates=True)
        data2 = self._get_processed_group_data(group2_key, use_selected_replicates=True)

        # Create new Toplevel window for the overlay plot
        overlay_window = tk.Toplevel(self.master)
        overlay_window.title(f"Overlay: {group1_key} vs {group2_key}")
        overlay_window.geometry("900x700")  # Smaller size for a pop-up plot

        fig_overlay, ax_overlay = plt.subplots(figsize=(7, 5))
        canvas_overlay = FigureCanvasTkAgg(fig_overlay, master=overlay_window)
        canvas_widget_overlay = canvas_overlay.get_tk_widget()
        canvas_widget_overlay.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar_overlay = NavigationToolbar2Tk(canvas_overlay, overlay_window)
        toolbar_overlay.update()
        canvas_widget_overlay.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        ax_overlay.set_xscale('log')
        ax_overlay.set_xlabel('Ligand Concentration [M]')
        ax_overlay.set_ylabel('Fnorm')
        ax_overlay.set_title(f'Binding Curve Overlay\n({self.model_choice.get()} Model)')
        ax_overlay.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.7)

        # Plot Group 1
        if data1['has_data']:
            ax_overlay.errorbar(data1['avg_concs'], data1['avg_fnorms'], yerr=data1['avg_fnorms_errors'],
                                fmt='o', color='blue', capsize=5, zorder=5, label=f'{group1_key} (Avg ± SEM)')
            if data1['fit_success'] and data1['fit_popt'] is not None:
                min_log_c = np.min(np.log10(data1['avg_concs']))
                max_log_c = np.max(np.log10(data1['avg_concs']))
                fine_log_concs = np.linspace(min_log_c - (max_log_c - min_log_c) * 0.1,
                                             max_log_c + (max_log_c - min_log_c) * 0.1, 500)
                fine_concs = 10 ** fine_log_concs

                if self.model_choice.get() == "Hill":
                    fit_vals = binding_curve(fine_log_concs, *data1['fit_popt'])
                    kd_text = self.format_molar_conc(data1['kd_or_ec50_result'], "Kd")
                elif self.model_choice.get() == "Kd":
                    fit_vals = kd_model(fine_concs, data1['fit_popt'][0], data1['fit_popt'][1], data1['fit_popt'][2],
                                        self.c_target)
                    kd_text = self.format_molar_conc(data1['kd_or_ec50_result'], "Kd")
                elif self.model_choice.get() == "Dose-Response":
                    fit_vals = dose_response_curve(fine_log_concs, *data1['fit_popt'])
                    kd_text = self.format_molar_conc(data1['kd_or_ec50_result'], "EC50")

                ax_overlay.plot(fine_concs, fit_vals, 'b-', linewidth=2, label=f'{group1_key} Fit ({kd_text})')
            else:
                # Display fit failure message in plot
                ax_overlay.text(0.05, 0.95, f"{group1_key} Fit: {data1['fit_comment']}",
                                horizontalalignment='left', verticalalignment='top',
                                transform=ax_overlay.transAxes, color='blue', fontsize=8)
        else:
            # Display no data message in plot
            ax_overlay.text(0.05, 0.95, f"{group1_key}: {data1['fit_comment']}",
                            horizontalalignment='left', verticalalignment='top',
                            transform=ax_overlay.transAxes, color='blue', fontsize=8)

        # Plot Group 2
        if data2['has_data']:
            ax_overlay.errorbar(data2['avg_concs'], data2['avg_fnorms'], yerr=data2['avg_fnorms_errors'],
                                fmt='s', color='red', capsize=5, zorder=5, label=f'{group2_key} (Avg ± SEM)')
            if data2['fit_success'] and data2['fit_popt'] is not None:
                min_log_c = np.min(np.log10(data2['avg_concs']))
                max_log_c = np.max(np.log10(data2['avg_concs']))
                fine_log_concs = np.linspace(min_log_c - (max_log_c - min_log_c) * 0.1,
                                             max_log_c + (max_log_c - min_log_c) * 0.1, 500)
                fine_concs = 10 ** fine_log_concs

                if self.model_choice.get() == "Hill":
                    fit_vals = binding_curve(fine_log_concs, *data2['fit_popt'])
                    kd_text = self.format_molar_conc(data2['kd_or_ec50_result'], "Kd")
                elif self.model_choice.get() == "Kd":
                    fit_vals = kd_model(fine_concs, data2['fit_popt'][0], data2['fit_popt'][1], data2['fit_popt'][2],
                                        self.c_target)
                    kd_text = self.format_molar_conc(data2['kd_or_ec50_result'], "Kd")
                elif self.model_choice.get() == "Dose-Response":
                    fit_vals = dose_response_curve(fine_log_concs, *data2['fit_popt'])
                    kd_text = self.format_molar_conc(data2['kd_or_ec50_result'], "EC50")

                ax_overlay.plot(fine_concs, fit_vals, 'r-', linewidth=2, label=f'{group2_key} Fit ({kd_text})')
            else:
                # Display fit failure message in plot
                ax_overlay.text(0.05, 0.90, f"{group2_key} Fit: {data2['fit_comment']}",
                                horizontalalignment='left', verticalalignment='top',
                                transform=ax_overlay.transAxes, color='red', fontsize=8)
        else:
            # Display no data message in plot
            ax_overlay.text(0.05, 0.90, f"{group2_key}: {data2['fit_comment']}",
                            horizontalalignment='left', verticalalignment='top',
                            transform=ax_overlay.transAxes, color='red', fontsize=8)

        ax_overlay.legend(fontsize=8, loc='best')
        fig_overlay.tight_layout()
        canvas_overlay.draw()

        # Add export button to the overlay window
        export_overlay_button = ttk.Button(overlay_window, text="Export Merged Data",
                                           command=lambda: self.export_overlay_plot_data(group1_key, data1, group2_key,
                                                                                         data2))
        export_overlay_button.pack(pady=5)

    def format_molar_conc(self, conc_molar, label):
        if np.isnan(conc_molar):
            return f'{label} = NaN'

        unit = self.display_unit.get()
        if unit == "M":
            return f'{label} = {conc_molar:.2e} M'
        elif unit == "mM":
            return f'{label} = {conc_molar * 1e3:.2f} mM'
        elif unit == "µM":
            return f'{label} = {conc_molar * 1e6:.2f} µM'
        elif unit == "nM":
            return f'{label} = {conc_molar * 1e9:.2f} nM'
        elif unit == "pM":
            return f'{label} = {conc_molar * 1e12:.2f} pM'
        else:  # Fallback to M if unit is invalid
            return f'{label} = {conc_molar:.2e} M'

    def export_main_plot_data(self):
        current_ligand_group = self.ligand_group_selector.get()
        if not current_ligand_group or current_ligand_group not in self.all_data_points:
            messagebox.showwarning("Export Error", "No ligand group selected or data available for export.")
            return

        processed_data = self._get_processed_group_data(current_ligand_group, use_selected_replicates=True)

        if not processed_data['has_data']:
            messagebox.showwarning("Export Error",
                                   f"No valid data to export for '{current_ligand_group}'. {processed_data['fit_comment']}")
            return

        file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx")],
            title=f"Export Data for {current_ligand_group}"
        )
        if not file_path:
            return  # User cancelled

        with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
            # 1. Averaged and Fitted Data
            averaged_df = pd.DataFrame({
                'Concentration (M)': processed_data['avg_concs'],
                'Fnorm_Mean': processed_data['avg_fnorms'],
                'Fnorm_SEM': processed_data['avg_fnorms_errors']
            })
            averaged_df.to_excel(writer, sheet_name='Averaged Data', index=False)

            # 2. Fit Parameters
            fit_params_df = pd.DataFrame()
            if processed_data['fit_success'] and processed_data['fit_popt'] is not None:
                if self.model_choice.get() == "Hill":
                    param_names = ["Bottom", "Top", "logKd", "Hill Slope"]
                elif self.model_choice.get() == "Kd":
                    param_names = ["Unbound (Bottom)", "Bound (Top)", "Kd (M)"]
                    fit_params_df.loc["Target Concentration (M)", "Value"] = self.c_target
                elif self.model_choice.get() == "Dose-Response":
                    param_names = ["Bottom", "Top", "logEC50", "Hill Slope"]

                for i, param_name in enumerate(param_names):
                    fit_params_df.loc[param_name, "Value"] = processed_data['fit_popt'][i]

                # Add Kd/EC50 result
                result_label = "Kd (M)" if self.model_choice.get() in ["Hill", "Kd"] else "EC50 (M)"
                fit_params_df.loc[result_label, "Value"] = processed_data['kd_or_ec50_result']
            else:
                fit_params_df.loc["Fit Status", "Value"] = processed_data['fit_comment']

            fit_params_df.to_excel(writer, sheet_name='Fit Parameters', index=True)

            # 3. Individual Replicate Data (Included & Excluded)
            all_reps_data = []
            for rep_data in processed_data['individual_replicates_data']:
                rep_df = pd.DataFrame({
                    'File Name': rep_data['file_name'],
                    'Concentration (M)': rep_data['Concentration'],
                    'Fnorm': rep_data['Fnorm'],
                    'Status': 'Included'
                })
                # Add excluded points
                if len(rep_data['Excluded_Concentration']) > 0:
                    excluded_df = pd.DataFrame({
                        'File Name': rep_data['file_name'],
                        'Concentration (M)': rep_data['Excluded_Concentration'],
                        'Fnorm': rep_data['Excluded_Fnorm'],
                        'Status': 'Excluded'
                    })
                    rep_df = pd.concat([rep_df, excluded_df], ignore_index=True)
                all_reps_data.append(rep_df)

            if all_reps_data:
                combined_reps_df = pd.concat(all_reps_data, ignore_index=True)
                combined_reps_df.sort_values(by=['File Name', 'Concentration (M)'], inplace=True)
                combined_reps_df.to_excel(writer, sheet_name='Individual Replicates', index=False)

        messagebox.showinfo("Export Successful", f"Data exported to:\n{file_path}")

    def export_overlay_plot_data(self, group1_key, data1, group2_key, data2):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".xlsx",
            filetypes=[("Excel files", "*.xlsx")],
            title=f"Export Overlay Data ({group1_key} vs {group2_key})"
        )
        if not file_path:
            return  # User cancelled

        with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
            # Sheet for Group 1 Averaged Data and Fit
            if data1['has_data']:
                df1 = pd.DataFrame({
                    'Concentration (M)': data1['avg_concs'],
                    'Fnorm_Mean': data1['avg_fnorms'],
                    'Fnorm_SEM': data1['avg_fnorms_errors']
                })
                df1.to_excel(writer, sheet_name=f'{group1_key}_Averaged', index=False)

                fit_params_df1 = pd.DataFrame()
                if data1['fit_success']:
                    if self.model_choice.get() == "Hill":
                        param_names = ["Bottom", "Top", "logKd", "Hill Slope"]
                    elif self.model_choice.get() == "Kd":
                        param_names = ["Unbound", "Bound", "Kd (M)"]
                    elif self.model_choice.get() == "Dose-Response":
                        param_names = ["Bottom", "Top", "logEC50", "Hill Slope"]
                    for i, p_name in enumerate(param_names): fit_params_df1.loc[p_name, "Value"] = data1['fit_popt'][i]
                    result_label = "Kd (M)" if self.model_choice.get() in ["Hill", "Kd"] else "EC50 (M)"
                    fit_params_df1.loc[result_label, "Value"] = data1['kd_or_ec50_result']
                    if self.model_choice.get() == "Kd": fit_params_df1.loc[
                        "Target Concentration (M)", "Value"] = self.c_target
                else:
                    fit_params_df1.loc["Fit Status", "Value"] = data1['fit_comment']
                fit_params_df1.to_excel(writer, sheet_name=f'{group1_key}_FitParams', index=True)

            # Sheet for Group 2 Averaged Data and Fit
            if data2['has_data']:
                df2 = pd.DataFrame({
                    'Concentration (M)': data2['avg_concs'],
                    'Fnorm_Mean': data2['avg_fnorms'],
                    'Fnorm_SEM': data2['avg_fnorms_errors']
                })
                df2.to_excel(writer, sheet_name=f'{group2_key}_Averaged', index=False)

                fit_params_df2 = pd.DataFrame()
                if data2['fit_success']:
                    if self.model_choice.get() == "Hill":
                        param_names = ["Bottom", "Top", "logKd", "Hill Slope"]
                    elif self.model_choice.get() == "Kd":
                        param_names = ["Unbound", "Bound", "Kd (M)"]
                    elif self.model_choice.get() == "Dose-Response":
                        param_names = ["Bottom", "Top", "logEC50", "Hill Slope"]
                    for i, p_name in enumerate(param_names): fit_params_df2.loc[p_name, "Value"] = data2['fit_popt'][i]
                    result_label = "Kd (M)" if self.model_choice.get() in ["Hill", "Kd"] else "EC50 (M)"
                    fit_params_df2.loc[result_label, "Value"] = data2['kd_or_ec50_result']
                    if self.model_choice.get() == "Kd": fit_params_df2.loc[
                        "Target Concentration (M)", "Value"] = self.c_target
                else:
                    fit_params_df2.loc["Fit Status", "Value"] = data2['fit_comment']
                fit_params_df2.to_excel(writer, sheet_name=f'{group2_key}_FitParams', index=True)

            # Combined Averaged Data (for easier comparison)
            if data1['has_data'] and data2['has_data']:
                combined_avg_df = pd.DataFrame({
                    'Concentration (M)': data1['avg_concs'],
                    f'{group1_key}_Fnorm_Mean': data1['avg_fnorms'],
                    f'{group1_key}_Fnorm_SEM': data1['avg_fnorms_errors'],
                    f'{group2_key}_Fnorm_Mean': data2['avg_fnorms'],
                    # Assumes same concentrations, might need merge if not
                    f'{group2_key}_Fnorm_SEM': data2['avg_fnorms_errors']
                })
                # If concentrations are not identical, you'd need to use pd.merge
                # For simplicity, assuming they are the same unique concentrations for now.
                # If not, a more robust merge based on concentration would be needed.
                combined_avg_df.to_excel(writer, sheet_name='Combined Averaged Data', index=False)
            elif data1['has_data']:
                messagebox.showwarning("Export Warning", "Only Group 1 data was available for combined export.")
                df1.to_excel(writer, sheet_name='Combined Averaged Data (Group 1 Only)', index=False)
            elif data2['has_data']:
                messagebox.showwarning("Export Warning", "Only Group 2 data was available for combined export.")
                df2.to_excel(writer, sheet_name='Combined Averaged Data (Group 2 Only)', index=False)
            else:
                messagebox.showwarning("Export Warning", "No data available for combined export.")

        messagebox.showinfo("Export Successful", f"Overlay data exported to:\n{file_path}")


if __name__ == "__main__":
    root = tk.Tk()
    app = MSTAnalyzerGUI(root)
    root.mainloop()