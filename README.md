MST Ultimate GUI Data Analyzer

A user-friendly graphical interface (GUI) application for analyzing MicroScale Thermophoresis (MST) spectral shift data. This tool streamlines the processing of raw .xlsx MST files, performs binding curve fitting (Hill and Quadratic models), generates interactive plots, and provides comprehensive data export in Excel format.

✨ Features

    Easy Data Import: Load raw MST .xlsx files (from NanoTemper Monolith instruments).

    Automated Processing: Automatically extracts and organizes concentrations, Fnorm, Spectral Shift (670/650 ratio), and initial fluorescence from your raw data.

    Replicate Handling: Seamlessly averages data from multiple replicates for robust analysis, including calculation of standard error of the mean (SEM).

    Flexible Binding Models:

        Hill Binding Model: Standard model for calculating binding affinity (K_d).

        Quadratic Binding Model: Accounts for ligand depletion, providing more accurate K_d values, especially crucial for high-affinity interactions or when target concentration is close to or higher than the K_d.

    Interactive Visualizations: Generates dynamic plots of Fnorm vs. Concentration and Spectral Shift vs. Concentration, complete with fitted curves, residual plots, and adjustable axes.

    Data Overlay: Compare and visualize binding curves from two distinct data groups (e.g., different protein variants, buffer conditions) on the same plot.

    Comprehensive Excel Export: All processed data, averaged results, fit parameters, and combined overlay data are meticulously organized and exported into a single multi-sheet Excel workbook.

    Intuitive GUI: Built with Tkinter for a straightforward user experience.

🚀 Getting Started

Prerequisites

    Python 3.8+ (The script uses f-strings and other modern Python features.)

    Required Python packages:

        pandas

        numpy

        matplotlib

        scipy

        xlsxwriter (for Excel export)

        openpyxl (often used by pandas for reading/writing Excel, good to have)

You can install these packages using pip:
Bash

pip install pandas numpy matplotlib scipy xlsxwriter openpyxl

Installation & Usage (From Source)

    Clone the repository:
    Bash

git clone https://github.com/YourUsername/YourRepoName.git
cd YourRepoName

Run the application:
Bash

    python MST_ultimateGUI_NT_V5.py

Installation & Usage (Executable)

For users who prefer not to install Python or dependencies, a standalone executable for Windows is available.

    Download the latest release:
    Visit the Releases page and download the MST_ultimateGUI_NT_V5.exe file (or the .zip containing it).

    Run the executable:
    Double-click the MST_ultimateGUI_NT_V5.exe file to launch the application.

📈 How to Use the GUI

    Launch Application: Run the Python script or the executable.

    Load Data:

        Click "Load Data Group 1" to select one or more raw MST .xlsx files for your first experimental group.

        Click "Load Data Group 2" to select files for a second group (optional, for comparison).

        The application will display the number of files loaded for each group.

    Set Parameters:

        Target Concentration (M): Enter the concentration of your labeled target molecule in Molar (e.g., 5e-9 for 5 nM).

        Initial Kd Guess (M): Provide an initial estimate for the dissociation constant in Molar. This helps the fitting algorithm.

        Initial Hill Slope Guess: Enter an initial guess for the Hill slope (typically around 1.0).

    Process & Plot:

        Click the "Process & Plot Data" button.

        The application will process the data, perform curve fitting, and display two plots: Fnorm vs. Concentration and Spectral Shift vs. Concentration.

        Fitting parameters (Kd, Hill Slope, etc.) and their standard errors will be displayed.

        Residual plots are also generated to help assess fit quality.

    Export Data:

        Click the "Export All Data" button.

        A "Save As" dialog will appear, allowing you to choose the name and location for your comprehensive Excel output file. This file will contain multiple sheets with raw data, averaged data, fit results, and combined overlay data.

🛠️ Development

Building the Executable (using PyInstaller)

If you modify the source code and want to create your own executable:

    Install PyInstaller:
    Bash

pip install pyinstaller

Build the executable:
Navigate to the directory containing MST_ultimateGUI_NT_V5.py in your terminal and run:
Bash

pyinstaller --onefile --windowed --icon=MST.ico --hidden-import=xlsxwriter --hidden-import=openpyxl MST_ultimateGUI_NT_V5.py

    --onefile: Creates a single executable file.

    --windowed or -w: Prevents a console window from appearing when the GUI starts.

    --icon=MST.ico: Sets the icon of the executable. (Ensure MST.ico is in the same directory as your script).

    --hidden-import=xlsxwriter --hidden-import=openpyxl: Explicitly includes these crucial libraries for Excel functionality, preventing common packaging issues.

The executable will be generated in the dist folder.

Okay, let's craft a Release Note and a README file for your MST data analysis tool, "MST_ultimateGUI_NT_V5.py," based on the code snippets and functionality observed.

Release Note Example

This would typically be used on GitHub Releases.

Release v1.0.0 (or appropriate version) - MST Ultimate GUI Data Analyzer

This is the initial public release of the MST Ultimate GUI Data Analyzer, a powerful and user-friendly tool for processing, analyzing, and visualizing MicroScale Thermophoresis (MST) spectral shift data. This application provides robust features for handling multiple MST experiments, calculating key parameters, and generating high-quality plots and reports.

Key Features:

    Load and Process Raw MST Data: Easily import raw .xlsx files from NanoTemper Monolith instruments.

    Automatic Data Extraction: Extracts crucial data points including concentrations, Fnorm, spectral shift (670/650 ratio), and initial fluorescence.

    Replicate Averaging: Automatically calculates average Fnorm and spectral shift values, along with standard errors, for replicate experiments.

    Binding Curve Fitting:

        Utilizes the Hill binding model for direct binding affinity (Kd) determination.

        Includes a Quadratic binding model (Kd model) to account for target depletion, providing more accurate Kd values, especially for high-affinity interactions or when target concentration is significant relative to Kd.

    Interactive Plotting: Generates clear and customizable plots of Fnorm vs. Concentration and Spectral Shift vs. Concentration with fitted curves.

    Overlay Functionality: Compare up to two different data groups (e.g., wild-type vs. mutant, or different conditions) on the same plot for easy comparison and export.

    Comprehensive Data Export:

        Export raw processed data for each group.

        Export averaged data tables.

        Export fit parameters (Kd, Hill Slope, etc.).

        Export combined overlay data (averaged Fnorm, spectral shift, and SEM for both groups).

        All exports are conveniently generated in a single Excel workbook with multiple sheets.

    User-Friendly Tkinter GUI: Intuitive interface for easy navigation and operation.

How to Use:

    Launch the application: Run the MST_ultimateGUI_NT_V5.exe file.

    Load Data: Use the "Load Data Group 1" and "Load Data Group 2" buttons to select your raw MST .xlsx files. You can load single or multiple files for each group.

    Set Parameters: Input your target concentration and initial guess for Kd and Hill Slope.

    Process & Plot: Click "Process & Plot Data" to generate analyses and plots.

    Export: Use the "Export All Data" button to save your results to a comprehensive Excel file.

Fixes & Improvements in this Release:

    Initial release, establishing core functionality.

Known Issues / Limitations:

    Assumes consistent concentration ranges for combined overlay plots and data export.

    (Add any other known limitations or bugs here, e.g., "Error handling for malformed input files could be improved.")

Feedback & Support:

For any issues, suggestions, or feedback, please open an issue on the [suspicious link removed].

README.md Example

This would be your main README.md file in the root of your GitHub repository.

MST Ultimate GUI Data Analyzer

A user-friendly graphical interface (GUI) application for analyzing MicroScale Thermophoresis (MST) spectral shift data. This tool streamlines the processing of raw .xlsx MST files, performs binding curve fitting (Hill and Quadratic models), generates interactive plots, and provides comprehensive data export in Excel format.

✨ Features

    Easy Data Import: Load raw MST .xlsx files (from NanoTemper Monolith instruments).

    Automated Processing: Automatically extracts and organizes concentrations, Fnorm, Spectral Shift (670/650 ratio), and initial fluorescence from your raw data.

    Replicate Handling: Seamlessly averages data from multiple replicates for robust analysis, including calculation of standard error of the mean (SEM).

    Flexible Binding Models:

        Hill Binding Model: Standard model for calculating binding affinity (K_d).

        Quadratic Binding Model: Accounts for ligand depletion, providing more accurate K_d values, especially crucial for high-affinity interactions or when target concentration is close to or higher than the K_d.

    Interactive Visualizations: Generates dynamic plots of Fnorm vs. Concentration and Spectral Shift vs. Concentration, complete with fitted curves, residual plots, and adjustable axes.

    Data Overlay: Compare and visualize binding curves from two distinct data groups (e.g., different protein variants, buffer conditions) on the same plot.

    Comprehensive Excel Export: All processed data, averaged results, fit parameters, and combined overlay data are meticulously organized and exported into a single multi-sheet Excel workbook.

    Intuitive GUI: Built with Tkinter for a straightforward user experience.

🚀 Getting Started

Prerequisites

    Python 3.8+ (The script uses f-strings and other modern Python features.)

    Required Python packages:

        pandas

        numpy

        matplotlib

        scipy

        xlsxwriter (for Excel export)

        openpyxl (often used by pandas for reading/writing Excel, good to have)

You can install these packages using pip:
Bash

pip install pandas numpy matplotlib scipy xlsxwriter openpyxl

Installation & Usage (From Source)

    Clone the repository:
    Bash

git clone https://github.com/YourUsername/YourRepoName.git
cd YourRepoName

Run the application:
Bash

    python MST_ultimateGUI_NT_V5.py

Installation & Usage (Executable)

For users who prefer not to install Python or dependencies, a standalone executable for Windows is available.

    Download the latest release:
    Visit the Releases page and download the MST_ultimateGUI_NT_V5.exe file (or the .zip containing it).

    Run the executable:
    Double-click the MST_ultimateGUI_NT_V5.exe file to launch the application.

📈 How to Use the GUI

    Launch Application: Run the Python script or the executable.

    Load Data:

        Click "Load Data Group 1" to select one or more raw MST .xlsx files for your first experimental group.

        Click "Load Data Group 2" to select files for a second group (optional, for comparison).

        The application will display the number of files loaded for each group.

    Set Parameters:

        Target Concentration (M): Enter the concentration of your labeled target molecule in Molar (e.g., 5e-9 for 5 nM).

        Initial Kd Guess (M): Provide an initial estimate for the dissociation constant in Molar. This helps the fitting algorithm.

        Initial Hill Slope Guess: Enter an initial guess for the Hill slope (typically around 1.0).

    Process & Plot:

        Click the "Process & Plot Data" button.

        The application will process the data, perform curve fitting, and display two plots: Fnorm vs. Concentration and Spectral Shift vs. Concentration.

        Fitting parameters (Kd, Hill Slope, etc.) and their standard errors will be displayed.

        Residual plots are also generated to help assess fit quality.

    Export Data:

        Click the "Export All Data" button.

        A "Save As" dialog will appear, allowing you to choose the name and location for your comprehensive Excel output file. This file will contain multiple sheets with raw data, averaged data, fit results, and combined overlay data.

🛠️ Development

Building the Executable (using PyInstaller)

If you modify the source code and want to create your own executable:

    Install PyInstaller:
    Bash

pip install pyinstaller

Build the executable:
Navigate to the directory containing MST_ultimateGUI_NT_V5.py in your terminal and run:
Bash

    pyinstaller --onefile --windowed --icon=MST.ico --hidden-import=xlsxwriter --hidden-import=openpyxl MST_ultimateGUI_NT_V5.py

        --onefile: Creates a single executable file.

        --windowed or -w: Prevents a console window from appearing when the GUI starts.

        --icon=MST.ico: Sets the icon of the executable. (Ensure MST.ico is in the same directory as your script).

        --hidden-import=xlsxwriter --hidden-import=openpyxl: Explicitly includes these crucial libraries for Excel functionality, preventing common packaging issues.

    The executable will be generated in the dist folder.

🤝 Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the issues page.

📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

🙏 Acknowledgments

    NanoTemper Technologies for developing MicroScale Thermophoresis technology.

    The developers of pandas, numpy, matplotlib, and scipy for their invaluable scientific computing libraries.

    The PyInstaller community for enabling easy distribution of Python applications.
