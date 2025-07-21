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
