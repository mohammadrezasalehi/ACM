# HCP Time-Series Parcellation Module

This module is responsible for reducing the dimensionality of high-resolution HCP dense time-series (`.dtseries.nii`) into parcel-based time-series using the **Connectome Workbench**.

## 📂 File Structure

*   `main_parcellation_pipeline.m`: **(Run this)** The master script. It handles configuration, session management (REST1/REST2), and batch processing loops.
*   `parcellate_cifti_wrapper.m`: **(Do not run)** A helper function that interfaces with `wb_command` to perform the actual parcellation and conversion to `.mat`.

## ⚙️ Prerequisites

1.  **MATLAB** (Tested on R2018b+).
2.  **Connectome Workbench**: Must be installed. Update `wb_cmd_path` in the main script.
3.  **HCP Data**: Cleaned surface data (`*Atlas_MSMAll_hp2000_clean.dtseries.nii`).
4.  **Atlas**: A `.dlabel.nii` file (e.g., Schaefer 200).

## 🚀 How to Run

1.  Open `main_parcellation_pipeline.m`.
2.  **Edit Global Configuration**:
    *   Set `project_root` to your project folder.
    *   Set `wb_cmd_path` to your local `wb_command.exe`.
3.  **Select Session**:
    *   Change the variable `CURRENT_SESSION` to one of the predefined options:
        ```matlab
        CURRENT_SESSION = 'REST1_LR'; % Options: 'REST1_LR', 'REST1_RL', 'REST2_LR', 'REST2_RL'
        ```
    *   *Note:* The script automatically maps the correct raw data paths based on this variable.
4.  **Run the script**.

## 📤 Outputs

The script automatically creates a structured output directory based on the Atlas name and Session:

```text
data/inter_parcellated/
└── Schaefer200_Kong17/         <-- Atlas Name
    ├── REST1_LR/               <-- Session Name
    │   ├── 100307.mat          <-- Result (N_Parcels x Time)
    │   ├── 100408.mat
    │   └── ...
    └── REST1_RL/
        └── 
```

## 🛠 Troubleshooting
*   **wb_command failed:** Ensure the path to `wb_command.exe` is correct and does not contain spaces/symbols that might confuse the shell, or that the input `.dtseries` file exists.
*   **Dimension Mismatch:** Ensure your Atlas (`.dlabel.nii`) matches the space of the data (usually fs_LR 32k or 91k grayordinates).
