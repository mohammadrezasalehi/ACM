# HCP Denoising Module (Python)

This module performs the final cleaning steps (GSR, Detrending, Filtering) on parcellated HCP data.

## 📦 Installation

This project is designed to run in an isolated environment to avoid conflicts.

### Prerequisites
*   **Python 3.8+** installed and added to System PATH.

### Quick Setup (Windows)
1.  Double-click `setup_env.bat`.
    *   This will create a folder named `env_hcp`.
    *   It will install `numpy`, `scipy`, `nilearn`, etc. automatically.

## 🚀 Usage

1.  Ensure you have processed the MATLAB step first (Outputting .mat files).
2.  Open `hcp_denoising_dual.py` and verify the `PROJECT_ROOT` path in the Config section.
3.  Double-click `run_pipeline.bat`.

## 📁 Output
The script generates two folders in `data/processed_clean/`:
*   `With_GSR/`: For intrinsic topology analysis.
*   `No_GSR/`: For investigating global fluctuations and Simpson's Paradox.