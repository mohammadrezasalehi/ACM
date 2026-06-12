# HCP Data Acquisition

This module handles the automated download of resting-state fMRI data from the Human Connectome Project (HCP) S3 buckets.

## 📋 Overview
We retrieve the `MSMAll` registered, `ICA-FIX` denoised dense time-series for 100 unrelated subjects. The data is hosted on Amazon S3 (`s3://hcp-openaccess`).

## 🛠 Prerequisites

1.  **AWS CLI**: Must be installed on the system.
    *   Download: [AWS CLI Installer](https://aws.amazon.com/cli/)
2.  **HCP Credentials**:
    *   You must have an account at [db.humanconnectome.org](https://db.humanconnectome.org/).
    *   You need to accept the data usage terms for the HCP 1200 release.
    *   Obtain your **Access Key ID** and **Secret Access Key** from the HCP database profile or AWS IAM (if using your own setup, though HCP usually provides specific instructions for S3 access).
3.  **Configuration**:
    *   Run `aws configure` in your terminal and enter your credentials.

## 📂 File Description

*   `download_hcp_raw.bat`: The main script. It loops through the subject list and downloads all 4 sessions (`REST1_LR`, `REST1_RL`, `REST2_LR`, `REST2_RL`).
*   `hcp_100_unrelated.txt`: A text file containing the Subject IDs (one per line).

## ⚙️ Usage

1.  Open `download_hcp_raw.bat`.
2.  Update the `TARGET_ROOT` variable to point to your hard drive (e.g., `D:\HCP\Data`).
3.  Ensure `hcp_100_unrelated.txt` is in the same folder.
4.  Run the batch file.

## 📦 Data Specifications
*   **Source Bucket:** `s3://hcp-openaccess/HCP_1200/`
*   **File Path:** `<SubjectID>/MNINonLinear/Results/<Session>/<Session>_Atlas_MSMAll_hp2000_clean.dtseries.nii`
*   **Size:** Approx. 0.5 GB per file (Total ~200 GB for 100 subjects × 4 sessions).