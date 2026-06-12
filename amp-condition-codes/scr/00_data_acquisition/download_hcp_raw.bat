@echo off
:: =========================================================================
:: SCRIPT: download_hcp_raw.bat
:: AUTHOR: Mohammad Reza Salehi
:: DESC:   Batch downloads MSMAll clean dtseries files from HCP S3 bucket
::         for a list of subjects and all 4 resting-state sessions.
:: REQ:    AWS CLI installed and configured ('aws configure').
:: =========================================================================

setlocal EnableDelayedExpansion

:: --- CONFIGURATION -------------------------------------------------------
:: List of subjects (Ensure this file is next to the script)
set SUB_LIST=hcp_100_unrelated.txt

:: Destination Root Directory (Where to save data)
:: Example: D:\HCP\Raw_Data
set TARGET_ROOT=D:\HCP\Raw_Data

:: List of Sessions to download
:: We iterate over all 4 sessions needed for the thesis
set SESSIONS=rfMRI_REST1_LR rfMRI_REST1_RL rfMRI_REST2_LR rfMRI_REST2_RL
:: -------------------------------------------------------------------------

echo ==========================================================
echo STARTING HCP DATA DOWNLOAD
echo Subject List: %SUB_LIST%
echo Target Dir:   %TARGET_ROOT%
echo ==========================================================

:: Check if subject list exists
if not exist %SUB_LIST% (
    echo [ERROR] Subject list file not found: %SUB_LIST%
    pause
    exit /b
)

:: Loop through Subjects
for /F "tokens=* delims=" %%S in (%SUB_LIST%) do (
    echo.
    echo ----------------------------------------------------------
    echo Processing Subject: %%S
    echo ----------------------------------------------------------

    :: Loop through Sessions
    for %%V in (%SESSIONS%) do (
        
        :: Define local path structure
        set LOCAL_FILE=%TARGET_ROOT%\%%S\%%V\%%V_Atlas_MSMAll_hp2000_clean.dtseries.nii
        
        :: Define S3 path
        set S3_URI=s3://hcp-openaccess/HCP_1200/%%S/MNINonLinear/Results/%%V/%%V_Atlas_MSMAll_hp2000_clean.dtseries.nii

        echo [DOWNLOADING] %%V ...
        
        :: Execute AWS CP command
        :: Note: aws cp automatically creates the file, but we ensure parent dir exists implies organization
        aws s3 cp "!S3_URI!" "!LOCAL_FILE!"
        
        if errorlevel 1 (
            echo [ERROR] Failed to download %%V for subject %%S
        ) else (
            echo [OK] Success.
        )
    )
)

echo.
echo ==========================================================
echo ALL DOWNLOADS COMPLETED.
echo ==========================================================
pause