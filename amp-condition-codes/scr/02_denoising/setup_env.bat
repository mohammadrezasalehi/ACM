@echo off
:: ============================================================
:: SETUP SCRIPT for HCP Denoising Module
:: Creates a local Python virtual environment and installs deps.
:: ============================================================

set ENV_NAME=env_hcp

echo [INFO] Checking Python installation...
python --version
if %errorlevel% neq 0 (
    echo [ERROR] Python is not installed or not in PATH.
    echo Please install Python 3.8+ from python.org
    pause
    exit /b
)

echo.
echo [INFO] Creating Virtual Environment: %ENV_NAME% ...
if exist %ENV_NAME% (
    echo [INFO] Environment already exists. Skipping creation.
) else (
    python -m venv %ENV_NAME%
)

echo.
echo [INFO] Activating environment and installing libraries...
call %ENV_NAME%\Scripts\activate

echo [INFO] Upgrading pip...
python -m pip install --upgrade pip

echo [INFO] Installing requirements from requirements.txt...
pip install -r requirements.txt

echo.
echo [SUCCESS] Setup complete!
echo To run the pipeline, use 'run_pipeline.bat'
pause