@echo off
:: ============================================================
:: RUN SCRIPT
:: Activates the environment and runs the python pipeline.
:: ============================================================

set ENV_NAME=env_hcp
set SCRIPT_NAME=hcp_denoising_dual.py

if not exist %ENV_NAME%\Scripts\activate.bat (
    echo [ERROR] Virtual environment not found.
    echo Please run 'setup_env.bat' first.
    pause
    exit /b
)

echo [INFO] Activating environment...
call %ENV_NAME%\Scripts\activate

echo [INFO] Running Pipeline: %SCRIPT_NAME%
python %SCRIPT_NAME%

echo.
echo [INFO] Pipeline finished.
pause