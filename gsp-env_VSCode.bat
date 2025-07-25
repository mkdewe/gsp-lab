@echo off
setlocal enabledelayedexpansion

:: === Conda detection ===
set "CONDA_ROOT="

:: Try conda info --base first
for /f "delims=" %%i in ('conda info --base 2^>nul') do set "CONDA_ROOT=%%i"

:: If not found, search PATH
if not defined CONDA_ROOT (
    for /f "delims=" %%i in ('where conda 2^>nul') do (
        set "CONDA_PATH=%%~dpi"
        if exist "!CONDA_PATH!\..\..\condabin\conda.bat" set "CONDA_ROOT=!CONDA_PATH!\..\.."
        if exist "!CONDA_PATH!\..\..\Scripts\conda.exe" set "CONDA_ROOT=!CONDA_PATH!\..\.."
    )
)

:: If still not found, check common locations
if not defined CONDA_ROOT (
    if exist "C:\ProgramData\Miniconda3\" set "CONDA_ROOT=C:\ProgramData\Miniconda3"
    if exist "C:\ProgramData\Anaconda3\" set "CONDA_ROOT=C:\ProgramData\Anaconda3"
    if exist "!USERPROFILE!\Miniconda3\" set "CONDA_ROOT=!USERPROFILE!\Miniconda3"
    if exist "!USERPROFILE!\Anaconda3\" set "CONDA_ROOT=!USERPROFILE!\Anaconda3"
)

:: Error if conda not found
if not defined CONDA_ROOT (
    echo [ERROR] Conda installation not found
    pause
    exit /b
)

:: Get absolute path
for %%I in ("!CONDA_ROOT!") do set "CONDA_ROOT=%%~fI"

:: === Find activate.bat ===
set "ACTIVATE_PATH="
if exist "!CONDA_ROOT!\condabin\activate.bat" set "ACTIVATE_PATH=!CONDA_ROOT!\condabin\activate.bat"
if exist "!CONDA_ROOT!\Scripts\activate.bat" set "ACTIVATE_PATH=!CONDA_ROOT!\Scripts\activate.bat"
if exist "!CONDA_ROOT!\Library\bin\activate.bat" set "ACTIVATE_PATH=!CONDA_ROOT!\Library\bin\activate.bat"

if not defined ACTIVATE_PATH (
    echo [ERROR] activate.bat not found in Conda installation
    pause
    exit /b
)

:: === Set environment variables ===
set "ENV_NAME=gsp-env"
set "PROJECT_DIR=%~dp0"
set "VSCODE_SETTINGS=%PROJECT_DIR%.vscode\settings.json"

:: === Locate VS Code ===
set "VSCODE_EXE="
where code >nul 2>nul
if %ERRORLEVEL%==0 (
    set "VSCODE_EXE=code"
) else (
    if exist "%USERPROFILE%\AppData\Local\Programs\Microsoft VS Code\Code.exe" (
        set "VSCODE_EXE=""%USERPROFILE%\AppData\Local\Programs\Microsoft VS Code\Code.exe"""
    )
)

if not defined VSCODE_EXE (
    echo [ERROR] Visual Studio Code not found
    pause
    exit /b
)

:: === Create .vscode folder ===
if not exist "%PROJECT_DIR%.vscode" mkdir "%PROJECT_DIR%.vscode"

:: === Generate VS Code settings ===
set "CMD=!ACTIVATE_PATH! && conda activate !ENV_NAME!"
set "CMD=!CMD:\=/!"

(
echo {
echo   "terminal.integrated.profiles.windows": {
echo     "Conda (auto-activate)": {
echo       "path": "cmd.exe",
echo       "args": [
echo         "/k",
echo         "!CMD!"
echo       ]
echo     }
echo   },
echo   "terminal.integrated.defaultProfile.windows": "Conda (auto-activate)"
echo }
) > "%VSCODE_SETTINGS%"

echo [SUCCESS] VS Code configured with Conda auto-activation

:: === Launch VS Code ===
start "" /min cmd /c "!VSCODE_EXE! "%PROJECT_DIR%" & exit"
exit
