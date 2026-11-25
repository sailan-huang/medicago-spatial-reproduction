param(
    [switch]$SkipDownload,
    [switch]$SkipNotebook
)

$ErrorActionPreference = "Stop"
$PSNativeCommandUseErrorActionPreference = $true
$projectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $projectRoot

if (-not (Test-Path ".venv\Scripts\python.exe")) {
    uv venv .venv --python 3.12
}

uv sync --frozen
& ".venv\Scripts\python.exe" -m ipykernel install --prefix .venv --name medicago-cds --display-name "Python (Medicago CDS)" | Out-Null

if (-not $SkipDownload) {
    & ".venv\Scripts\python.exe" "scripts\download_data.py"
}

& ".venv\Scripts\python.exe" -m pytest -q

if (-not $SkipNotebook) {
    & ".venv\Scripts\jupytext.exe" --to notebook --execute --set-kernel medicago-cds --output "notebooks\medicago_spatial_reproduction.ipynb" "notebooks\medicago_spatial_reproduction.py"
}

& ".venv\Scripts\python.exe" "scripts\build_report.py"
& ".venv\Scripts\python.exe" "scripts\validate_results.py"

