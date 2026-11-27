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
if (-not (Test-Path "node_modules\docx")) {
    npm install --silent --omit=dev
}
& ".venv\Scripts\python.exe" -m ipykernel install --prefix .venv --name medicago-cds --display-name "Python (Medicago CDS)" | Out-Null

if (-not $SkipDownload) {
    & ".venv\Scripts\python.exe" "scripts\download_data.py"
}

& ".venv\Scripts\python.exe" -m pytest -q

if (-not $SkipNotebook) {
    & ".venv\Scripts\jupytext.exe" --to notebook --execute --set-kernel medicago-cds --output "notebooks\medicago_spatial_reproduction.ipynb" "notebooks\medicago_spatial_reproduction.py"
}

& ".venv\Scripts\python.exe" "scripts\build_report.py"
& ".venv\Scripts\python.exe" "scripts\render_pdf_previews.py"
$env:NODE_PATH = (Join-Path $projectRoot "node_modules")
& node "scripts\build_docx_report.js"
& ".venv\Scripts\python.exe" "scripts\validate_docx_report.py"
& ".venv\Scripts\python.exe" "scripts\validate_results.py"

