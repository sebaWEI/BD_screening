# Install bsst dependencies, download pinned databases, then check_requirements.
# Run from anywhere: .\scripts\setup.ps1
$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
Set-Location $Root

if (-not (Test-Path "pyproject.toml")) {
    Write-Error "Run this script from a bsst checkout (missing pyproject.toml)."
}

function Refresh-Path {
    $machine = [Environment]::GetEnvironmentVariable("Path", "Machine")
    $user = [Environment]::GetEnvironmentVariable("Path", "User")
    $env:Path = "$user;$machine;$env:USERPROFILE\.local\bin"
}

Write-Host "==> uv"
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    irm https://astral.sh/uv/install.ps1 | iex
    Refresh-Path
}
if (-not (Get-Command uv -ErrorAction SilentlyContinue)) {
    Write-Error "uv was installed but is not on PATH. Open a new PowerShell and re-run .\scripts\setup.ps1."
}

Write-Host "==> Python project"
uv python install 3.12
uv sync

$missing = @("RNAup", "blastn", "makeblastdb", "tabix") | Where-Object {
    -not (Get-Command $_ -ErrorAction SilentlyContinue)
}

Write-Host "==> RNAup, BLAST+, tabix"
if (-not $missing) {
    Write-Host "Already on PATH."
}
elseif (Get-Command conda -ErrorAction SilentlyContinue) {
    conda install -y -c conda-forge -c bioconda viennarna blast htslib
}
else {
    Write-Error @"
conda is not on PATH, so RNAup, BLAST+, and tabix were not installed.

Install Miniforge (https://github.com/conda-forge/miniforge#miniforge3), open a conda PowerShell, then:
  conda install -y -c conda-forge -c bioconda viennarna blast htslib

Re-run .\scripts\setup.ps1 with that environment activated.
"@
}

Write-Host "==> Databases (~2 GB) and BLAST index"
uv run bsst db init --dbsnp-common-all --gencode-v45-transcripts

Write-Host "==> check_requirements"
uv run bsst check_requirements
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}
Write-Host "Requirements passed. Continue with the Usage section of README.md."
