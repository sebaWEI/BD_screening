# Install bssf dependencies, download pinned databases, then check_requirements.
# Run from anywhere: .\scripts\setup.ps1
$ErrorActionPreference = "Stop"

$Root = Split-Path -Parent $PSScriptRoot
Set-Location $Root

if (-not (Test-Path "pyproject.toml")) {
    Write-Error "Run this script from a bssf checkout (missing pyproject.toml)."
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

$missing = @("blastn", "makeblastdb", "tabix") | Where-Object {
    -not (Get-Command $_ -ErrorAction SilentlyContinue)
}

Write-Host "==> BLAST+, tabix"
if (-not $missing) {
    Write-Host "Already on PATH."
}
elseif (Get-Command conda -ErrorAction SilentlyContinue) {
    conda install -y -c bioconda blast htslib
}
else {
    Write-Error @"
conda is not on PATH, so BLAST+ and tabix were not installed.

Install Miniforge (https://github.com/conda-forge/miniforge#miniforge3), open a conda PowerShell, then:
  conda install -y -c bioconda blast htslib

Re-run .\scripts\setup.ps1 with that environment activated.
"@
}

Write-Host "==> Databases (GENCODE ~1 GB; dbSNP tabix index ~3 MB) and BLAST index"
if (-not $env:https_proxy -and -not $env:HTTPS_PROXY) {
    try {
        $tcp = New-Object System.Net.Sockets.TcpClient
        $iar = $tcp.BeginConnect("127.0.0.1", 7897, $null, $null)
        if ($iar.AsyncWaitHandle.WaitOne(500) -and $tcp.Connected) {
            $env:http_proxy = "http://127.0.0.1:7897"
            $env:https_proxy = "http://127.0.0.1:7897"
            $env:HTTP_PROXY = $env:http_proxy
            $env:HTTPS_PROXY = $env:https_proxy
            Write-Host "    Using local proxy http://127.0.0.1:7897"
        }
        $tcp.Close()
    } catch {}
}
uv run bssf db init --dbsnp-common-all --gencode-v45-transcripts

Write-Host "==> check_requirements"
uv run bssf check_requirements
if ($LASTEXITCODE -ne 0) {
    exit $LASTEXITCODE
}
Write-Host "Requirements passed. Continue with the Usage section of README.md."
