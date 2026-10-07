#!/usr/bin/env bash
# Install bssf dependencies, download pinned databases, then check_requirements.
# Run from anywhere: ./scripts/setup.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if [[ ! -f pyproject.toml ]]; then
  echo "Run this script from a bssf checkout (missing pyproject.toml)." >&2
  exit 1
fi

echo "==> uv"
if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
  export PATH="${HOME}/.local/bin:${PATH}"
fi
if ! command -v uv >/dev/null 2>&1; then
  echo "uv was installed but is not on PATH. Open a new terminal and re-run ./scripts/setup.sh." >&2
  exit 1
fi

echo "==> Python project"
uv python install 3.12
uv sync

need_native=0
for cmd in blastn makeblastdb tabix; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    need_native=1
  fi
done

echo "==> BLAST+, tabix"
if [[ "$need_native" -eq 0 ]]; then
  echo "Already on PATH."
elif command -v brew >/dev/null 2>&1; then
  brew install blast htslib
elif command -v conda >/dev/null 2>&1; then
  conda install -y -c bioconda blast htslib
else
  cat >&2 <<'EOF'
Neither Homebrew nor conda is on PATH, so BLAST+ and tabix were not installed.

macOS (Homebrew):
  brew install blast htslib

Linux / macOS without Homebrew:
  conda install -y -c bioconda blast htslib

Then re-run ./scripts/setup.sh
EOF
  exit 1
fi

echo "==> Databases (GENCODE ~1 GB; dbSNP tabix index ~3 MB) and BLAST index"
# Proxy helps NCBI HTTPS tabix / index fetch from CN networks.
if [[ -z "${https_proxy:-}${HTTPS_PROXY:-}" ]] && nc -z -G 1 127.0.0.1 7897 2>/dev/null; then
  export http_proxy=http://127.0.0.1:7897
  export https_proxy=http://127.0.0.1:7897
  export HTTP_PROXY="$http_proxy"
  export HTTPS_PROXY="$https_proxy"
  echo "    Using local proxy http://127.0.0.1:7897"
fi
uv run bssf db init --dbsnp-common-all --gencode-v45-transcripts

echo "==> check_requirements"
uv run bssf check_requirements
echo "Requirements passed. Continue with the Usage section of README.md."
