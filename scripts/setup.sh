#!/usr/bin/env bash
# Install bsst dependencies, download pinned databases, then check_requirements.
# Run from anywhere: ./scripts/setup.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if [[ ! -f pyproject.toml ]]; then
  echo "Run this script from a bsst checkout (missing pyproject.toml)." >&2
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

echo "==> Databases (~2 GB) and BLAST index"
uv run bsst db init --dbsnp-common-all --gencode-v45-transcripts

echo "==> check_requirements"
uv run bsst check_requirements
echo "Requirements passed. Continue with the Usage section of README.md."
