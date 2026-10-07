#!/usr/bin/env bash
# Fetch the dbSNP b157 GRCh38.p14 tabix index (~3 MB) and point config at the
# HTTPS GCF URL. --variants then queries regions on demand (no 28 GB download).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

if [[ -z "${https_proxy:-}${HTTPS_PROXY:-}" ]] && nc -z -G 1 127.0.0.1 7897 2>/dev/null; then
  export http_proxy=http://127.0.0.1:7897
  export https_proxy=http://127.0.0.1:7897
  export HTTP_PROXY="$http_proxy"
  export HTTPS_PROXY="$https_proxy"
  echo "using proxy $https_proxy"
fi

.venv/bin/bssf db init --dbsnp-common-all
.venv/bin/bssf check_requirements
