#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SELECTOR="${R_SELECTOR:-release}"
PROJECT_LIB_ROOT="${PROJECT_LIB_ROOT:-$ROOT/r-lib}"

if ! command -v rig >/dev/null 2>&1; then
  echo "rig not found. Run ./010_install_rig_ubuntu.sh first." >&2
  exit 1
fi

echo "Adding R via rig: ${R_SELECTOR}"
rig add "$R_SELECTOR"

# Avoid rig resolve here: it may require a live network lookup.
rig default "$R_SELECTOR"
rig system make-links || true
rig system setup-user-lib || true

RESOLVED_R="$(python3 - <<'PY'
import subprocess
lines = subprocess.check_output(['rig', 'list'], text=True).splitlines()
for raw in lines:
    line = raw.strip()
    if not line or line.startswith('-'):
        continue
    parts = line.split()
    if len(parts) >= 2 and parts[0] == '*' and parts[1] != 'name':
        print(parts[1])
        break
PY
)"

if [[ -z "$RESOLVED_R" ]]; then
  echo "Failed to determine the default installed R version from 'rig list'." >&2
  exit 1
fi

echo "Resolved/default R version: $RESOLVED_R"

PROJECT_LIB="${PROJECT_LIB_ROOT}/${RESOLVED_R}"
mkdir -p "$PROJECT_LIB"

echo "Installing R packages into: $PROJECT_LIB"
R_LIBS_USER="$PROJECT_LIB" rig run -r "$RESOLVED_R" -f "$ROOT/install_r_packages.R" -- "$PROJECT_LIB"

echo
echo "Setup complete."
echo "R version: $RESOLVED_R"
echo "Project library: $PROJECT_LIB"
echo "Use this R for the current shell with one of:"
echo "  rig default $RESOLVED_R"
echo "  rig run -r $RESOLVED_R -f your_script.R"
echo "For this project, export:"
echo "  export R_LIBS_USER=$PROJECT_LIB"
