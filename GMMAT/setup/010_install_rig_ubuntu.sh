#!/usr/bin/env bash
set -euo pipefail

if [[ "${EUID}" -eq 0 ]]; then
  echo "Run this script as a normal user with sudo available, not as root." >&2
  exit 1
fi

if [[ ! -f /etc/os-release ]]; then
  echo "Cannot detect OS: /etc/os-release not found." >&2
  exit 1
fi

. /etc/os-release

if [[ "${ID:-}" != "ubuntu" && "${ID_LIKE:-}" != *debian* ]]; then
  echo "This script targets Ubuntu/Debian systems with apt." >&2
  exit 1
fi

if ! command -v curl >/dev/null 2>&1; then
  echo "curl is required but not found." >&2
  exit 1
fi

if ! command -v sudo >/dev/null 2>&1; then
  echo "sudo is required but not found." >&2
  exit 1
fi

KEY_PATH="/etc/apt/trusted.gpg.d/rig.gpg"
LIST_PATH="/etc/apt/sources.list.d/rig.list"
REPO_LINE='deb http://rig.r-pkg.org/deb rig main'

if [[ ! -f "$KEY_PATH" ]]; then
  echo "Adding rig apt key..."
  sudo curl -L https://rig.r-pkg.org/deb/rig.gpg -o "$KEY_PATH"
else
  echo "rig apt key already present: $KEY_PATH"
fi

if [[ ! -f "$LIST_PATH" ]] || ! grep -Fqx "$REPO_LINE" "$LIST_PATH" 2>/dev/null; then
  echo "Adding rig apt repository..."
  echo "$REPO_LINE" | sudo tee "$LIST_PATH" >/dev/null
else
  echo "rig apt repository already present: $LIST_PATH"
fi

echo "Installing rig..."
sudo apt-get update
sudo apt-get install -y r-rig

echo "rig installation complete."
rig --version
rig system detect-platform || true
