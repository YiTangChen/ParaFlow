#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

find "$SCRIPT_DIR" -name "CMakeCache.txt" -delete
find "$SCRIPT_DIR" -name "CMakeFiles" -type d -exec rm -rf {} +

echo "==> Clean done."