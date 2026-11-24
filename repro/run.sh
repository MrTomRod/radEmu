#!/usr/bin/env bash
# Run one or all reproduction scripts with plain Rscript (requires installed radEmu and fastEmu).
# Usage:
#   ./run.sh                 # run all *.R scripts in this directory
#   ./run.sh <script.R> ...  # run specific script(s)
set -euo pipefail
DIR="$(cd "$(dirname "$0")" && pwd)"
if [ "$#" -gt 0 ]; then
  scripts=("$@")
else
  scripts=("$DIR"/*.R)
fi
for s in "${scripts[@]}"; do
  [ -f "$s" ] || s="$DIR/$s"
  printf '=== %s ===\n' "$(basename "$s")"
  Rscript "$s"
done
