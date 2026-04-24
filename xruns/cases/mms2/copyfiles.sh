#!/usr/bin/env bash
set -euo pipefail

FILES_TO_COPY=(
  "config.txt"
  "*.py"
  "*.sh"
  "*.geo"
  "*.msh"
  "*.petsc"
)

usage() {
  echo "Usage: $0 <dest_dir>"
  exit 1
}

warn() {
  if [[ -t 2 ]]; then
    printf '\033[1;91mWARNING:\033[0m %s\n' "$*" >&2
  else
    printf 'Warning: %s\n' "$*" >&2
  fi
}

[[ $# -eq 1 ]] || usage
dest="$1"

src_dir="$(pwd -P)"

# ---- check "same as current directory" BEFORE mkdir ----
case "$dest" in
  "."|"./"|"$src_dir")
    echo "Destination is the current directory. Nothing to do."
    exit 0
    ;;
esac

if command -v realpath >/dev/null 2>&1; then
  dest_abs="$(realpath -m -- "$dest")"
  if [[ "$dest_abs" == "$src_dir" ]]; then
    echo "Destination is the current directory. Nothing to do."
    exit 0
  fi
fi
# --------------------------------------------------------

mkdir -p -- "$dest"
dest_dir="$(cd "$dest" && pwd -P)"

# Expand patterns; WARN (don't error) if any pattern matches nothing
shopt -s nullglob
expanded=()
for pat in "${FILES_TO_COPY[@]}"; do
  matches=( $pat )
  if [[ ${#matches[@]} -eq 0 ]]; then
    warn "pattern matched nothing: $pat"
    continue
  fi
  expanded+=( "${matches[@]}" )
done
shopt -u nullglob

# Print expanded array
echo "Expanded files (${#expanded[@]}):"
for f in "${expanded[@]}"; do
  printf '  %q\n' "$f"
done

if [[ ${#expanded[@]} -eq 0 ]]; then
  echo "No files matched. Nothing to copy."
  exit 0
fi

# Copy (files only)
for f in "${expanded[@]}"; do
  [[ -f "$f" ]] || continue
  cp -p -- "$f" "$dest_dir/"
  echo "Copied: $f -> $dest_dir/"
done
