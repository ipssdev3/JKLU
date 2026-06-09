#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="${1:-$ROOT_DIR/scripts/tamu-power-matrices.txt}"
OUT_DIR="${2:-$ROOT_DIR/target/tamu-power-matrices}"

mkdir -p "$OUT_DIR"

while IFS= read -r name; do
  case "$name" in
    ""|\#*) continue ;;
  esac
  group="${name%%/*}"
  matrix="${name##*/}"
  matrix_dir="$OUT_DIR/$matrix"
  matrix_file="$matrix_dir/${matrix}.mtx"
  if [ -f "$matrix_file" ]; then
    echo "present $name $matrix_file"
    continue
  fi
  mkdir -p "$matrix_dir"
  archive="$OUT_DIR/${group}_${matrix}.tar.gz"
  url="https://sparse.tamu.edu/MM/$group/$matrix.tar.gz"
  echo "download $url"
  curl -L --fail --retry 3 --output "$archive" "$url"
  tar -xzf "$archive" -C "$OUT_DIR"
  if [ ! -f "$matrix_file" ]; then
    found="$(find "$OUT_DIR" -path "*/${matrix}/${matrix}.mtx" -print -quit)"
    if [ -z "$found" ]; then
      echo "missing extracted matrix for $name" >&2
      exit 1
    fi
  fi
done < "$MANIFEST"
