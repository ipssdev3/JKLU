#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
MATRIX_ROOT="${1:-$ROOT_DIR/target/tamu-power-matrices}"
MANIFEST="${2:-$ROOT_DIR/scripts/tamu-power-matrices.txt}"
WARMUPS="${JKLU_BENCH_WARMUPS:-3}"
ITERATIONS="${JKLU_BENCH_ITERATIONS:-10}"
NATIVE_BIN="${NATIVE_KLU_BENCH:-/private/tmp/native_klu_z_bench}"
CP="$ROOT_DIR/target/classes:$ROOT_DIR/target/test-classes:$(cat "$ROOT_DIR/target/benchmark-classpath.txt")"

tmp_csv="$(mktemp)"
trap 'rm -f "$tmp_csv"' EXIT
echo "case,status,javaFactorMs,nativeFactorMs,factorPct,javaRefactorMs,nativeRefactorMs,refactorPct,weightedJavaMs,weightedNativeMs,weightedPct,javaSolveMs,nativeSolveMs,solvePct" | tee "$tmp_csv"

while IFS= read -r name; do
  case "$name" in
    ""|\#*) continue ;;
  esac
  group="${name%%/*}"
  matrix="${name##*/}"
  matrix_file="$MATRIX_ROOT/$matrix/${matrix}.mtx"
  if [ ! -f "$matrix_file" ]; then
    matrix_file="$MATRIX_ROOT/$group/$matrix/${matrix}.mtx"
  fi
  if [ ! -f "$matrix_file" ]; then
    echo "missing matrix: $matrix_file" >&2
    exit 1
  fi

  if ! java_out="$(java -Xmx24g -Djklu.benchmark.reuseSymbolic=true \
    -cp "$CP" edu.ufl.cise.klu.bench.ZkluMatrixMarketBenchmark \
    "$matrix_file" - "$WARMUPS" "$ITERATIONS" 2>&1)"; then
    echo "$name,FAIL_JKLU,,,,,,,,," | tee -a "$tmp_csv"
    echo "$java_out" >&2
    continue
  fi
  if ! native_out="$(NATIVE_KLU_REUSE_SYMBOLIC=1 "$NATIVE_BIN" \
    "$matrix_file" - "$WARMUPS" "$ITERATIONS" 2>&1)"; then
    echo "$name,FAIL_NATIVE,,,,,,,,," | tee -a "$tmp_csv"
    echo "$native_out" >&2
    continue
  fi

  python3 - "$name" "$java_out" "$native_out" <<'PY' | tee -a "$tmp_csv"
import sys

def parse(line):
    out = {}
    for part in line.strip().split(","):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out

name, java_line, native_line = sys.argv[1:4]
j = parse(java_line)
n = parse(native_line)
def pct(key):
    return 100.0 * float(n[key]) / float(j[key])
weighted_java = 0.2 * float(j["factorMs"]) + 0.8 * float(j["refactorMs"])
weighted_native = 0.2 * float(n["factorMs"]) + 0.8 * float(n["refactorMs"])
print(",".join([
    name, "OK",
    j["factorMs"], n["factorMs"], f"{pct('factorMs'):.2f}",
    j["refactorMs"], n["refactorMs"], f"{pct('refactorMs'):.2f}",
    f"{weighted_java:.6f}", f"{weighted_native:.6f}",
    f"{100.0 * weighted_native / weighted_java:.2f}",
    j["solveMs"], n["solveMs"], f"{pct('solveMs'):.2f}",
]))
PY
done < "$MANIFEST"

python3 - "$tmp_csv" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as handle:
    rows = [row for row in csv.DictReader(handle) if row["status"] == "OK"]
if not rows:
    raise SystemExit

def avg(name):
    return sum(float(row[name]) for row in rows) / len(rows)

print(",".join([
    "AVERAGE", str(len(rows)),
    f"{avg('javaFactorMs'):.6f}", f"{avg('nativeFactorMs'):.6f}", f"{avg('factorPct'):.2f}",
    f"{avg('javaRefactorMs'):.6f}", f"{avg('nativeRefactorMs'):.6f}", f"{avg('refactorPct'):.2f}",
    f"{avg('weightedJavaMs'):.6f}", f"{avg('weightedNativeMs'):.6f}", f"{avg('weightedPct'):.2f}",
    f"{avg('javaSolveMs'):.6f}", f"{avg('nativeSolveMs'):.6f}", f"{avg('solvePct'):.2f}",
]))
PY
