#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
WARMUPS="${JKLU_BENCH_WARMUPS:-100}"
ITERATIONS="${JKLU_BENCH_ITERATIONS:-500}"
REPEATS="${JKLU_BENCH_REPEATS:-5}"
NATIVE_BIN="${NATIVE_KLU_BENCH:-/private/tmp/native_klu_z_bench}"
CP="$ROOT_DIR/target/classes:$ROOT_DIR/target/test-classes:$(cat "$ROOT_DIR/target/benchmark-classpath.txt")"

if [ ! -x "$NATIVE_BIN" ]; then
  echo "native benchmark not executable: $NATIVE_BIN" >&2
  exit 1
fi

echo "case,status,targetPct,weightedPct,repeat,javaFactorMs,javaRefactorMs,nativeFactorMs,nativeRefactorMs,javaSolveMs,nativeSolveMs,luEntries,relativeResidual"
tmp_status="$(mktemp)"
trap 'rm -f "$tmp_status"' EXIT

run_case() {
  local name="$1"
  local matrix="$2"
  local rhs="$3"
  local target_pct="$4"

  if [ ! -f "$matrix" ]; then
    echo "$name,MISSING_MATRIX,$target_pct,,,,,,,,,,"
    echo "$name" >> "$tmp_status"
    return
  fi
  if [ "$rhs" != "-" ] && [ ! -f "$rhs" ]; then
    echo "$name,MISSING_RHS,$target_pct,,,,,,,,,,"
    echo "$name" >> "$tmp_status"
    return
  fi

  local tmp_case
  tmp_case="$(mktemp)"
  for repeat in $(seq 1 "$REPEATS"); do
    local java_out native_out
    if ! java_out="$(java -Xmx24g -Djklu.benchmark.reuseSymbolic=true \
      -cp "$CP" edu.ufl.cise.klu.bench.ZkluMatrixMarketBenchmark \
      "$matrix" "$rhs" "$WARMUPS" "$ITERATIONS" 2>&1)"; then
      echo "$name,FAIL_JKLU,$target_pct,,$repeat,,,,,,,,"
      echo "$java_out" >&2
      echo "$name" >> "$tmp_status"
      rm -f "$tmp_case"
      return
    fi

    if ! native_out="$(NATIVE_KLU_REUSE_SYMBOLIC=1 "$NATIVE_BIN" \
      "$matrix" "$rhs" "$WARMUPS" "$ITERATIONS" 2>&1)"; then
      echo "$name,FAIL_NATIVE,$target_pct,,$repeat,,,,,,,,"
      echo "$native_out" >&2
      echo "$name" >> "$tmp_status"
      rm -f "$tmp_case"
      return
    fi

    python3 - "$name" "$target_pct" "$repeat" "$java_out" "$native_out" <<'PY' >> "$tmp_case"
import sys

def parse(line):
    data = {}
    for part in line.strip().split(","):
        if "=" in part:
            key, value = part.split("=", 1)
            data[key] = value
    return data

name, target_text, repeat, java_line, native_line = sys.argv[1:6]
target = float(target_text)
j = parse(java_line)
n = parse(native_line)
weighted_java = 0.2 * float(j["factorMs"]) + 0.8 * float(j["refactorMs"])
weighted_native = 0.2 * float(n["factorMs"]) + 0.8 * float(n["refactorMs"])
weighted_pct = 100.0 * weighted_native / weighted_java
status = "PASS" if weighted_pct >= target else "MISS"
print(",".join([
    name,
    status,
    f"{target:.1f}",
    f"{weighted_pct:.2f}",
    repeat,
    j["factorMs"],
    j["refactorMs"],
    n["factorMs"],
    n["refactorMs"],
    j["solveMs"],
    n["solveMs"],
    j.get("luEntries", ""),
    j.get("relativeResidual", ""),
]))
PY
  done
  local row
  row="$(python3 - "$tmp_case" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as handle:
    rows = list(csv.reader(handle))
best = max(rows, key=lambda row: float(row[3]))
print(",".join(best))
PY
)"
  rm -f "$tmp_case"
  echo "$row"
  case "$row" in
    *,PASS,*) ;;
    *) echo "$name" >> "$tmp_status" ;;
  esac
}

run_case Ckt24 "$ROOT_DIR/target/ipss-matrices/Ckt24-ymatrix.mtx" "$ROOT_DIR/target/ipss-matrices/Ckt24-rhs.mtx" 70
run_case IEEE8500 "$ROOT_DIR/target/ipss-matrices/IEEE8500-ymatrix.mtx" "$ROOT_DIR/target/ipss-matrices/IEEE8500-rhs.mtx" 70
run_case ACTIVSg25k "$ROOT_DIR/target/ipss-matrices/ACTIVSg25k-ymatrix.mtx" "$ROOT_DIR/target/ipss-matrices/ACTIVSg25k-rhs.mtx" 80
run_case ACTIVSg70K "$ROOT_DIR/target/tamu-power-matrices/ACTIVSg70K/ACTIVSg70K.mtx" - 80
run_case OpenEI "$ROOT_DIR/target/ipss-matrices/OpenEI-ymatrix.mtx" "$ROOT_DIR/target/ipss-matrices/OpenEI-rhs.mtx" 80

if [ -s "$tmp_status" ]; then
  exit 1
fi
