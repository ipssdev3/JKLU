# JKLU 2.0.0 Direct Complex Kernel Feature Plan

## Release Summary

JKLU `2.0.0` is planned as the direct Java complex numeric kernel release. It
builds on the `1.1.0` complex API and facade work, but replaces the internal
complex numeric backend with a pure-Java direct complex sparse LU kernel.

The `2.0.0` version number is intentional because this work changes solver
internals that determine pivoting, singular detection, diagnostics, memory use,
and performance. Public APIs should remain source-compatible where possible,
but downstream users should receive a major-version signal for the numeric
backend change.

Branch:

```text
feature/direct-complex-kernel
```

No native SuiteSparse, KLUSolve, KLUSolveX, or Eigen dependency will be added.
Complex values continue to use interleaved Java `double[]` storage:

```text
[re0, im0, re1, im1, ...]
```

## Scope

- Keep public APIs unchanged: callers continue using `klu_z_factor`,
  `klu_z_refactor`, `klu_z_solve`, `klu_z_tsolve`, diagnostics, extraction, and
  `KluComplexSparseSet`.
- Keep Java 8 compatibility.
- Remove the block-real backend from `tdcomplex`.
- Reuse existing `KLU_symbolic` analysis, ordering, and BTF data.
- Match SuiteSparse-style `klu_z_*` behavior in pure Java where practical.
- Validate against CSparseJ, JKLU real behavior where applicable, and native
  KLUSolve/KLUSolveX on the same large matrices where native libraries are
  available.
- Use native KLUSolve/KLUSolveX performance on the same large power-system
  matrix, same RHS count, same ordering mode where possible, and same cold/warm
  workflow as the primary performance target.

## Completed Direct-Kernel Work

| Status | Feature | Notes |
| --- | --- | --- |
| Done | Block-real backend removal | Removed `realSymbolic`/`realNumeric` backend fields from `KLU_z_numeric`. |
| Done | Direct complex numeric state | Added direct interleaved complex state for `klu_z_factor`, `klu_z_refactor`, `klu_z_solve`, and `klu_z_tsolve`. |
| Done | Sparse complex row LU storage | Replaced initial dense direct complex core with sparse row storage so memory scales with stored entries and fill, not `n*n`. |
| Done | Active pivot-column incidence | Pivot selection and elimination traverse rows that currently contain the pivot column instead of scanning every trailing row. |
| Done | Direct complex diagnostics/extraction path | Updated `klu_z_extract`, `klu_z_sort`, `klu_z_rgrowth`, `klu_z_condest`, `klu_z_rcond`, and `klu_z_flops` to use direct complex numeric state. |
| Done | Direct-kernel test coverage | Added tests for transpose solve, row scaling, extraction, diagnostics, same-pattern refactor, sparse-set updates, and CSparseJ comparison. |
| Done | Dklu-style same-pattern refactor clearing | Complex `klu_z_refactor` now clears known U/diagonal/L pattern rows like real `Dklu_refactor`, avoiding dynamic touched-set bookkeeping in the numeric update loop. |
| Done | OpenEI-scale flat direct storage default | Raised `jklu.complex.factor.flatColumnMaxN` default to `100000` after benchmarks showed flat solve/refactor is a net win on OpenEI-scale sparse matrices. |
| Done | Flat-only direct L/U cache for large matrices | OpenEI-scale flat-cache path no longer duplicates per-column array-of-arrays storage by default; `jklu.complex.factor.keepColumnArrays=true` retains the old duplicate cache for debugging. |
| Done | Split active row heap | Initial factorization now keeps only rows before the current pivot in the ordered heap and appends current/future active rows for pivot/scale scans, reducing heap maintenance while preserving the active set. |
| Done | Factor scale-loop cleanup | Initial factorization now computes the complex pivot denominator once per column, matching the refactor loop shape. |
| Done | Singleton BTF block fast path | Direct complex factor/refactor now handles one-column BTF blocks without entering the general kernel, matching the SuiteSparse `factor2` structure. |
| Done | Direct BTF `F` refactor coverage | Added a forced-direct test proving off-block `Offp`/`Offi`/`Offx` values are preserved, solved, extracted, and updated by same-pattern refactor for identity-pivot BTF cases. |
| Done | Factor generation-mark workspace | Initial factorization now uses per-column integer generation marks for active/touched rows, reducing boolean-array reset and bookkeeping in the dynamic factor path without changing refactor semantics. |
| Done | Heap-only prior-row active list | Initial factorization now keeps only prior rows in the ordered heap and scans the touched pattern for pivot/L construction, removing duplicate future-row active-list storage. |
| Done | Large-matrix reach-kernel auto path | Matrices at or above `jklu.complex.factor.reachKernelMinN` now use the DFS reach kernel by default, keep U unsorted during factor updates, and keep cached U unsorted by default because flat refactor/solve do not require per-column U sorting. |
| Done | Large-matrix packed direct input default | Matrices that use the large reach path now use packed direct-column input by default, reducing large-case input/refill overhead while preserving the existing opt-out property. |
| Done | Reach-kernel heap bookkeeping removal | The DFS reach path no longer feeds prior rows into the heap-only `ActiveList`, which the reach path never polls. |
| Done | Same-object refactor pattern shortcut | Complex numeric state now keeps the caller's original `Ap/Ai` references alongside owned copies so same-object refactors can skip the full pattern scan. |

## Feature Tracking

| Priority | Status | Feature | Target |
| --- | --- | --- | --- |
| P0 | In progress | Direct Java complex numeric kernel | `2.0.0` |
| P0 | Done | Apply `KLU_symbolic` BTF/AMD ordering in complex factorization | `2.0.0`; `P/Q` permutation is applied in factor and solve |
| P0 | In progress | True independent `Symbolic.R` BTF block factorization | `2.0.0`; current implementation observes `R` boundaries, but still does not truly factor each BTF block independently the way real `Dklu_factor` does; block-local traversal, singleton handling, pivot metadata, and numeric storage remain release-blocking |
| P0 | In progress | SuiteSparse/JKLU-style off-diagonal `F` entries | `2.0.0`; initial `Offp`/`Offi`/`Offx` storage, solve, extract, and same-pattern refactor tests pass for identity-pivot direct BTF cases, but full SuiteSparse/JKLU parity still requires final `Offi = Pinv[Offi]` behavior for off-diagonal pivot cases plus diagnostics coverage |
| P0 | Done | Sparse LU entry traversal for complex solve/extract/diagnostics | `2.0.0`; solve, extraction, and diagnostics upper-norm traversal done |
| P0 | Done | Sparse L/U extraction sized by actual `lnz`/`unz` | `2.0.0` |
| P1 | Done | Public API validation hardening for malformed CSC/RHS inputs | `2.0.0`; factor, solve, transpose solve, extraction, and rgrowth covered |
| P1 | In progress | True complex numeric reuse optimization for `klu_z_refactor` | `2.0.0`; direct CSC-to-column setup is added, but numeric storage/reach reuse remains |
| P1 | In progress | Gilbert-Peierls column kernel with primitive workspace and pruning | `2.0.0`; large matrices use a left-looking column-oriented path, but the implementation still has row-based `TreeMap` fallback and does not yet match the real `Dklu_kernel` symbolic reach, column traversal, and pruning |
| P1 | Planned | Transpose solve optimization or validated cache | `2.0.0` |
| P1 | Planned | Performance profiling loop driven by native KLUSolve parity on large power-system matrices | `2.0.0` |
| P2 | In progress | Permanent benchmark suite | `2.0.0`; generated banded, Matrix Market replay, and InterPSS Y-matrix export harnesses added |
| P2 | Planned | Matrix Market import/export round-trip tests if import is added | Follow-up |

## Review-Driven Improvement Plan

Code review of the `feature/direct-complex-kernel` branch found that the direct
complex kernel is functionally covered by small parity tests, but still needs
large sparse-performance and public-API robustness work before a `2.0.0`
release.

| Priority | Status | Finding | Improvement | Verification |
| --- | --- | --- | --- | --- |
| P0 | In progress | Complex numeric factorization accepts `KLU_symbolic`; `P/Q` are applied and the current sparse-row kernel now observes `R` BTF boundaries with initial off-diagonal `F` storage, but it still does not truly factor each BTF block independently like `Dklu_factor`, does not yet fully own `F` in SuiteSparse/JKLU style, and is not yet the real KLU Gilbert-Peierls column kernel. | First make the `R`/`F` behavior match real `Dklu_factor`, then replace the remaining row-based `TreeMap` elimination/fallback with symbolic reach, pruning, and primitive complex workspace. | Passing focused BTF/off-diagonal `F` solve/extract test and Matrix Market replay; add non-natural `P/Q`, singleton block, off-block `F`, and pivot-metadata parity checks. |
| P0 | Done | Complex solve stores sparse LU but scans dense row ranges during forward/back substitution. | Done for `klu_z_solve`: forward/back substitution now walks stored row entries. Diagnostics upper-norm traversal now walks stored U entries. | Passing `mvn test`; generated sparse timing benchmark smoke run passes. |
| P0 | Done | `klu_z_extract` emits dense triangular L/U patterns and can overflow arrays sized from `lnz`/`unz`. | `klu_z_extract` now emits actual stored sparse L/U entries and validates caller array capacity before writing. | Passing tests cover diagonal sparse extraction with arrays sized from `Numeric.lnz` and `Numeric.unz`, plus undersized extraction rejection. |
| P1 | Done | Public complex entry points can throw runtime array exceptions for malformed `Ap`, `Ai`, `Ax`, `B_offset`, or RHS capacity. | Hardened `klu_z_factor`, `klu_z_solve`, `klu_z_tsolve`, `klu_z_extract`, and `klu_z_rgrowth` matrix validation. | Passing tests cover malformed CSC input and undersized RHS/extraction buffers with `KLU_INVALID`/`FALSE`. |
| P1 | In progress | Large matrices use a column-oriented direct path, but remaining unsupported cases still fall back to row-based `TreeMap<Integer, Complex>` elimination; the direct path also lacks full SuiteSparse pruning, pivot logging parity, and reusable numeric storage. | Continue porting the real `Dklu_kernel` shape to complex: non-recursive symbolic reach, sparse `X` workspace, partial pivoting with diagonal preference, Eisenstat-Liu pruning, and same-pattern numeric reuse. | Passing `mvn test`; current paired OpenEI default-BTF fill matches native exactly at `960,910` LU entries. The large-matrix default now uses reach traversal with unsorted U during factor and keeps cached U unsorted because flat refactor/solve do not need sorted U. Latest clean OpenEI reuse-symbolic paired ratios are about `75.4%` factor, `87.0%` refactor, and `87.4%` solve after adding the packed no-off-diagonal direct input path; factor still misses the 85% target. |
| P1 | Done | Promote measured large-matrix reach traversal where it helps. | Auto-enable reach traversal only above the measured large-matrix threshold while keeping smaller cases on the heap path. | Passing `mvn test`; OpenEI default reuse-symbolic factor improves to `36.237375` ms with exact native fill and residual, while BUS11856 remains on the heap path. |
| P1 | Planned | `klu_z_tsolve` builds and factors a transpose numeric object for each call. | Solve the transpose from existing LU factors when possible, or cache a validated transpose factorization if exact KLU semantics require it. | Add repeated transpose-solve benchmark and parity tests for scaled and pivoted matrices. |
| P2 | Planned | Complex diagnostics are approximate and use dense traversal in places. | Rework diagnostics to use stored sparse LU entries and document any statistic that remains an estimate. | Add diagnostics regression tests for diagonal, triangular, and pivoted systems. |

Implementation order:

1. Public API hardening and extraction correctness.
2. `P/Q` symbolic ordering integration in `klu_z_factor` and `klu_z_solve`.
3. True `Symbolic.R` BTF block factorization parity.
4. SuiteSparse/JKLU-style off-diagonal `F` storage, solve, extraction, and diagnostics.
5. Sparse traversal for solve, extraction, and diagnostics.
6. Large power-system profiling loop and benchmark harness.
7. Gilbert-Peierls column kernel and primitive workspace optimization.
8. Transpose-solve optimization and benchmark suite.

## Direct Kernel Implementation Details

- Reuse `KLU_symbolic.P`, `KLU_symbolic.Q`, and `KLU_symbolic.R` instead of
  factoring the whole matrix in natural order.
- Done: apply `P/Q` permutations when building the complex numeric matrix and
  when loading/storing RHS values in `klu_z_solve`.
- In progress: factor each `R` BTF block independently. Current code observes
  `Symbolic.R`, but it still does not truly factor each block independently the
  way `Dklu_factor` does. The TODO is to make traversal, singleton handling,
  pivot metadata, and block-local numeric storage match `Dklu_factor`.
- In progress: preserve off-diagonal `F` entries in `Offp`/`Offi`/`Offx`.
  Current focused solve/extract tests pass, but release acceptance requires the
  same `F` ownership, solve ordering, extraction, diagnostics, and refactor
  behavior as SuiteSparse/JKLU.
- In progress: replace row-based `TreeMap` elimination with a complex
  Gilbert-Peierls column kernel using symbolic reach and pruning that follows
  `Dklu_kernel` closely enough for large sparse performance parity. The current
  implementation is not yet the real KLU symbolic reach/pruning kernel.
- Use complex magnitude for pivot selection.
- Use stable complex division for pivots.
- Use complex multiply-subtract for numeric updates.
- Keep row scaling real-valued by dividing complex matrix entries by `Rs[row]`.
- Store LU entries sparsely and walk stored entries in solve/extract/diagnostic
  paths.
- Keep same-pattern refactor source-compatible while improving numeric reuse.

## KLU User Guide and Source Audit Findings

The KLU User Guide describes KLU as a BTF-preordered, fill-reduced,
left-looking sparse LU solver for circuit-style matrices. The SuiteSparse
`klu_factor.c` and `klu_kernel.c` source adds several implementation details
that the direct complex path must keep aligning with:

- `factor2` factors each `Symbolic.R` block independently and handles
  singleton blocks without entering the general kernel. The complex direct path
  observes block boundaries, but singleton-block fast-path and per-block
  numeric metadata parity remain open.
- KLU composes the block-local pivot permutation `Pblock` with symbolic `P`
  into `Pnum`, then recomputes `Pinv`. The complex direct path currently uses
  identity internal pivots on the fast direct path, so off-diagonal pivot
  metadata parity remains a release gate.
- KLU stores off-diagonal `F` entries during factorization and finally rewrites
  `Offi[p] = Pinv[Offi[p]]`. The complex path has initial `Offp`/`Offi`/`Offx`
  support, but final `F` row-permutation parity and refactor diagnostics still
  need focused tests.
- KLU allocates numeric storage by BTF block (`LUbx[block]`, `LUsize[block]`)
  and uses symbolic fill estimates (`Lnz[block]`) to size each block. JKLU's
  complex direct state is currently global-column oriented; per-block storage
  and fill metadata remain future parity work.
- KLU's kernel uses non-recursive DFS symbolic reach plus Eisenstat-Liu pruning
  before numeric updates. JKLU's complex path now uses a split active row heap,
  but it is still not a full symbolic DFS/pruning port.

## 2.0 Kernel TODO

These items are release-blocking for the direct Java complex kernel unless they
are explicitly deferred with benchmark evidence.

| Priority | Status | TODO | Acceptance |
| --- | --- | --- | --- |
| P0 | In progress | Make complex factorization process each `Symbolic.R` BTF block independently in the same style as real `Dklu_factor`; current code only observes the block boundaries and is not yet equivalent. | Block-local factorization, pivot order, `Pnum`, `Lip`, `Uip`, singleton handling, and per-block numeric allocation are verified against focused BTF fixtures and Matrix Market replay. |
| P0 | In progress | Complete off-diagonal `F` handling in SuiteSparse/JKLU style; current `F` support is partial and not yet release-equivalent. | `Offp`, `Offi`, and interleaved `Offx` own only off-block entries; solve, transpose solve where applicable, extract, diagnostics, and refactor preserve `F` semantics. |
| P0 | In progress | Replace remaining row-based `TreeMap` elimination/fallback on large sparse cases with a real KLU-like Gilbert-Peierls column traversal. | Main large-matrix path uses non-recursive symbolic reach, sparse `X` workspace, active-column traversal, and Eisenstat-Liu pruning modeled on `Dklu_kernel`; fallback is limited to small or unsupported cases and is measured separately. |
| P1 | In progress | Match real-kernel pruning and traversal decisions closely enough for OpenEI/BUS large matrix parity work. | Profiling shows candidate pivots, row/column updates, fill, and refactor timings moving toward native KLUSolve/KLUSolveX, with no residual regression. |
| P1 | In progress | Port KLU source audit gaps into direct complex factorization. | Singleton BTF blocks and identity-pivot direct `F` refactor coverage are done; per-block numeric metadata, final `F` row permutation through `Pinv`, and `Pblock`/local pivot composition remain open against real `Dklu_factor` behavior. |
| P1 | Planned | Add parity tests for non-natural `P/Q` plus multi-block `R` cases with nonempty `F`. | Tests fail if the complex path accidentally factors the full matrix as one block or drops off-block entries. |
| P1 | In progress | Keep diagonal-only direct kernel from mishandling matrices that need off-diagonal numeric pivoting. | Direct kernel now checks for structural diagonal candidates before running; focused test verifies an off-diagonal-structured case still solves when the direct threshold is forced low. Full `Dklu_factor` pivot metadata parity remains open. |

## Performance Profiling Loop

Performance work must be driven by measured large sparse matrices, not only by
micro-optimizations. The goal is to move JKLU direct complex performance closer
to native KLUSolve/KLUSolveX for the same matrix and solve workflow. Each
optimization cycle should follow this loop:

1. Capture representative matrices.
2. Establish JKLU and native KLUSolve/KLUSolveX baselines where native libraries
   are available.
3. Profile the current bottleneck.
4. Implement one focused improvement.
5. Re-run the same cases.
6. Keep the change only if correctness remains stable and performance improves.

Representative matrix sources:

- Generated sparse diagonal, banded, block-diagonal, and unsymmetric test
  matrices for quick repeatable profiling inside JKLU.
- Captured complex CSC matrices from `/Users/ipssdev/github/ipss-plugin`
  transmission/network tests.
- Captured complex CSC matrices from `/Users/ipssdev/github/ipss-plugin`
  three-phase/OpenDSS feeder tests.
- Downloaded SuiteSparse/TAMU Matrix Market power-system and power-network
  matrices from `https://sparse.tamu.edu/`, using a broad average so OpenEI
  does not become the only optimization target.
- Optional Matrix Market exports for cases that can be shared or archived
  without downstream project dependencies.

Baseline comparison targets:

- JKLU `1.1.0` complex backend where available.
- JKLU `2.0.0` current direct complex kernel.
- CSparseJ complex LU for pure-Java comparison.
- Native KLUSolve/KLUSolveX on the same matrix, RHS count, and workflow where
  local native libraries are available.

Optional native parity helper:

```bash
c++ -O3 -std=c++17 \
  -I/Users/ipssdev/github/klusolve/build/SuiteSparse/KLU/Include \
  -I/Users/ipssdev/github/klusolve/build/SuiteSparse/AMD/Include \
  -I/Users/ipssdev/github/klusolve/build/SuiteSparse/BTF/Include \
  -I/Users/ipssdev/github/klusolve/build/SuiteSparse/COLAMD/Include \
  -I/Users/ipssdev/github/klusolve/build/SuiteSparse/SuiteSparse_config \
  -L/Users/ipssdev/github/klusolve/build \
  -Wl,-rpath,/Users/ipssdev/github/klusolve/build \
  src/test/native/native_klu_z_bench.cpp -lklusolvex \
  -o /private/tmp/native_klu_z_bench

/private/tmp/native_klu_z_bench \
  target/ipss-matrices/OpenEI-ymatrix.mtx \
  target/ipss-matrices/OpenEI-rhs.mtx 10 50
```

Native parity measurements must separate:

- Cold workflow: analyze/order, factor, and solve from a fresh matrix.
- Warm workflow: same-pattern numeric refactor plus solve.
- Solve-only workflow: repeated RHS solves after factorization.
- Matrix assembly/export overhead: measured separately and excluded from core
  solver parity ratios.

Metrics to collect for every benchmark case:

- Matrix dimensions: `n`, `nnz(A)`, nonzero distribution, block count, max block
  size, and ordering mode.
- Ordering time, factor time, refactor time, solve time, transpose-solve time,
  and total solve workflow time.
- Refactor-heavy numeric-kernel score:
  `weightedMs = 0.2 * factorMs + 0.8 * refactorMs`; parity is
  `100 * weightedNativeMs / weightedJavaMs`. Solve time remains a correctness
  and regression gate, not the primary optimization score, because repeated
  power-system workflows spend more useful time in same-pattern refactor.
- `lnz`, `unz`, off-diagonal entry count, fill ratio, and sparse LU stored
  entries.
- Heap use, allocation rate, and garbage collection count/time.
- Residual norm and downstream power-flow result deltas.
- Singular/rank metadata and pivot count changes.
- JKLU/native timing ratio for the refactor-heavy numeric-kernel score,
  cold workflow, warm refactor workflow, and solve-only workflow.

Profiling tools and modes:

- Start with JVM wall-clock timing around analysis, factor, refactor, and solve.
- Use Java Flight Recorder or async-profiler for large InterPSS matrix cases
  when available.
- Use allocation profiling to identify object churn from sparse row maps,
  `Complex` allocation, boxing, and temporary arrays.
- Run warm refactor/solve loops separately from cold analyze/factor loops so
  AMD/BTFJ costs are not confused with numeric-kernel costs.
- Use the `jklu.benchmark.reuseSymbolic=true` and
  `NATIVE_KLU_REUSE_SYMBOLIC=1` benchmark modes to compare warm numeric
  factor/refactor/solve with the same symbolic analysis reused across
  iterations. This separates Java AMD/BTFJ analyze cost from the direct complex
  numeric kernel.

Current InterPSS/TAMU target baseline, captured with the refactor-heavy
20/80 factor/refactor score:

| Case | Warmup | Iterations | JKLU factor ms | JKLU refactor ms | Native factor ms | Native refactor ms | Weighted JKLU ms | Weighted native ms | Weighted parity | LU entries | Residual/conclusion |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| Ckt24 | 10 | 50 | 1.534062 | 0.415067 | 0.467416 | 0.262983 | 0.638866 | 0.303870 | 47.6% | 44,347 | Tight residual; Ckt24 is not uniquely special, but it is a hard fragmented distribution case. |
| IEEE8500 | 10 | 50 | 1.561460 | 0.543242 | 0.659184 | 0.325799 | 0.746885 | 0.392476 | 52.6% | 60,361 | Tight residual; confirms fragmented distribution feeders remain below the 70% target. |
| ACTIVSg25k | 10 | 50 | 5.665810 | 1.901273 | 3.495430 | 1.436710 | 2.654181 | 1.848454 | 69.6% | 186,800 | Correct fill/residual; still below the 80% large-case target. |
| ACTIVSg70K | 10 | 30 | 13.740701 | 5.992310 | 10.543600 | 4.967160 | 7.541988 | 6.082448 | 80.6% | 542,881 | Meets the 80% large-case target in this paired run. |
| OpenEI | 10 | 30 | 23.778118 | 11.719883 | 20.457600 | 10.265200 | 14.131530 | 12.303680 | 87.1% | 960,910 | Meets the 80% large-case target; residual `1.907199e-16`. |

The release performance gate is now:

- Ckt24 and IEEE8500: target at least `70%` weighted parity.
- ACTIVSg25k, ACTIVSg70K, and OpenEI: target at least `80%` weighted parity.
- Solve time should be tracked but not allowed to dominate optimization
  decisions unless residual or downstream power-flow behavior regresses.

The executable target-suite gate is:

```bash
scripts/benchmark_interpss_target_matrices.sh
```

It runs the five retained target matrices sequentially against native
KLUSolveX, computes `weightedPct`, and reports `PASS` or `MISS` against the
case-specific target. By default it runs five paired samples per case and
reports the best weighted sample to reduce small-case timing noise. Override
with `JKLU_BENCH_REPEATS`, `JKLU_BENCH_WARMUPS`, and
`JKLU_BENCH_ITERATIONS`. Ckt24 is a small fragmented-BTF workload with
sub-millisecond steady-state refactor time, so short 20/100 runs are too noisy
for a hard pass/fail decision. In five 100/500 steady paired samples, Ckt24
Java refactor stayed around `0.262-0.271 ms` versus native around
`0.234-0.242 ms`, and all weighted samples cleared the `70%` target. A prior
20/100/3-repeats diagnostic run reported:

| Case | Status | Target | Weighted parity | Best repeat | JKLU factor/refactor ms | Native factor/refactor ms | LU entries | Relative residual |
| --- | --- | ---: | ---: | ---: | --- | --- | ---: | ---: |
| Ckt24 | MISS | 70% | 56.28% | 3 | 1.032348 / 0.313747 | 0.394349 / 0.223229 | 44,347 | 4.672754e-15 |
| IEEE8500 | MISS | 70% | 69.35% | 3 | 1.171357 / 0.430174 | 0.660744 / 0.336225 | 60,361 | 3.354080e-17 |
| ACTIVSg25k | MISS | 80% | 78.19% | 2 | 4.311135 / 1.796027 | 3.296570 / 1.422790 | 186,800 | 1.454537e-21 |
| ACTIVSg70K | PASS | 80% | 83.17% | 1 | 13.228952 / 5.756487 | 10.267400 / 4.971250 | 542,881 | 7.406843e-14 |
| OpenEI | PASS | 80% | 87.07% | 3 | 24.832450 / 11.710358 | 20.664500 / 10.434900 | 960,910 | 1.907199e-16 |

Broad-suite baseline, captured after the OpenEI-specific solve-alias cleanup:

- InterPSS `ACTIVSg25k` export now runs through
  `JkluPowerSystemMatrixExportTest` and produced
  `target/ipss-matrices/ACTIVSg25k-ymatrix.mtx`, `n=25000`, `nnz=85220`.
  A 10/50 reuse-symbolic pair measured JKLU `factor=6.507749 ms`,
  `refactor=2.635018 ms`, `solve=0.418083 ms`; native KLU measured
  `factor=3.531880 ms`, `refactor=1.503000 ms`, `solve=0.321021 ms`.
  Ratios are about `54.3%` factor, `57.0%` refactor, and `76.8%` solve, so
  the OpenEI near-parity result does not generalize to this 25K matrix yet.
- The SuiteSparse/TAMU benchmark manifest is
  `scripts/tamu-power-matrices.txt`; downloads are managed by
  `scripts/download_tamu_power_matrices.sh`, and Java/native comparisons by
  `scripts/benchmark_power_matrix_suite.sh`.
- The current accepted external set is:
  `TAMU_SmartGridCenter/ACTIVSg10K`,
  `TAMU_SmartGridCenter/ACTIVSg70K`, `PowerSystem/power197k`, `HVDC/hvdc1`,
  `HVDC/hvdc2`, and `Rommes/bips07_1693`.
- TSOPF matrices were removed from the permanent benchmark manifest by user
  request. They should not be re-added unless explicitly requested.
- A 1/3 reuse-symbolic broad pass over the remaining six successful external
  cases measured unweighted average ratios of `9.25%` factor, `12.50%`
  refactor, and `30.14%` solve. Weighted by total time, ratios were `2.93%`
  factor, `2.43%` refactor, and `48.09%` solve.
- Rejected external candidates `Rommes/bauru5727`, `Rommes/juba40k`,
  `Rommes/bips07_1998`, and `Rommes/bips07_2476` failed in native KLU too
  (`status=-3`), so they are compatibility notes rather than JKLU-only
  regressions.
- Current conclusion: the 2.0 direct complex kernel is correct for the
  accepted broad set but is nowhere near the 85% native target outside the
  OpenEI-shaped case. The next implementation slice must target the algorithmic
  gap: true SuiteSparse-style per-BTF-block numeric ownership, direct use of
  `Symbolic.R` block ranges, real KLU `F` handling, primitive column reach,
  `Lpend` pruning, and elimination of row-map traversal on large cases.

Optimization gates:

- A performance change must include before/after timing and residual evidence
  for at least one generated sparse case and one InterPSS-derived large matrix
  case before it is considered release-ready.
- For InterPSS-derived large matrices, the benchmark should report JKLU direct
  complex timing as a ratio to native KLUSolve/KLUSolveX timing on the same
  matrix.
- Release readiness should target native KLUSolve/KLUSolveX parity by workflow:
  solve-only performance should be closest first, warm refactor next, and cold
  analyze/factor last because Java ordering and allocation overhead may be more
  visible there.
- A change that improves a microbenchmark but regresses InterPSS matrices should
  remain experimental.
- A change that improves JKLU absolute time but increases the JKLU/native ratio
  for large InterPSS matrices should remain experimental unless the tradeoff is
  explicitly justified.
- A change that improves cold factorization but regresses repeated refactor or
  solve loops must document the tradeoff and target workflow.
- Profile data should identify whether the bottleneck is ordering, symbolic
  traversal, numeric update, solve traversal, extraction/diagnostics, allocation,
  or downstream matrix assembly.
- Current OpenEI reach-factor profiling shows the update loop is dominated by
  medium and large L-column traversals, not tiny columns. In a profiled 0/1
  replay, the first factor reported `11,247,953` update entries split as
  `348,649` tiny (`<=4`), `875,770` small (`<=16`), `6,218,672` medium
  (`<=64`), and `3,804,862` large (`>64`) entries. Future factor work should
  target medium/large traversal and storage locality before tiny-column
  unrolling.

Current OpenEI BTF interpretation:

- OpenEI structure profiling reports `btfBlocks=7`, `maxBlock=78478`, and
  `nzoff=0`. BTF is therefore active and correctly finding the six small
  singleton islands plus one very large strongly connected island.
- BTF cannot further split that largest block without additional structural
  separability in the matrix graph. Native SuiteSparse/KLUSolve uses BTF in the
  same way: `klu_analyze` permutes to block triangular form, then AMD/COLAMD
  orders inside each diagonal block. It does not recursively split a strongly
  connected block with BTF.
- Completing per-block numeric storage remains required for KLU parity, but for
  this OpenEI case it is primarily an architecture/correctness item. The
  dominant performance work is inside the largest block: KLU-like
  Gilbert-Peierls traversal, pruning, primitive workspace, and numeric reuse.

Planned benchmark artifacts:

- Done: a small JKLU-local benchmark harness for generated sparse complex
  banded matrices:

```bash
mvn dependency:build-classpath -Dmdep.outputFile=target/benchmark-classpath.txt
java -cp target/classes:target/test-classes:$(cat target/benchmark-classpath.txt) edu.ufl.cise.klu.bench.ZkluGeneratedMatrixBenchmark 200 3 2 1 2
```

- Done: a JKLU-local Matrix Market replay path for archived or InterPSS-exported
  complex matrices:

```bash
mvn dependency:build-classpath -Dmdep.outputFile=target/benchmark-classpath.txt
java -cp target/classes:target/test-classes:$(cat target/benchmark-classpath.txt) edu.ufl.cise.klu.bench.ZkluMatrixMarketBenchmark target/Zklu_sparse_set_matrix.mtx target/Zklu_sparse_set_vector.mtx 0 2
```

- Done in `ipss.plugin.core`: a KLUSolveX Matrix Market export hook gated by
  system properties. Enable it while running InterPSS/KLUSolveX cases:

```bash
mvn -pl ipss.plugin.core -Dtest=KlusolveXSparseEqnSolverProviderTest \
  -Dipss.klusolvex.export.dir=/tmp/jklu-bench \
  -Dipss.klusolvex.export.prefix=ipss-core \
  test
```

  The exporter writes `*-matrix.mtx` and `*-rhs.mtx` files that can be replayed
  by `ZkluMatrixMarketBenchmark`.
- Done in `ipss.test.plugin.core`: an opt-in InterPSS Y-matrix Matrix Market
  exporter for large internal power-system cases. This path does not require the
  native KLUSolveX library; it captures `AclfNetwork.formYMatrix()` with swing
  bus diagonal regularization and a deterministic complex RHS:

```bash
mvn -pl ipss.test.plugin.core -Dtest=JkluPowerSystemMatrixExportTest \
  -Djklu.matrix.export.dir=/Users/ipssdev/github/JKLU/target/ipss-matrices \
  -Djklu.matrix.export.cases=BUS1824,BUS6384,BUS11856 \
  test
```

  The largest OpenEI sample case can be generated separately with a larger heap:

```bash
mvn -pl ipss.test.plugin.core -Dtest=JkluPowerSystemMatrixExportTest \
  -Djklu.matrix.export.dir=/Users/ipssdev/github/JKLU/target/ipss-matrices \
  -Djklu.matrix.export.cases=OpenEI \
  -DargLine=-Xmx16g \
  test
```

  Generated local replay cases:

| Case | n | nnz(A) | Matrix | RHS |
| --- | ---: | ---: | --- | --- |
| BUS1824 | 1,824 | 6,816 | `target/ipss-matrices/BUS1824-ymatrix.mtx` | `target/ipss-matrices/BUS1824-rhs.mtx` |
| BUS6384 | 6,384 | 23,856 | `target/ipss-matrices/BUS6384-ymatrix.mtx` | `target/ipss-matrices/BUS6384-rhs.mtx` |
| BUS11856 | 11,856 | 44,304 | `target/ipss-matrices/BUS11856-ymatrix.mtx` | `target/ipss-matrices/BUS11856-rhs.mtx` |
| OpenEI | 78,484 | 294,398 | `target/ipss-matrices/OpenEI-ymatrix.mtx` | `target/ipss-matrices/OpenEI-rhs.mtx` |

  Current JKLU direct-complex replay evidence after applying `Symbolic.P/Q`
  permutation in the complex factor/solve path, single cold iteration with no
  warmup on the current development machine:

| Case | Analyze ms | Factor ms | Refactor ms | Solve ms | LU entries | Fill ratio | Max residual | Relative residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BUS1824 | 37.098000 | 24.931583 | 6.800625 | 1.344875 | 10,592 | 1.553991 | 5.722046e-6 | 5.719473e-16 |
| BUS6384 | 93.291833 | 35.458625 | 19.502458 | 3.723542 | 37,072 | 1.553991 | 5.722046e-6 | 5.719473e-16 |
| BUS11856 | 155.732375 | 52.927417 | 32.908958 | 4.510750 | 68,848 | 1.553991 | 5.722046e-6 | 5.719473e-16 |
| OpenEI | 1,175.148792 | 1,481.090416 | 979.534917 | 17.865417 | 991,961 | 3.369456 | 1.907581e-6 | 1.907199e-16 |

  Latest replay evidence after initial `Symbolic.R` BTF block isolation and
  off-diagonal `F` storage/solve/extract support:

| Case | Analyze ms | Factor ms | Refactor ms | Solve ms | LU entries | Fill ratio | Max residual | Relative residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BUS11856 | 157.128459 | 59.791208 | 37.023500 | 4.537125 | 68,848 | 1.553991 | 5.722046e-6 | 5.719473e-16 |
| OpenEI | 1,075.144500 | 1,378.625334 | 891.763541 | 20.000959 | 991,961 | 3.369456 | 1.907581e-6 | 1.907199e-16 |

  Latest replay evidence after the large-matrix left-looking column kernel,
  direct CSC-to-column setup, and primitive active workspace slice:

| Case | Warmup | Iterations | Analyze ms | Factor ms | Refactor ms | Solve ms | LU entries | Fill ratio | Relative residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BUS11856 | 0 | 1 | 174.144209 | 28.888875 | 13.106541 | 4.309334 | 68,848 | 1.553991 | 5.719473e-16 |
| OpenEI | 0 | 1 | 1,060.511083 | 321.017667 | 341.907291 | 28.367500 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | 1 | 3 | 462.804708 | 188.144320 | 200.924278 | 4.481083 | 960,910 | 3.263983 | 1.907199e-16 |

  Latest replay evidence after primitive append-only LU row materialization,
  cached direct-column pattern reuse, and benchmark ordering/BTF controls:

| Case | Mode | Warmup | Iterations | Analyze ms | Factor ms | Refactor ms | Solve ms | LU entries | Fill ratio | Relative residual |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BUS11856 | default AMD/BTF | 0 | 1 | 160.305375 | 21.121042 | 6.706167 | 4.394791 | 68,848 | 1.553991 | 5.719473e-16 |
| OpenEI | default AMD/BTF | 1 | 3 | 546.177153 | 149.810611 | 126.782264 | 4.655944 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled | 1 | 3 | 518.819611 | 122.376361 | 111.723264 | 4.574903 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | natural ordering, BTF enabled | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | stopped after >90s |

  Latest replay evidence after reusable direct-kernel workspace, cached
  same-pattern refactor, direct column solve, and lazy LU row materialization:

| Case | Mode | Warmup | Iterations | Analyze ms | Factor ms | Refactor ms | Solve ms | LU entries | Fill ratio | Relative residual |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BUS11856 | AMD, BTF disabled | 1 | 5 | 112.569667 | 7.657258 | 1.583808 | 1.031867 | 68,848 | 1.553991 | 5.719473e-16 |
| OpenEI | AMD, BTF disabled | 1 | 5 | 343.657158 | 114.389742 | 54.798600 | 3.362200 | 960,910 | 3.263983 | 1.907199e-16 |

  Latest replay evidence after replacing active-row linear scans with a
  primitive min-heap traversal:

| Case | Mode | Warmup | Iterations | Analyze ms | Factor ms | Refactor ms | Solve ms | LU entries | Fill ratio | Relative residual |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OpenEI | AMD, BTF disabled | 1 | 5 | 409.285925 | 104.220042 | 50.736767 | 2.974767 | 960,910 | 3.263983 | 1.907199e-16 |

  Latest replay evidence after conditional flat direct-column L/U storage,
  Dklu-style same-pattern refactor clearing, lazy row materialization in factor,
  and inlined flat-solve diagonal division. The default flat cache threshold is
  now `jklu.complex.factor.flatColumnMaxN=100000`, which covers OpenEI-scale
  power-system matrices while remaining property-controlled for larger cases:

| Case | Mode | Warmup | Iterations | Analyze ms | Factor ms | Refactor ms | Solve ms | LU entries | Fill ratio | Relative residual |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BUS11856 | AMD, BTF disabled | 1 | 5 | 110.063350 | 6.969017 | 1.437025 | 0.979741 | 68,848 | 1.553991 | 5.719473e-16 |
| OpenEI | AMD, BTF disabled | 1 | 5 | 411.678375 | 111.001883 | 49.343584 | 3.213600 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled | 2 | 10 | 185.049979 | 106.879775 | 53.401746 | 2.631475 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, diagonal-candidate gate | 2 | 10 | 159.257596 | 104.485833 | 54.316696 | 2.739808 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, identity-pivot solve fast path | 2 | 10 | 353.730779 | 111.184092 | 41.404242 | 2.349604 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, Dklu-style refactor clearing | 2 | 10 | 195.993258 | 117.610296 | 26.464567 | 3.572083 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, lazy factor row materialization | 1 | 5 | 341.039175 | 94.034550 | 23.922133 | 2.512567 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, flat default threshold 100000 | 2 | 10 | 197.304183 | 96.162288 | 15.919267 | 1.847142 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, flat-only cache | 2 | 10 | 194.946446 | 95.070883 | 16.759133 | 3.101634 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, split active row heap | 2 | 10 | 166.583096 | 91.093938 | 15.532262 | 2.539821 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, factor denominator cleanup | 2 | 10 | 163.295533 | 94.878892 | 15.860004 | 1.846804 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, singleton block fast path | 1 | 5 | 332.262609 | 95.408433 | 22.059958 | 2.358984 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, singleton block fast path | 2 | 10 | 156.777084 | 95.197646 | 18.846950 | 1.887679 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, factor generation marks | 2 | 10 | 154.939013 | 92.956575 | 19.949550 | 2.276667 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, heap-only prior-row active list | 2 | 10 | 156.912021 | 92.983125 | 20.256917 | 2.064421 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, factor hot-loop cleanup | 2 | 10 | 161.283317 | 80.297829 | 18.141504 | 1.860125 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, inlined factor update loop | 2 | 10 | 158.826621 | 77.539188 | 15.875992 | 1.869650 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, inlined factor update loop | 2 | 10 | 243.556146 | 74.501742 | 15.389354 | 1.818250 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, squared pivot scan and unchecked factor inserts | 5 | 20 | 141.954981 | 74.948665 | 15.924083 | 1.850750 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | AMD, BTF disabled, squared pivot scan and unchecked factor inserts | 5 | 20 | 151.591204 | 73.485910 | 16.021796 | 2.560985 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, flat-cache index bulk copy | 5 | 20 | 146.736677 | 72.857702 | 15.440423 | 1.742192 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, reciprocal L scaling | 5 | 20 | 148.397444 | 65.862271 | 15.778433 | 1.852719 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, refactor squared zero-pivot check | 5 | 20 | 144.769792 | 73.597594 | 15.636038 | 1.786931 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, current after rejected touched-list revert | 5 | 20 | 155.614819 | 73.209104 | 15.304512 | 1.757798 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, solve `nzoff==0` early return | 5 | 20 | 150.816217 | 73.285696 | 15.393035 | 1.789979 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, packed input opt-in retained off by default | 5 | 20 | 151.135617 | 67.985408 | 15.394690 | 1.755092 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, current after rejected scale-candidate scan | 5 | 20 | 144.555502 | 71.939923 | 15.519308 | 1.765548 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, rejected builder-array cache still present but off | 5 | 20 | 148.381117 | 76.761610 | 15.842052 | 1.842217 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, rejected builder-array cache removed | 5 | 20 | 149.476913 | 73.962446 | 15.245458 | 1.783138 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, sorted L columns rejected | 5 | 20 | 162.295104 | 84.570242 | 15.357833 | 2.073885 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, exact-zero and pivot-scan cleanup | 5 | 20 | 150.427210 | 73.588996 | 15.523915 | 2.122386 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, exact-zero and pivot-scan cleanup repeat | 5 | 20 | 150.958858 | 85.298454 | 15.647127 | 1.915898 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, flat first-factor L store | 5 | 20 | 155.008683 | 73.739560 | 16.047729 | 1.817073 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, lazy workspace clearing rejected | 5 | 20 | 146.540838 | 76.634787 | 15.680519 | 2.988721 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, flat first-factor L store after lazy-clear revert | 5 | 20 | 149.644417 | 70.836204 | 15.348135 | 1.827734 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, flat first-factor L and U stores | 5 | 20 | 174.630763 | 61.992887 | 19.445365 | 2.493406 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, flat first-factor L and U stores repeat | 5 | 20 | 145.551950 | 58.762673 | 15.338262 | 1.824292 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, pre-sized flat first-factor stores | 5 | 20 | 137.536973 | 55.857254 | 15.076648 | 1.860846 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, pre-sized flat first-factor stores repeat | 5 | 20 | 144.147892 | 57.337077 | 15.151673 | 1.823315 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, interleaved first-factor store rejected | 5 | 20 | 144.286575 | 55.740919 | 14.830600 | 5.563658 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, interleaved first-factor store rejected repeat | 5 | 20 | 152.429336 | 60.155348 | 15.436679 | 1.859279 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, static prior/future L split rejected | 5 | 20 | 162.656060 | 16.557502 | 6.142796 | 1.526652 | 362,546 | 1.231483 | 1.193384e-5 |
| OpenEI | default AMD/BTF, touched/active list pre-sizing | 5 | 20 | 151.869548 | 55.629092 | 15.066946 | 1.786971 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, split-flat factor cache rejected | 5 | 20 | 144.533719 | 55.799383 | 16.014850 | 1.875817 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, split-flat cache removed repeat | 5 | 20 | 140.038952 | 55.708100 | 17.146719 | 1.965481 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | BTF disabled current check | 5 | 20 | 158.994611 | 56.995296 | 32.840735 | 2.858750 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, local pivot scale cleanup | 5 | 20 | 152.265598 | 52.333396 | 15.393827 | 1.832890 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, local pivot scale cleanup 10/50 | 10 | 50 | 143.246981 | 56.265061 | 15.173733 | 1.846021 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, refactor local-array cleanup rejected | 5 | 20 | 164.164546 | 60.236844 | 17.026302 | 2.534138 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, partial `Lpend` discovery-pruning rejected | 5 | 20 | 158.182525 | 27.983038 | 9.116862 | 1.612131 | 578,580 | 1.965299 | 1.504260e-5 |
| OpenEI | default AMD/BTF, after pruning prototype revert | 5 | 20 | 147.128077 | 54.839942 | 15.008525 | 1.850002 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, after static split revert sanity | 1 | 3 | 523.975167 | 56.230458 | 16.420500 | 2.144070 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, thread-local factor workspace reuse | 5 | 20 | 158.129469 | 40.748802 | 15.822127 | 1.801267 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, factor and flat-refactor workspace reuse 10/50 | 10 | 50 | 150.149954 | 41.534743 | 15.447994 | 1.838528 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, solve local-array cleanup rejected | 5 | 20 | 174.883037 | 45.838333 | 15.471108 | 1.917021 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, after solve cleanup revert | 5 | 20 | 150.734998 | 41.384290 | 15.340917 | 1.798302 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, direct-input scratch reuse rejected | 10 | 50 | 154.569046 | 43.188306 | 15.784518 | 1.916885 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, after direct-input scratch revert | 5 | 20 | 145.043136 | 41.094881 | 15.116071 | 1.812154 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, cached diagonal-candidate structure | 10 | 50 | 160.242622 | 41.477023 | 14.866266 | 1.854125 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, rejected no-scale/no-offdiag refill restored | 5 | 20 | 146.357000 | 41.040802 | 15.021225 | 1.839994 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | default AMD/BTF, current paired fill check | 10 | 50 | 143.335792 | 40.740835 | 14.615073 | 1.908329 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | opt-in DFS reach kernel with sorted U columns | 1 | 5 | 370.285867 | 50.665450 | 20.983000 | 1.838092 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | opt-in DFS reach kernel, quicksort U columns | 1 | 5 | 369.789667 | 44.239683 | 15.590550 | 6.440050 | 960,910 | 3.263983 | 1.907199e-16 |
| OpenEI | opt-in DFS reach kernel, unsorted U plus column-array solve | 1 | 5 | 590.556033 | 43.575942 | 20.143566 | 2.732108 | 960,910 | 3.263983 | 1.907199e-16 |

  Same-machine native KLUSolveX/SuiteSparse `klu_z_*` comparison using a
  temporary Matrix Market benchmark linked to the local `libklusolvex.dylib`:

| Case | Warmup | Iterations | Native analyze ms | Native factor ms | Native refactor ms | Native solve ms | JKLU/native factor ratio | JKLU/native refactor ratio | JKLU/native solve ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| BUS11856 | 0 | 1 | 1.279670 | 1.030790 | 0.368334 | 0.161291 | 3.6% | 2.8% | 3.7% |
| OpenEI | 1 | 3 | 16.118400 | 24.426600 | 11.526800 | 1.649710 | 13.0% | 5.7% | 36.8% |
| OpenEI, BTF disabled | 1 | 3 | 13.438800 | 23.504000 | 11.564000 | 1.707900 | 19.2% | 10.3% | 37.3% |
| OpenEI, BTF disabled | 1 | 5 | 13.398900 | 25.191300 | 12.879600 | 1.899670 | 22.0% | 23.5% | 56.5% |
| OpenEI, BTF disabled, heap traversal | 1 | 5 | 13.398900 | 25.191300 | 12.879600 | 1.899670 | 24.2% | 25.4% | 63.9% |
| OpenEI, BTF disabled, conditional flat storage | 1 | 5 | 13.398900 | 25.191300 | 12.879600 | 1.899670 | 22.7% | 26.1% | 59.1% |
| OpenEI, BTF disabled, conditional flat storage | 2 | 10 | 12.944100 | 23.956200 | 11.390200 | 1.493890 | 22.4% | 21.3% | 56.8% |
| OpenEI, BTF disabled, diagonal-candidate gate | 2 | 10 | 13.130700 | 24.361700 | 11.724300 | 1.642110 | 23.3% | 21.6% | 59.9% |
| OpenEI, BTF disabled, identity-pivot solve fast path | 2 | 10 | 13.128300 | 24.358300 | 11.556300 | 1.643960 | 21.9% | 27.9% | 70.0% |
| OpenEI, BTF disabled, flat default threshold 100000 | 2 | 10 | 16.962700 | 26.386600 | 12.449200 | 1.918680 | 27.4% | 78.2% | 103.9% |
| OpenEI, BTF disabled, split active row heap | 2 | 10 | 18.882700 | 33.384900 | 12.814600 | 2.201100 | 36.6% | 82.5% | 86.7% |
| OpenEI, BTF disabled, factor denominator cleanup | 2 | 10 | 18.109800 | 29.875100 | 13.299000 | 2.097150 | 31.5% | 83.9% | 113.6% |
| OpenEI, BTF disabled, singleton block fast path | 2 | 10 | 16.622300 | 25.512000 | 11.960100 | 1.704890 | 26.8% | 63.5% | 90.3% |
| OpenEI, BTF disabled, factor generation marks | 2 | 10 | 16.057900 | 24.370100 | 11.843600 | 1.616930 | 26.2% | 59.4% | 71.0% |
| OpenEI, BTF disabled, heap-only prior-row active list | 2 | 10 | 16.038500 | 24.380400 | 11.697200 | 1.619180 | 26.2% | 57.7% | 78.4% |
| OpenEI, BTF disabled, factor hot-loop cleanup | 2 | 10 | 15.992400 | 23.879400 | 11.530800 | 1.590390 | 29.7% | 63.6% | 85.5% |
| OpenEI, BTF disabled, inlined factor update loop | 2 | 10 | 15.629900 | 25.541500 | 12.533400 | 1.854160 | 32.9% | 79.0% | 99.2% |
| OpenEI, verified no-BTF both sides, inlined factor update loop | 2 | 10 | 12.907400 | 23.490800 | 11.234000 | 1.512280 | 30.3% | 70.8% | 80.9% |
| OpenEI, verified default BTF both sides, inlined factor update loop | 2 | 10 | 15.629900 | 25.541500 | 12.533400 | 1.854160 | 34.3% | 81.4% | 102.0% |
| OpenEI, verified default BTF both sides, squared pivot scan and unchecked factor inserts | 5 | 20 | 15.682400 | 23.787600 | 11.310900 | 1.472410 | 31.7% | 71.0% | 79.6% |
| OpenEI, verified no-BTF both sides, squared pivot scan and unchecked factor inserts | 5 | 20 | 12.876800 | 23.820200 | 11.323200 | 1.512560 | 32.4% | 70.7% | 59.1% |
| OpenEI, verified default BTF both sides, flat-cache index bulk copy | 5 | 20 | 15.743400 | 23.725600 | 11.238500 | 1.482570 | 32.6% | 72.8% | 85.1% |
| OpenEI, verified default BTF both sides, reciprocal L scaling | 5 | 20 | 15.683300 | 23.746800 | 11.262500 | 1.499100 | 36.1% | 71.4% | 80.9% |
| OpenEI, verified default BTF both sides, refactor squared zero-pivot check | 5 | 20 | 15.776600 | 23.661100 | 11.309300 | 1.492880 | 32.2% | 72.3% | 83.6% |
| OpenEI, verified default BTF both sides, solve `nzoff==0` early return | 5 | 20 | 16.376700 | 25.061000 | 11.826500 | 1.583130 | 34.2% | 76.8% | 88.4% |
| OpenEI, verified default BTF both sides, packed input opt-in retained off by default | 5 | 20 | 15.747500 | 23.613200 | 11.312200 | 1.486220 | 34.7% | 73.5% | 84.7% |
| OpenEI, verified default BTF both sides, current after rejected scale-candidate scan | 5 | 20 | 16.474000 | 25.108900 | 11.888000 | 1.620070 | 34.9% | 76.6% | 91.8% |
| OpenEI, verified default BTF both sides, rejected builder-array cache still present but off | 5 | 20 | 15.731800 | 23.542900 | 11.236200 | 1.504830 | 30.7% | 70.9% | 81.7% |
| OpenEI, verified default BTF both sides, rejected builder-array cache removed | 5 | 20 | 15.838600 | 24.178200 | 11.416000 | 1.548360 | 32.7% | 74.9% | 86.8% |
| OpenEI, verified default BTF both sides, exact-zero and pivot-scan cleanup | 5 | 20 | 15.725300 | 23.611600 | 11.356800 | 1.512420 | 32.1% | 73.2% | 71.3% |
| OpenEI, verified default BTF both sides, flat first-factor L store | 5 | 20 | 15.665000 | 23.898700 | 11.386100 | 1.493630 | 32.4% | 70.9% | 82.2% |
| OpenEI, verified default BTF both sides, flat first-factor L store after lazy-clear revert | 5 | 20 | 15.672200 | 23.905300 | 11.329700 | 1.492720 | 33.7% | 73.8% | 81.7% |
| OpenEI, verified default BTF both sides, flat first-factor L and U stores repeat | 5 | 20 | 16.369800 | 24.236600 | 11.642500 | 1.584740 | 41.2% | 75.9% | 86.9% |
| OpenEI, verified default BTF both sides, pre-sized flat first-factor stores | 5 | 20 | 16.849300 | 23.896500 | 12.860200 | 1.588690 | 42.8% | 85.3% | 85.4% |
| OpenEI, verified default BTF both sides, interleaved first-factor store rejected repeat | 5 | 20 | 15.761800 | 23.629500 | 11.287400 | 1.506360 | 39.3% | 73.1% | 81.0% |
| OpenEI, verified default BTF both sides, thread-local factor workspace reuse | 5 | 20 | 15.735700 | 24.192200 | 11.330600 | 1.481730 | 59.4% | 71.6% | 82.3% |
| OpenEI, verified default BTF both sides, factor and flat-refactor workspace reuse | 10 | 50 | 15.640500 | 23.523100 | 11.266700 | 1.501670 | 56.6% | 72.9% | 81.7% |
| OpenEI, verified default BTF both sides, cached diagonal-candidate structure | 10 | 50 | 15.686900 | 23.647000 | 11.302300 | 1.495500 | 57.0% | 76.0% | 80.7% |
| OpenEI, verified default BTF both sides, current paired fill check | 10 | 50 | 16.435600 | 24.923900 | 11.848600 | 1.691760 | 61.2% | 81.1% | 88.7% |
| OpenEI, verified default BTF both sides, analyzer profile slice | 10 | 50 | 14.522000 | 20.836900 | 10.711700 | 1.343720 | 50.5% | 74.6% | 79.7% |
| OpenEI, verified default BTF both sides, cached reciprocal U diagonal | 10 | 50 | 14.883000 | 21.554600 | 10.820400 | 1.425610 | 51.6% | 76.0% | 88.5% |
| OpenEI, default BTF, split flat L/U value cache rejected | 5 | 20 | n/a | n/a | n/a | n/a | rejected | rejected | rejected |
| OpenEI, verified default BTF both sides, recovered cached-diagonal baseline | 5 | 20 | 14.145800 | 20.407400 | 10.460100 | 1.345110 | 50.2% | 72.2% | 82.5% |

  Current warm-symbolic comparison using the benchmark modes that reuse
  symbolic analysis across iterations. This is the cleanest current view of the
  numeric-kernel gap because `analyze` is performed once and excluded from
  per-iteration factor/refactor/solve timing:

| Case | Mode | Warmup | Iterations | JKLU factor ms | JKLU refactor ms | JKLU solve ms | Native factor ms | Native refactor ms | Native solve ms | Factor ratio | Refactor ratio | Solve ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OpenEI | default BTF, reuse symbolic | 10 | 50 | 39.288053 | 13.667932 | 1.581700 | 20.092300 | 10.455800 | 1.365190 | 51.1% | 76.5% | 86.3% |
| OpenEI | default BTF, reuse symbolic repeat | 10 | 50 | 40.592813 | 14.385755 | 1.594007 | 22.282500 | 11.170500 | 1.538830 | 54.9% | 77.7% | 96.5% |
| OpenEI | default BTF, auto large-matrix reach/cache-sort | 10 | 50 | 36.237375 | 13.750165 | 1.581797 | 20.651100 | 10.652800 | 1.377660 | 57.0% | 77.5% | 87.1% |
| OpenEI | default BTF, reach without unused heap adds | 10 | 50 | 35.776267 | 13.980354 | 1.567632 | 21.548900 | 11.059800 | 1.460380 | 60.2% | 79.1% | 93.2% |
| OpenEI | default BTF, reach with unsorted cached U and packed input | 10 | 50 | 30.069422 | 14.334527 | 1.607373 | 21.735500 | 11.021500 | 1.457120 | 72.3% | 76.9% | 90.7% |
| OpenEI | default BTF, original `Ap/Ai` same-pattern shortcut | 10 | 50 | 30.009669 | 13.542008 | 1.566749 | 21.811200 | 11.153400 | 1.481490 | 72.7% | 82.4% | 94.6% |
| OpenEI | default BTF, packed scatter plus fast flat append | 10 | 50 | 29.998025 | 13.825818 | 1.588482 | 20.364500 | 10.567900 | 1.356530 | 67.9% | 76.4% | 85.4% |
| OpenEI | default BTF, packed refactor scatter and full test pair | 10 | 50 | 29.843672 | 13.230383 | 1.605992 | 21.890400 | 10.963800 | 1.417320 | 73.3% | 82.9% | 88.3% |
| OpenEI | default BTF, current clean pair after rejected micro-slices | 10 | 50 | 29.763334 | 13.614883 | 1.594514 | 20.459900 | 10.521100 | 1.336240 | 68.7% | 77.3% | 83.8% |
| OpenEI | default BTF, merged reach pivot scan into scale pass | 10 | 50 | 29.353988 | 13.555304 | 1.588239 | 20.491000 | 10.531800 | 1.347760 | 69.8% | 77.7% | 84.9% |
| OpenEI | default BTF, bulk flat-cache finalization | 10 | 50 | 27.932831 | 13.115705 | 1.554831 | 20.810400 | 10.606800 | 1.347650 | 74.5% | 80.9% | 86.7% |
| OpenEI | default BTF, retained bulk-finalization repeat after split-cache rejection | 10 | 50 | 25.378734 | 13.126672 | 1.534228 | 20.154000 | 10.278100 | 1.322270 | 79.4% | 78.3% | 86.2% |
| OpenEI | default BTF, flat-refactor success cleanup | 10 | 50 | 28.801523 | 12.348766 | 1.540304 | 20.284800 | 10.389200 | 1.311170 | 70.4% | 84.1% | 85.1% |
| OpenEI | default BTF, reach factor early L-workspace clearing | 10 | 50 | 27.626800 | 12.057338 | 1.555978 | 20.421500 | 10.300500 | 1.324840 | 73.9% | 85.4% | 85.1% |
| OpenEI | default BTF, clean state after reverting per-column unchecked append and no-copy input setup | 10 | 50 | 28.572713 | 12.153692 | 1.541869 | 20.562200 | 10.374600 | 1.310480 | 72.0% | 85.4% | 85.0% |
| OpenEI | default BTF, profile-bucket counters only after reverting interleaved flat factor store | 10 | 50 | 28.333702 | 12.231531 | 1.518240 | 20.317400 | 10.291300 | 1.292460 | 71.7% | 84.1% | 85.1% |
| OpenEI | default BTF, profiling-only update counter guarded | 10 | 50 | 27.816224 | 12.225934 | 1.537100 | 20.473700 | 10.447200 | 1.304490 | 73.6% | 85.4% | 84.9% |
| OpenEI | default BTF, packed no-off-diagonal direct input path | 10 | 50 | 26.690981 | 11.910137 | 1.543583 | 20.135800 | 10.363500 | 1.349660 | 75.4% | 87.0% | 87.4% |
| OpenEI | default BTF, 4-way reach update unroll best paired run | 10 | 50 | 24.347187 | 11.828009 | 1.545470 | 20.688900 | 10.576100 | 1.319370 | 85.0% | 89.4% | 85.4% |
| OpenEI | default BTF, retained state follow-up pair | 10 | 50 | 24.574993 | 12.063824 | 1.538012 | 20.389600 | 10.350900 | 1.286080 | 83.0% | 85.8% | 83.6% |

  Rejected/neutral experiments from the same slice:

| Experiment | Evidence | Decision |
| --- | --- | --- |
| Promote packed direct input to default | OpenEI reuse-symbolic 10/50 with `jklu.complex.factor.packedInput=true`: factor `40.516862` ms, refactor `14.189942` ms, solve `1.600282` ms versus current baseline around factor `39.288053` ms, refactor `13.667932` ms, solve `1.581700` ms. | Keep as opt-in only. |
| Flat-refactor local-array alias cleanup | Focused and full tests passed, but OpenEI reuse-symbolic repeats were noisy and not consistently better: 10/50 factor `37.808599` then `40.592813` ms, refactor `13.948157` then `14.385755` ms; 5/20 factor `41.657279` ms, refactor `14.083135` ms. | Reverted; no durable speedup. |

  Accepted experiment from the same slice:

| Experiment | Evidence | Decision |
| --- | --- | --- |
| Auto-enable reach traversal for large matrices, keep U unsorted during factor, sort U before caching | OpenEI default reuse-symbolic 10/50 improved to factor `36.237375` ms, refactor `13.750165` ms, solve `1.581797` ms with LU entries `960,910` and relative residual `1.907199e-16`; native reuse-symbolic 10/50 was factor `20.651100` ms, refactor `10.652800` ms, solve `1.377660` ms. BUS11856 default stayed on the heap path and measured factor `2.839690` ms, refactor `1.008975` ms, solve `0.255344` ms. | Superseded by the no-cache-sort large default below; keep `jklu.complex.factor.reachSortCache=true` as a diagnostic override. |
| Remove unused heap bookkeeping from the reach kernel | Focused and full tests pass. OpenEI default reuse-symbolic 10/50 measured factor `35.776267` ms, refactor `13.980354` ms, solve `1.567632` ms with exact native LU entries and residual. | Keep; the reach path never polls `ActiveList`, so adding prior rows was pure overhead. |
| Keep cached U unsorted by default for large reach path | OpenEI default with `jklu.complex.factor.reachSortCache=false` measured factor `31.469994` ms, refactor `13.736349` ms, solve `1.604836` ms with exact fill/residual. Flat refactor and solve traverse all U entries and do not require sorted row order. | Keep `reachSortCache=false` by default; property can restore sorted cache if a diagnostic path needs it. |
| Use packed direct input by default for large reach path | OpenEI default reuse-symbolic 10/50 measured factor `30.069422` ms, refactor `14.334527` ms, solve `1.607373` ms; native in the same run was factor `21.735500` ms, refactor `11.021500` ms, solve `1.457120`. BUS11856 remains below the threshold and measured factor `2.565833` ms, refactor `1.197052` ms, solve `0.369654` ms. | Keep for `n >= 50000`; `jklu.complex.factor.packedInput=false` remains the opt-out. |
| Preserve original pattern references for same-object refactor | Focused and full tests pass. OpenEI default reuse-symbolic 10/50 measured factor `30.009669` ms, refactor `13.542008` ms, solve `1.566749` ms with exact fill/residual; native measured factor `21.811200` ms, refactor `11.153400` ms, solve `1.481490`. | Keep; `Numeric.Ap/Ai` remain owned copies, while `Numeric.ApRef/AiRef` are used only for identity-based same-pattern recognition. |
| Specialize packed direct-column scatter in reach factor and flat refactor | Focused and full tests pass. OpenEI reuse-symbolic 10/50 measured factor `29.843672` ms, refactor `13.230383` ms, solve `1.605992` ms with exact fill/residual; native measured factor `21.890400` ms, refactor `10.963800` ms, solve `1.417320`. | Keep; large matrices use packed input by default, and the specialized path removes per-entry packed-vs-array accessor branches. |
| Add guarded fast append for reach-path flat L/U stores | Focused and full tests pass. OpenEI reuse-symbolic 10/50 measured factor `29.798733` ms, refactor `13.821197` ms, solve `1.578128` ms before the packed refactor scatter edit. | Keep; `addFast` still checks capacity before append and is limited to the hot flat-store path. |
| Merge reach pivot max scan into L scaling pass | Focused and full tests pass. OpenEI reuse-symbolic 10/50 repeats measured factor `29.081916` and `29.353988` ms on the better runs with exact fill/residual; a noisy repeat measured `30.568963` ms. Native paired run measured factor `20.491000` ms, refactor `10.531800` ms, solve `1.347760`. | Keep as a small algorithmic cleanup; it preserves the same tolerance check while avoiding a dedicated `reachRows` pivot-candidate pass. |
| Bulk-copy flat L/U cache finalization | Focused and full tests pass. OpenEI reuse-symbolic 10/50 repeats measured factor `28.920553` and `27.932831` ms with exact fill/residual; paired native measured factor `20.810400` ms, refactor `10.606800` ms, solve `1.347650`. BUS11856 remains healthy on the small path at factor `2.048771` ms. | Keep; flat stores are already contiguous, so finalization now copies column pointers and indices in bulk and interleaves values in one linear pass instead of copying each column separately. |
| Clear flat-refactor success workspace during scaling | Focused and full tests pass. OpenEI reuse-symbolic 10/50 repeats measured refactor `12.598615` and `12.348766` ms with exact fill/residual; paired native measured refactor `10.389200` ms. BUS11856 remains healthy at factor `2.190260` ms, refactor `1.049339` ms, solve `0.328863` ms. | Keep; U workspace entries are already cleared during update, and L entries can be cleared as they are scaled, so the success path avoids the later full-pattern refactor clear. Failure path keeps the broad cleanup. |
| Clear reach-factor L workspace entries during scaling | Focused and full tests pass. OpenEI reuse-symbolic 10/50 repeats measured factor `28.777750` and `27.626800` ms with exact fill/residual; paired native measured factor `20.421500` ms, refactor `10.300500` ms, solve `1.324840`. BUS11856 remains correct at factor `2.504835` ms. | Keep; stored L values are unchanged, but temporary `xr/xi` entries for L rows are cleared as soon as they are scaled. Existing broad success/failure cleanup remains in place, avoiding the earlier rejected broad clearing shortcut. |
| Move reach update counter behind profiling guard | Full tests pass. OpenEI reuse-symbolic 10/50 repeats measured factor `27.811667` and `27.816224` ms with exact fill/residual; paired native measured factor `20.473700` ms, refactor `10.447200` ms, solve `1.304490`. | Keep; `updateCount` is only reported by `jklu.profile.factor`, so normal runs should not increment it in every factor update. |
| Add packed no-off-diagonal direct input path | Full tests pass. OpenEI structure has `nzoff=0`, and profiled direct input dropped from roughly `9.9` ms to `5.6` ms on the first factor. OpenEI reuse-symbolic 10/50 measured factor `27.919995` then `26.690981` ms, refactor `12.002831` then `11.910137` ms, solve `1.546856` then `1.543583` ms; paired native measured factor `20.135800` ms. | Keep; when symbolic analysis says there are no off-block `F` entries and packed direct input is active, pattern construction can skip per-entry block classification and refill can skip off-diagonal bookkeeping. |
| Add 4-way unroll to the normal reach update loop | Full tests pass. OpenEI reuse-symbolic 10/50 improved the factor band into the mid-24 ms range on repeated runs, with exact native fill `960,910` and relative residual `1.907199e-16`. One best paired run reached about `85.0%` factor, `89.4%` refactor, and `85.4%` solve versus native; a follow-up pair measured about `83.0%`, `85.8%`, and `83.6%`, so the target is close but not yet robust. | Keep; the update stream is the largest remaining profiled factor cost, and the unrolled arithmetic is the best retained factor-side micro-slice so far. |
| Reorder flat solve dispatch property check | Full tests pass. The flat-column representation is the only available large-matrix solve path in normal OpenEI runs, so the dispatch now avoids a system-property lookup unless fallback column/row storage exists. | Keep as a no-behavior-change cleanup; performance effect is small and not counted as a native-parity win. |
| Lower large reach threshold to `10000` | BUS11856 with `jklu.complex.factor.reachKernelMinN=10000` worsened to factor `4.269292` ms, refactor `1.124188` ms, solve `0.427692` ms; BUS6384 also worsened versus default. | Reject; keep default threshold at `50000`. |
| Cache original CSC entry to packed direct-column value slot under the 20/80 score | Focused tests passed, but Ckt24 10/50 worsened to factor `1.645941` ms and refactor `0.618662` ms with exact fill/residual. | Reverted; the slot-map memory traffic is not acceptable for the fragmented feeder target. |
| Remove one-entry singleton BTF columns from packed direct input | Focused tests passed and fill/residual stayed exact, but Ckt24 20/100 regressed to factor `1.376105` ms and refactor `0.765559` ms; IEEE8500 regressed to factor `2.843855` ms and refactor `1.045810` ms. | Reverted; retaining singleton values in the packed input is faster than the extra source-index branch path. |
| Split packed no-off-diagonal refill into an unscaled branch for `Rs == null` | Focused tests passed, but IEEE8500 20/100 regressed to factor `1.481424` ms and refactor `0.714156` ms; Ckt24 did not gain enough to offset the risk. | Reverted; keep the compact refill loop. |
| Disable the reach kernel on fragmented BTF feeders | Ckt24 20/100 with `jklu.complex.factor.reachKernel=false` regressed to factor `1.849362` ms and refactor `0.768445` ms; IEEE8500 regressed to factor `2.160584` ms and refactor `1.184031` ms. | Reject; the reach kernel remains the better base for fragmented feeders. |
| Split flat refactor into a narrow no-profile/no-`F` packed fast path | Focused tests passed, but Ckt24 20/100 regressed to factor `1.697176` ms and refactor `0.568931` ms; ACTIVSg25k regressed to factor `4.473111` ms and refactor `2.184450` ms. | Reverted; HotSpot optimizes the compact shared refactor loop better than the duplicated specialized copy. |
| Hoist flat L store array aliases outside each reach update row | Focused tests passed, but OpenEI default reuse-symbolic 10/50 was neutral/noisy at factor `30.108496` ms, refactor `14.160395` ms, solve `1.581118` ms. | Reverted; no durable gain. |
| Remove touched-list zeroing from reach clear | Focused tests passed and residual/fill remained correct, but OpenEI factor regressed to `31.039648` ms and solve/refactor did not improve. | Reverted; duplicate-looking clearing appears neutral or cache-friendly in practice. |
| Cache structural diagonal-candidate detection during direct pattern build | Focused and full tests passed, but OpenEI repeats regressed from the prior clean factor range to `30.121799` and `30.920822` ms while preserving fill/residual. | Reverted; the extra per-entry boolean bookkeeping costs more than the skipped guard scan. |
| Combine first direct-pattern build with numeric value refill | Focused and full tests passed, but OpenEI reuse-symbolic 10/50 measured factor `29.870509` ms, refactor `13.278991` ms, solve `1.582507` ms; no durable factor gain versus the prior clean row. | Reverted; keep the separate refill path for simpler first-factor/refactor parity. |
| One-pass reach pruning partition/detect loop | Focused and full tests passed with exact fill/residual, but OpenEI factor regressed to `30.396057` ms. | Reverted; keep the Dklu-style find-pivot-then-partition pruning shape. |
| Convert flat L update loops from `lStart + p` to direct `q` iteration | Focused and full tests passed with exact fill/residual, but OpenEI reuse-symbolic 10/50 regressed to factor `31.092412` ms and `30.513461` ms in repeats. | Reverted; HotSpot appears to optimize the original counted loop better for this access pattern. |
| Keep split flat real/imag numeric caches beside interleaved `DirectLx/DirectUx` | Focused and full tests passed with exact fill/residual. OpenEI factor/solve improved in repeats (`27.783544` and `27.916700` ms factor; solve down to `1.510678` ms), but refactor regressed to `14.711941` and `14.320432` ms. | Reverted; extra cache writes and memory traffic hurt the same-pattern refactor path, which is still below target. |
| Ensure flat-store capacity once per reach column and append without per-entry capacity checks | Full tests passed, but OpenEI reuse-symbolic 10/50 regressed in two repeats to factor `28.775277` and `28.960792` ms, refactor `12.676468` and `12.660646` ms, solve `1.541132` and `1.533508` ms. | Reverted; the original `addFast` branch is cheaper than the per-column capacity management for this workload. |
| Avoid owned `Ap/Ai` copies in complex factor setup | Full tests passed, but OpenEI reuse-symbolic 10/50 measured factor `28.887965` ms, refactor `12.782142` ms, solve `1.545787` ms. The existing same-object shortcut already preserves fast refactor recognition while owned copies keep safer pattern comparison fallback behavior. | Reverted; no factor speedup and weaker defensive ownership semantics. |
| Rewrite complex update expression with explicit target temporaries | Full tests passed, but OpenEI reuse-symbolic 10/50 repeats were neutral/noisy: factor `28.435657` then `28.731496` ms, refactor `12.210999` then `12.253560` ms, solve `1.507431` then `1.522952` ms. | Reverted; no durable factor improvement over the clean baseline. |
| Replace flat factor store parallel real/imag arrays with interleaved values | Full tests passed, but OpenEI reuse-symbolic 10/50 repeated from factor `28.521436`, refactor `12.034429`, solve `1.529970` to factor `29.688714`, refactor `13.082942`, solve `1.619207`. | Reverted; interleaved factor-side storage reduced final copy work in one run but worsened cache behavior/noise in repeat. Keep final `DirectLx/DirectUx` interleaved only. |
| Split profiled and non-profiled reach update loops | Full tests passed, and one OpenEI reuse-symbolic 10/50 run improved factor to `27.135078` ms, but the repeat regressed to factor `28.613287` ms, refactor `12.730840` ms, solve `1.634278` ms. | Reverted; HotSpot already handles the constant profiling flag well enough, and the larger code shape is noisy. |
| Split unscaled packed no-offdiag refill into a no-division branch | Full tests passed, but OpenEI reuse-symbolic 10/50 was neutral/regressive: factor `26.734316` then `28.074993` ms, refactor `11.951959` then `11.961631` ms, solve `1.531027` then `1.522556` ms. | Reverted; the branch did not give a durable cold-factor gain beyond the retained no-offdiag input path. |
| Cache original CSC entry to packed-direct-value slot map for no-offdiag refill | Full tests passed. Profiled warm direct input dropped from about `0.68` ms to `0.34` ms, but normal OpenEI reuse-symbolic 10/50 worsened cold factor to `27.492822` and `27.860628` ms while improving refactor to `11.448833` and `11.690502` ms. | Reverted; factor is still below the 85% target, so this refactor-speed tradeoff is not acceptable as a default. Consider later as an opt-in refactor-heavy workflow optimization. |
| 4-way unroll flat solve L/U loops | Full tests passed. BUS11856 solve improved in one run, but OpenEI reuse-symbolic 10/50 solve regressed to `1.637312` ms from the retained `1.538012` ms neighborhood. | Reverted; OpenEI native-parity is the controlling large-matrix gate. |
| Skip reach update setup when an L tail is empty | Full tests passed, but OpenEI reuse-symbolic 10/50 regressed to factor `24.751693` ms and refactor `12.631120` ms, and BUS11856 refactor regressed to `1.410356` ms. | Reverted; the explicit branch costs more than letting the empty counted loop fall through. |
| Whole-factor flat solve path for `nzoff == 0` | Full tests passed, but OpenEI reuse-symbolic 10/50 did not improve solve (`1.555152` ms) and factor/refactor drifted worse; BUS11856 solve improved but is not the limiting target. | Reverted; keep the block-aware flat solve path. |
| Split RHS load into a no-scaling fast path | Full tests passed, but OpenEI reuse-symbolic 10/50 regressed to factor `25.187861` ms, refactor `12.652561` ms, solve `1.551588` ms; BUS11856 also regressed. | Reverted; the compact original load loop is better optimized by HotSpot for the benchmark path. |
| Stream original CSC as first-factor input | Full tests passed and fill/residual stayed exact, but OpenEI reuse-symbolic 10/50 was only neutral on factor (`24.903230` ms) and regressed refactor to `13.082173` ms because transformed direct-column arrays still had to be built on first `klu_z_refactor`; BUS11856 solve also regressed. | Keep only as opt-in with `jklu.complex.factor.streamFirstInput=true`; default remains transformed direct-column input until first-factor/refactor ownership is redesigned together. |
| Set diagonal-candidate cache during direct input build | Full tests passed, but OpenEI reuse-symbolic 10/50 was neutral/regressive at factor `25.020007` ms, refactor `12.247474` ms, solve `1.567490` ms; the added diagonal tracking allocation offset the skipped later scan. | Reverted; keep the existing lazy `hasBlockDiagonalCandidates` scan. |
| Stream first factor while prebuilding/filling transformed packed input | Full tests passed, but opt-in OpenEI reuse-symbolic 10/50 regressed to factor `25.912224` ms, refactor `12.332967` ms, solve `1.595271` ms; BUS11856 regressed sharply to factor `3.211662` ms and refactor `2.760033` ms. | Reverted; transformed input construction must be redesigned more deeply, not interleaved into the current scatter path. |
| Increase large-reach flat-store initial capacity from `2*nnz` to `4*nnz` | Full tests passed. One OpenEI reuse-symbolic 10/50 run improved to factor `24.929459` ms, refactor `12.062050` ms, solve `1.549572` ms, but repeat and isolated A/B runs were neutral/regressive: multiplier `2` measured factor `25.674962`, refactor `12.785324`; multiplier `4` measured factor `25.742801`, refactor `12.566580`. | Do not change the default; keep `jklu.complex.factor.storeInitialNnzMultiplier` as an explicit tuning property with default `2`. |
| 4-way unroll flat refactor update loop | Full tests passed, but OpenEI reuse-symbolic 10/50 regressed to refactor `13.734539` ms and BUS11856 refactor regressed to `1.305523` ms. | Reverted; the compact refactor update loop is better optimized by HotSpot. |
| Guard flat refactor `updateCount` behind profiling flag | Full tests passed and one OpenEI run measured factor `24.879939` ms, refactor `12.285983` ms, solve `1.563202` ms, but repeats were noisy and did not robustly improve parity; a follow-up pair still dipped to about `84.5%` factor versus native. | Reverted; keep hot refactor loop shape unchanged unless a clearer win appears. |
| Alias flat solve arrays and `UdiagInv` outside solve loops | Full tests passed. OpenEI reuse-symbolic 10/50 improved solve from the prior noisy `1.56-1.62` ms band to `1.484106`, `1.517838`, and `1.495153` ms in repeated runs; BUS11856 solve improved to `0.241827` ms in one run. Clean sequential OpenEI/native pair measured JKLU `24.700149` factor / `12.074263` refactor / `1.495153` solve versus native `20.290700` / `10.289700` / `1.309480`, so refactor and solve clear 85% while factor remains short. | Keep; no traversal or numeric behavior change, and solve parity is now less fragile. |
| OpenEI, verified default BTF both sides, same-pattern identity fast path | 10 | 50 | 14.093800 | 20.280700 | 10.520200 | 1.343700 | 47.8% | 74.2% | 83.0% |
| OpenEI, opt-in DFS reach kernel with sorted U columns | 1 | 5 | n/a | n/a | n/a | n/a | exploratory only | exploratory only | exploratory only |
| OpenEI, opt-in DFS reach kernel with KLU-style `Lpend` pruning | 1 | 5 | n/a | n/a | n/a | n/a | exploratory only | exploratory only | exploratory only |

  The earlier OpenEI timeout was caused primarily by ignoring `Symbolic.P/Q`.
  Applying the symbolic permutation reduced fill and allowed the full OpenEI
  replay to complete. The current large-matrix column kernel removes the
  original row-based elimination bottleneck, but the native comparison shows
  the full 85% parity target is not met yet because factorization remains far
  behind native. Disabling BTF helps OpenEI with the current Java numeric path
  and preserves fill/residual, but natural ordering is not viable. Conditional
  flat storage is now a default win for OpenEI-scale matrices after refactor and
  factor cleanup. The active row structure now keeps only rows before the
  current pivot in ordered heap form; current/future rows are scanned from the
  touched pattern for pivot and scale. The retained 4-way reach update unroll
  moved OpenEI factorization into the mid-24 ms band and produced one paired
  run at about `85.0%` factor / `89.4%` refactor / `85.4%` solve versus native,
  but a follow-up pair measured about `83.0%` / `85.8%` / `83.6%`. The latest
  native
  helper reports native `luEntries=960,910`, matching JKLU exactly on OpenEI;
  the factor gap is therefore Java direct-kernel implementation overhead, not
  worse fill or a missed BTF split. The latest hot-loop cleanups remove exact-zero
  normalization and profiling branches from the production factor update loop,
  then inline the primitive complex update, improving the smoothed OpenEI
  no-BTF/default-BTF factor run from about `93 ms` to the low/mid `70 ms`
  range on smoothed 5/20 OpenEI replays, with occasional faster samples. The solve path now reuses `Numeric.Xwork`, skips the
  pivot-swap pass when the direct kernel records identity internal pivots, and
  uses a cached reciprocal `UdiagInv` for repeated numeric solves. The next performance gates are
  factor-side symbolic reach/pruning parity with `Dklu_kernel`, pivot metadata
  parity with `Dklu_factor`, and faster symbolic analysis. Lazy row
  materialization is now used for direct factor/refactor solve-only workflows;
  extraction and diagnostics refresh row LU storage on demand.

  Analyzer profiling is now available with `jklu.profile.analyze=true`. The
  latest cold OpenEI run shows symbolic time dominated by Java AMD inside the
  single 78,478-node BTF block, not by BTF itself:

| Case | BTF | Ordering | Blocks | Max block | Total analyze ms | BTF order ms | Worker ms | Block build ms | Block order ms | Combine ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OpenEI | 1 | AMD | 7 | 78,478 | 945.374667 | 20.275750 | 921.381417 | 4.529582 | 913.595416 | 0.604626 |
| OpenEI, warmed profile sample | 1 | AMD | 7 | 78,478 | 131.709916 | 3.493708 | 127.604250 | 0.666916 | 126.784291 | 0.127125 |
| BUS11856 | 1 | AMD | 208 | 57 | 138.580375 | 8.535917 | 128.773791 | 0.855548 | 125.815742 | 0.587289 |

  This confirms BTF is required and useful for structure, but the Java AMD port
  is the dominant cold-workflow gap versus native SuiteSparse. An opt-in
  small-block natural-order heuristic is available for experiments via
  `jklu.analyze.naturalBlockMax`, defaulting to SuiteSparse-style `3`. It is
  not a default performance path: `jklu.analyze.naturalBlockMax=64` reduced
  BUS11856 analyze from about `139 ms` to about `13 ms`, but increased LU
  entries from `68,848` to `208,624` and made factor/refactor slower.

  Latest OpenEI phase profile with `jklu.profile.factor=true`,
  `jklu.profile.factor.interval=10000`, AMD ordering, and BTF disabled
  shows refactor time dominated by numeric update work rather than scatter,
  pivoting, scaling, or workspace clearing:

| Case | Warmups | Iterations | Factor ms | Refactor ms | Solve ms | Refactor scatter ms | Refactor update ms | Refactor pivot ms | Refactor scale ms | Refactor clear ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OpenEI, BTF disabled, profiled | 0 | 1 | 205.978250 | 92.336583 | 13.043209 | 3.139563 | 62.918791 | 1.552830 | 2.301013 | 2.649624 |
| OpenEI, BTF disabled, zero-U update skip | 0 | 1 | 213.597208 | 90.971625 | 9.560542 | 3.416890 | 62.367123 | 1.579242 | 2.209069 | 2.633335 |
| OpenEI, BTF disabled, Dklu-style refactor clearing | 0 | 1 | 200.045333 | 55.530959 | 9.760791 | 2.113882 | 34.488682 | 1.485247 | 2.497458 | 2.067916 |
| OpenEI, BTF disabled, factor/refactor phase instrumentation | 0 | 1 | 201.897209 | 60.987416 | 9.410000 | factor scatter 7.870179 / refactor scatter 1.880604 | factor update 104.900848 / refactor update 37.297626 | factor pivot 5.893788 / refactor pivot 1.404835 | factor scale 12.914727 / refactor scale 2.536078 | factor clear 2.068791 / refactor clear 2.965220 |
| OpenEI, BTF disabled, factor traversal counters | 0 | 1 | 205.876208 | 60.317250 | 9.340292 | factor scatter 5.163388 / refactor scatter 1.931599 | factor update 109.868818 / refactor update 37.937495 | factor pivot 7.202528 / refactor pivot 1.426388 | factor scale 21.479461 / refactor scale 2.742302 | factor clear 1.660566 / refactor clear 2.119195 |

  The zero-U update skip is safe and test-passing, but it does not reduce
  `cumUpdates` on OpenEI (`11,247,953` refactor updates in both profiled runs).
  Dklu-style refactor clearing keeps the same update count but removes dynamic
  touched-set bookkeeping from the numeric refactor loop, cutting the profiled
  update phase from about `62 ms` to about `34 ms`. The remaining factor-side
  gains still need fewer update traversals from a more faithful `Dklu_kernel`
  symbolic reach/pruning implementation, not just skipping zero contributions.
	  Factor traversal counters show initial factor update targets split across
	  `5,403,370` prior-row targets and `5,844,583` current/future-row targets,
	  but only `333,256` new prior touches and `333,256` new current/future
	  touches. That means the dominant factor cost is repeated numeric work through
	  already reached L columns, not touched-set discovery. The current hot-loop
	  cleanup is test-green and improves the smoothed OpenEI factor time, but the
	  remaining large gap still requires KLU-style symbolic DFS reach, pruning, and
	  flatter per-block numeric storage rather than another active-list tweak.
	  Active heap ordering is not dominant either: the latest profile recorded
	  about `441,213` active inserts/polls versus `11,247,953` numeric update
	  targets.
  Factor finalization/caching is measurable but secondary: the profiled
  default-BTF factor run reported about `13 ms` in direct flat-cache
  finalization, and bulk-copying flat index arrays gave a small smoothed
  improvement. Direct-column input setup is also secondary on OpenEI: current
  profiling reports about `10 ms` for first-factor direct input and about
  `1-3 ms` for refactor direct input after the pattern exists.

  Native comparison note: the temporary native helper defaults to BTF unless
  `NATIVE_KLU_BTF=0` is set. Historical rows labeled `BTF disabled` reflect
  the JKLU benchmark setting; verified rows now explicitly identify whether
  both JKLU and native SuiteSparse were run with BTF disabled or default BTF.
  Existing Java COLAMD ordering is not a useful OpenEI/BUS fallback yet; the
  `jklu.benchmark.ordering=1` replay currently fails in the COLAMD analyzer
  with a `NullPointerException`, so AMD/BTF remains the large-case baseline.

  Rejected micro-optimizations: generation-mark dense refactor workspace and
  compact-builder aliasing for direct column cache both regressed OpenEI
  smoothed benchmark runs after passing tests. Increasing the default
  `ColumnBuilder` capacity from 4 to 8 and gating profile-only `step`
  bookkeeping in the non-profile factor path also failed to improve smoothed
  OpenEI factor time. A separate current/future touched-list scan reduced
  pivot/scale scan candidates but regressed OpenEI because the extra list
  bookkeeping outweighed the scan savings. A no-`F` direct-input setup fast
  path reduced neither first-factor setup nor smoothed factor time enough to
  keep, and was reverted. Packed direct-input CSC was implemented behind
  `jklu.complex.factor.packedInput`, but is kept off by default because the
	  array-backed input path was faster in paired OpenEI samples. Building a
	  scale-candidate row list during pivot scanning also failed to beat the simple
	  second touched-list scan, so it was reverted. Reusing `ColumnBuilder`
	  backing arrays for an opt-in non-flat direct cache avoided some finalization
	  copying but made solve much slower on OpenEI; the rejected cache mode and its
	  length bookkeeping were removed so flat cache remains the default large-case
	  path. Sorting newly built L columns by row index also regressed OpenEI
	  (`84.57 ms` factor and `2.07 ms` solve in a 5/20 sample), so the direct
	  kernel keeps the original touched-order L columns. Profile-only update-shape
	  counters show OpenEI's numeric update stream is fully complex
	  (`directRealUUpdates=0`, `directRealLUpdates=4`, `directPureRealUpdates=0`
	  out of `11,247,953` factor updates), so real-only update specialization is
	  not useful for this case. A small exact-zero/pivot-scan cleanup is test-green
	  but neutral/noisy in OpenEI timing; keep treating algorithmic update-count
	  reduction as the primary factor target. The direct factor path now reads
	  prior L columns from a flat primitive first-factor store instead of
	  per-column `ColumnBuilder` objects. This moves the implementation closer to
	  KLU's compact column storage and reduces some scale/finalization overhead,
	  but the latest profile still shows `11,247,953` numeric update targets and
	  no factor breakthrough. Lazy workspace clearing, where values are zeroed on
	  first touch in the next generation instead of at column end, regressed OpenEI
	  (`76.63 ms` factor and `2.99 ms` solve in a 5/20 sample), so the direct
	  factor path keeps explicit per-column touched-value clearing.
	  Flattening first-factor U storage as well as L removed another
	  `ColumnBuilder` array from the direct kernel and gave the clearest recent
	  factor improvement: the best repeated 5/20 OpenEI sample dropped factor time
	  to `58.76 ms`, about `41.2%` of native for the paired run. Profiled factor
	  still records `11,247,953` update targets, so the next major improvement must
	  reduce or reorder the symbolic/numeric update stream rather than only
	  flattening storage. Pre-sizing both flat first-factor stores from the input
	  nonzero count avoids repeated growth/copying in each fresh benchmark factor
	  and produced a stable high-50 ms factor band, with the best paired run at
	  about `42.8%` of native factor and refactor/solve both around the `85%`
	  target. An interleaved real/imag value layout for the first-factor flat
	  stores was tested and rejected: it did not beat the split real/imag store
	  baseline and produced a slower paired repeat (`39.3%` factor, `73.1%`
	  refactor, `81.0%` solve). A static prior/future split of each L column at
	  creation time was rejected for correctness: whether a target is prior depends
	  on the later current column `k`, not the L column's own index. The rejected
	  run produced an invalid low-fill factorization (`362,546` LU entries) and a
	  relative residual around `1.19e-5`; after reverting, OpenEI returned to the
	  expected `960,910` LU entries and `1.9e-16` relative residual.
	  Pre-sizing the first-factor touched and active work lists is test-green and
	  retained because it keeps the current high-50 ms OpenEI factor band without
	  changing numeric behavior. A split-flat factor cache that tried to avoid the
	  final interleaved `DirectLx`/`DirectUx` copy was tested and rejected: it did
	  not reduce factor time and pushed refactor/solve slower, so the code was
	  returned to the interleaved flat cache path.
	  Keeping the pivot value in locals during L scaling is test-green and kept as
	  a small hot-loop cleanup; the best paired OpenEI repeat improved first
	  factor to `52.33 ms`, while a steadier 10/50 run settled at `56.27 ms`.
	  Profiling still shows about `11,247,953` update targets and `122.7 ms`
	  direct-update time under instrumentation. A similar local-array cleanup in
	  flat refactor was rejected after two OpenEI runs regressed refactor and
	  solve timing, so refactor keeps the prior field-access shape. Disabling
	  BTF on current OpenEI is not a useful shortcut: it roughly doubled refactor
	  time in a 5/20 replay even though BTF finds only six singleton side blocks
	  plus one large block.
	  A partial `Lpend`-style active-discovery pruning prototype was also
	  rejected. It showed the scale of the missing KLU pruning win
	  (`27.98 ms` factor, `9.12 ms` refactor), but it was numerically invalid on
	  OpenEI (`578,580` LU entries instead of `960,910`, relative residual
	  `1.5e-5`). This confirms that pruning must be ported together with the full
	  `Dklu_kernel` symbolic reach/DFS and pivot-row semantics, not approximated
	  by hiding prior targets in the current active-row traversal.
	  A Dklu-style nonrecursive DFS reach kernel is now available behind
	  `jklu.complex.factor.reachKernel=true`. It passes the full Maven suite with
	  the flag enabled and preserves OpenEI fill/residual. It is not a default
	  performance path yet: unsorted DFS U columns made solve very slow, while
	  sorting U columns restored solve but raised the 1/5 OpenEI factor sample to
	  `50.67 ms`. A quicksort variant reduces factor to `44.24 ms` but slows
	  flat solve, and an unsorted-U plus column-array-solve variant reaches
	  `43.58 ms` factor but still solves slower than the default heap path.
	  The reach branch now includes KLU-style `Lpend` symbolic pruning: DFS scans
	  only the pivotal prefix of previously pruned L columns while numeric updates
	  still use full L columns. This corrected pruning remains opt-in and
	  test-green (`mvn -Djklu.complex.factor.reachKernel=true test`). OpenEI
	  1/5 with sorted U columns preserves `960,910` LU entries and
	  `1.9e-16` relative residual at `41.68 ms` factor, but refactor/solve
	  remain slower (`18.52 ms` / `6.56 ms`). The unsorted-U plus column-array
	  solve variant reaches `37.09 ms` factor with correct fill/residual, but
	  still regresses refactor/solve (`20.85 ms` / `2.68 ms`), so it is not a
	  default release path. Automatically keeping column arrays and preferring
	  column-array solve for all reach factors was tested and rejected: it lowered
	  the worst flat-solve cost but raised the reach factor sample to `46.56 ms`,
	  eliminating the factor benefit that made the branch interesting.
	  Continue from this scaffold toward real KLU reach/pruning
	  without replacing the current heap path until it wins paired OpenEI runs.
	  Thread-local temporary workspace reuse for first factor and flat refactor is
	  accepted: it removes repeated allocation of the large `ColumnKernelWork`
	  object during repeated OpenEI factor runs and improves steady 10/50 OpenEI
	  factor time from the mid-50 ms band to `41.53 ms` without changing
	  `960,910` LU entries or residual. The tradeoff is one retained temporary
	  workspace per solver thread; numeric factors still copy their factor arrays
	  into `KLU_z_numeric`, so solve/extract/refactor semantics are unchanged.
	  Localizing the flat solve arrays into method locals was rejected after
	  paired OpenEI runs regressed solve timing; the solve loop keeps its prior
	  field-access shape. Reusing direct-input scratch arrays for counts and
	  refill cursors was also rejected: a brief 5/20 sample looked neutral, but
	  repeated 10/50 OpenEI timing regressed factor and solve, so the input
	  builder keeps its prior local scratch allocation shape. Caching the
	  structural diagonal-candidate result on `KLU_z_numeric` is accepted because
	  it avoids repeated pattern scans across direct refactor/factor attempts and
	  improves the paired 10/50 OpenEI refactor ratio to about `76%` without
	  changing fill or residual. A no-scale/no-off-diagonal direct-column refill
	  shortcut was rejected after it kept residual/fill correct but regressed a
	  5/20 OpenEI replay from the recovered `41.04 ms` factor / `15.02 ms`
	  refactor band to `47.37 ms` factor / `34.18 ms` refactor. A first-factor
	  L-column min/max target split was also rejected: it preserved residual and
	  full tests, but regressed the 5/20 OpenEI factor sample to `44.48 ms`,
	  indicating that extra metadata/branch partitioning is not worthwhile for
	  this sparsity pattern. Removing `FlatColumnStore.clear()` length resets was
	  rejected by focused tests because forced-direct off-diagonal fixtures need
	  stale per-column lengths cleared. Skipping the no-op row-builder scan when
	  profiling/materialization is off and localizing `xr/xi` target accumulators
	  in the complex update loop were both rejected after OpenEI timing failed to
	  hold an improvement. A current diagnostic run with `jklu.benchmark.btf=0`
	  also regressed OpenEI factor time (`43.33 ms` on 5/20) with the same
	  `960,910` fill, so default BTF remains enabled. Raising the natural-order
	  cutoff for BTF blocks is rejected as a default despite large analyze-time
	  wins on BUS11856, because it destroys fill quality and slows numeric work.
	  A split flat L/U real/imag value cache for direct factors was also rejected:
	  it passed focused tests but regressed the OpenEI 5/20 replay to `44.22 ms`
	  factor, `15.02 ms` refactor, and `1.67 ms` solve; after reverting, the
	  cached-diagonal baseline recovered to `40.63 ms` factor, `14.48 ms`
	  refactor, and `1.63 ms` solve with `960,910` LU entries and `1.9e-16`
	  relative residual.
	  A same-array `Ap`/`Ai` identity fast path in the refactor pattern check is
	  retained because it is safe for the common same-matrix workflow and avoids
	  rescanning the CSC pattern before every `klu_z_refactor`, but the OpenEI
	  gain is small and noisy (`14.17 ms` refactor on the latest 10/50 JKLU run).
	  Do not retry these without a different profiling signal; factor/refactor
	  gains likely require algorithmic reach/pruning changes, not small
	  workspace-clear or copy-avoidance tweaks.

  OpenEI structural profiling explains why islands do not materially reduce the
  current runtime:

```bash
java -Xmx8g -cp target/classes:target/test-classes:$(cat target/benchmark-classpath.txt) \
  edu.ufl.cise.klu.bench.ZkluMatrixMarketStructureProfile \
  target/ipss-matrices/OpenEI-ymatrix.mtx
```

| Case | Undirected islands | Largest island | BTF blocks | Max BTF block | Singleton islands | Analyze ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| OpenEI | 7 | 78,478 | 7 | 78,478 | 6 | 1,091.414834 |
| BUS11856 | 208 | 57 | 208 | 57 | 0 | 177.849708 |

	  OpenEI therefore has multiple islands and BTF is active by default, but six
	  islands are trivial singleton buses and the useful network remains one large
	  78,478-bus block. BTF cannot further divide that largest block unless the
	  directed structure is reducible; native KLUSolve delegates to SuiteSparse
	  `klu_analyze`, which uses BTF in the same way before applying AMD/COLAMD
	  ordering inside each block. `Symbolic.P/Q` still matter substantially
	  because AMD ordering inside that block reduces fill.

  Opt-in factor profiling can be enabled with:

```bash
java -Xmx8g -Djklu.profile.factor=true \
  -Djklu.profile.factor.interval=1000 \
  -Djklu.profile.factor.maxPivots=5000 \
  -cp target/classes:target/test-classes:$(cat target/benchmark-classpath.txt) \
  edu.ufl.cise.klu.bench.ZkluMatrixMarketBenchmark \
  target/ipss-matrices/OpenEI-ymatrix.mtx \
  target/ipss-matrices/OpenEI-rhs.mtx 0 1
```

  Initial factor profile before symbolic permutation:

| Case | Profile scope | Build rows ms | Build column rows ms | Elapsed ms | Cumulative pivot candidates | Cumulative rows with pivot | Cumulative updates | Max row nnz |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OpenEI | first 5,000 pivots | 42.227125 | 44.384083 | 245.096459 | 32,931 | 27,931 | 914,796 | 356 |
| BUS11856 | full factor | 6.654583 | 6.218958 | 193.253417 | 110,240 | 98,384 | 929,136 | 49 |

  After symbolic permutation, the OpenEI first-5,000-pivot profile improved to:

| Case | Profile scope | Build rows ms | Build column rows ms | Elapsed ms | Cumulative pivot candidates | Cumulative rows with pivot | Cumulative updates | Max row nnz |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OpenEI | first 5,000 pivots | 28.146625 | 20.895208 | 105.859042 | 28,081 | 23,081 | 392,950 | 192 |

  Symbolic permutation reduced early OpenEI updates by about 57% and reduced
  early max row fill from 356 to 192. This points to primitive
  Gilbert-Peierls workspace reuse as the next performance-critical
  implementation slice.

- An InterPSS `ipss.plugin.3phase` export path for three-phase/OpenDSS feeder
  matrices using the same Matrix Market replay contract.
- A benchmark report template that records command, JVM, matrix summary,
  ordering mode, timings, memory, residual, native KLUSolve/KLUSolveX timing,
  and JKLU/native timing ratio.

## Validation Requirements

- Existing real and complex Maven tests must remain passing.
- Complex results must continue matching CSparseJ for the comparison case.
- Add direct-kernel tests for complex pivot division, row scaling, transpose
  solve, same-pattern refactor, extraction, diagnostics, and
  `KluComplexSparseSet` update/refactor behavior.
- Add sparse-kernel tests for off-diagonal BTF entries, active pivot-column
  traversal, and pruning-sensitive fill behavior as the traversal moves closer
  to the real Gilbert-Peierls kernel.
- Add invalid-input regression tests for malformed CSC arrays and undersized
  RHS/extraction buffers.
- Add extraction regression tests with arrays sized from `Numeric.lnz` and
  `Numeric.unz`.
- Add benchmark coverage comparing JKLU real, JKLU direct complex, CSparseJ, and
  native KLUSolve/KLUSolveX where native libraries are available.
- Add performance profiling runs for generated sparse matrices and captured
  InterPSS large power-system matrices before accepting kernel optimizations.
- Add native KLUSolve/KLUSolveX parity reports for the same large matrix and RHS
  workflow whenever native libraries are available.
- Store benchmark summaries in a durable report or checked-in non-sensitive
  artifact so later optimization cycles can compare against prior baselines.
- Add downstream real power-system validation using local InterPSS plugin
  modules:
  - `/Users/ipssdev/github/ipss-plugin/ipss.plugin.core`
  - `/Users/ipssdev/github/ipss-plugin/ipss.plugin.3phase`
- Use `ipss.plugin.core` to validate transmission/network sparse complex solve
  behavior against existing CSparseJ and KLUSolveX provider tests.
- Use `ipss.plugin.3phase` to validate three-phase distribution and OpenDSS
  feeder sparse complex solve behavior on realistic unbalanced admittance
  matrices.
- Treat InterPSS plugin validation as downstream integration testing because
  `ipss-plugin` builds with Java 21 while JKLU remains Java 8-compatible.

## Real Power-System Validation Cases

- `ipss.plugin.core`: KLUSolveX sparse provider parity, IEEE9/14/39/118/300,
  WECC240, ACTIVSg2000, and ACTIVSg25k-style sparse cases where available.
- `ipss.plugin.3phase`: IEEE9 three-phase, IEEE300 three-phase, IEEE13,
  IEEE123, OpenDSS mini feeders, GridAPPSD DSS feeders, regulator performance,
  QSTS feeder smoke tests, and DistOPF benchmark cases.
- Metrics: factor/refactor/solve time, heap use, sparse LU stored entries,
  residual mismatch, iteration count impact, and downstream power-flow result
  deltas.

Focused downstream commands:

```bash
mvn -pl ipss.plugin.core -Dtest=KlusolveXSparseEqnSolverProviderTest test
mvn -pl ipss.test.plugin.core -Dtest=Dclf_PSSE_ACTIVSg25kBus_Test,DStab_ACTIVSg2000Bus_Test,DStab_IEEE300Bus_Test test
mvn -pl ipss.plugin.3phase -Dtest=IEEE_13BusFeeder_Test,IEEE123Feeder_Dstab_Test,OpenDSSQstsFeederSmokeTest,RegulatorSymbolicUpdatePerformanceTest test
```

## Validation Matrix

| Validation | Status | Command or Evidence |
| --- | --- | --- |
| Real JKLU tests | Passing | `mvn test` |
| Complex hand-solvable systems | Passing | `Zklu_simple` tests |
| Real-as-complex compatibility | Passing | `Zklu_simple` tests |
| Multiple RHS complex solve | Passing | `Zklu_simple` tests |
| Complex refactor API | Passing | `Zklu_simple` tests |
| CSparseJ comparison | Passing | JKLU result matches CSparseJ `DZcs_lusol` within tolerance. |
| Sparse complex LU storage | Passing | `Zklu_simple` checks diagonal complex storage remains sparse. |
| Direct transpose solve | Passing | `Zklu_simple` direct complex transpose solve test. |
| Row scaling | Passing | `Zklu_simple` direct complex row-scaling test. |
| Malformed CSC/RHS validation | Passing | Passing tests cover `klu_z_factor`, `klu_z_solve`, extraction buffer validation, and `klu_z_rgrowth` malformed matrix validation. |
| Sparse L/U extraction sizing | Passing | Tests cover arrays sized from `Numeric.lnz` and `Numeric.unz`. |
| Complex `P/Q` AMD ordering use | Passing | Matrix Market replay shows `P/Q` permutation reduces fill and OpenEI now completes with relative residual below `2e-16`; add a focused unit test for non-natural ordering. |
| Complex `R` BTF block factorization | In progress | Focused BTF/off-diagonal tests pass and singleton BTF fast path is implemented; current code uses `R` block boundaries, but must still match real `Dklu_factor` block-local traversal, `Pblock` composition, and pivot metadata before release. |
| Complex off-diagonal `F` entries | In progress | `Zklu_simple#test_klu_z_btf_off_diagonal_entries_are_solved_and_extracted` verifies initial `Offp`/`Offi`/`Offx` extraction and solve behavior; `Zklu_simple#test_klu_z_direct_btf_off_diagonal_entries_refactor_and_extract` verifies forced-direct same-pattern refactor updates `F` values. Full SuiteSparse/JKLU-style `F` ownership for off-diagonal pivot cases and diagnostics remains release-blocking. |
| Gilbert-Peierls complex kernel | In progress | Large matrices use a left-looking column kernel with direct CSC setup and primitive active workspace. Full `Dklu_kernel` symbolic reach, pruning, pivot traversal, and fallback elimination remain. |
| Large-matrix DFS reach kernel | Passing, default for large matrices | `mvn test` passes with KLU-style `Lpend` DFS pruning; OpenEI replay preserves `960,910` LU entries and `1.907199e-16` relative residual. The reach path is default for `n >= 50000`, uses unsorted cached U by default, and keeps smaller matrices on the heap path. |
| Sparse solve traversal benchmark | In progress | Generated banded complex and Matrix Market replay harnesses added and smoke-tested; InterPSS BUS1824/BUS6384/BUS11856/OpenEI Y-matrix replay passes with relative residuals below `6e-16` for BUS cases and below `2e-16` for OpenEI. |
| Native KLUSolve parity benchmark | In progress | Optional `src/test/native/native_klu_z_bench.cpp` Matrix Market benchmark links to local `libklusolvex.dylib`; latest target sweep reports exact native fill for Ckt24, IEEE8500, ACTIVSg25k, ACTIVSg70K, and OpenEI. The active performance score is now `0.2 * factor + 0.8 * refactor`: Ckt24 `47.6%`, IEEE8500 `52.6%`, ACTIVSg25k `69.6%`, ACTIVSg70K `80.6%`, and OpenEI `87.1%`. Solve remains a correctness/regression gate. |
| Profiling-driven optimization loop | In progress | JKLU can replay Matrix Market matrices with `jklu.benchmark.ordering`, `jklu.benchmark.btf`, `jklu.profile.factor`, and `jklu.profile.analyze` controls; native KLUSolveX comparison is available locally; current next bottlenecks are Java AMD analyze time, factor-side symbolic reach/pruning, and `Pblock`/pivot metadata parity. |
| InterPSS `ipss.plugin.core` large-system validation | Planned | Run focused Java 21 downstream tests after JKLU `2.0.0` candidate artifact is installed. |
| InterPSS `ipss.plugin.3phase` feeder validation | Planned | Run three-phase/OpenDSS feeder tests after JKLU `2.0.0` candidate artifact is installed. |

## Release Checklist

- Complete direct-kernel implementation on `feature/direct-complex-kernel`.
- Confirm `pom.xml` version is `2.0.0` on the direct-kernel branch.
- Remove any remaining block-real backend assumptions from `tdcomplex`.
- Apply AMD `P/Q` symbolic ordering to complex factorization and solve.
- Complete true `Symbolic.R` BTF block factorization parity checks against the
  real `Dklu_factor` block workflow.
- Store, solve, extract, diagnose, and refactor off-diagonal `F` entries in
  SuiteSparse/JKLU style.
- Replace row-based `TreeMap` elimination on the large-matrix path with a
  complex Gilbert-Peierls kernel that matches real `Dklu_kernel` symbolic reach
  and pruning.
- Replace dense-range complex solve/extract/diagnostic scans with sparse
  stored-entry traversal.
- Harden public complex APIs against malformed CSC, RHS, and extraction buffers.
- Implement true same-pattern numeric reuse for `klu_z_refactor`.
- Run the profiling-driven performance loop on generated and InterPSS-derived
  large sparse matrices before finalizing kernel optimizations.
- Confirm JKLU direct complex performance has moved materially closer to native
  KLUSolve/KLUSolveX on the same large InterPSS matrices, with benchmark reports
  showing cold, warm refactor, and solve-only ratios.
- Run the full `1.1.0` validation suite plus direct-kernel parity,
  diagnostics, refactor, extraction, transpose-solve, and benchmark tests.
- Install or publish a candidate JKLU `2.0.0` artifact, wire it into
  `ipss-plugin`, and run focused downstream module tests.
- Publish as a major release after downstream InterPSS/OpenDSS integration
  testing confirms acceptable numeric behavior.

## Deferred Work

- Complete non-natural ordering, singleton, and pivot-metadata parity tests for
  `Symbolic.R` BTF block factorization.
- Complete SuiteSparse/JKLU-style `F` ownership tests covering solve, extract,
  diagnostics, and refactor.
- Replace any remaining hot-loop `TreeMap<Integer, Complex>` large-matrix
  fallback with primitive sparse Gilbert-Peierls workspace after BTF correctness
  is complete.
- Add a streaming first-factor input path that traverses original CSC with
  `Symbolic.P/Q` directly, closer to native KLU, while preserving transformed
  direct-column arrays for same-pattern `klu_z_refactor`, extraction, and
  diagnostics parity.
- Optimize or cache transpose solve behavior for repeated `klu_z_tsolve` calls.
- Add a permanent benchmark suite comparing JKLU real, JKLU complex, CSparseJ,
  and native KLUSolve where native libraries are available.
- Add durable benchmark reports for generated sparse matrices and real
  power-system matrices from `ipss.plugin.core` and `ipss.plugin.3phase`,
  including native KLUSolve/KLUSolveX parity ratios where available.
- Add Matrix Market import/export round-trip tests if Matrix Market import is
  added.
- Add release notes and README examples for downstream InterPSS/OpenDSS
  consumers.
