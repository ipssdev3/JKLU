JKLU: A sparse LU factorization algorithm suited to circuit simulation.
Copyright (C) 2006-2012, Timothy A. Davis.
Copyright (C) 2011-2012, Richard Lincoln.
http://www.cise.ufl.edu/research/sparse/klu/

------------------------------------------------------------------------------

Project home: https://github.com/ipssdev3/JKLU

JKLU 2.0.0 targets Java 8-compatible pure-Java complex sparse solve support
with a direct Java complex numeric kernel while preserving the existing
real-valued API. Complex values use interleaved double arrays in the
SuiteSparse/CSparseJ style:

    [re0, im0, re1, im1, ...]

The complex implementation provides low-level klu_z_* entry points backed by
KLU symbolic analysis, AMD/BTF ordering, same-pattern refactor support, sparse
L/U extraction, diagnostics, and a high-level
edu.ufl.cise.klu.solver.KluComplexSparseSet facade for circuit-style matrix
stamping, primitive admittance stamping, matrix extraction/export helpers, and
island detection.

The complex support is implemented in Java and does not bundle or require
native SuiteSparse, native KLUSolve, or KLUSolveX. Some high-level facade
features are inspired by KLUSolve/KLUSolveX behavior for sparse circuit
matrices, but no native KLUSolve source or binary dependency is included.
CSparseJ is used only as a test/comparison dependency for complex solve
validation. InterPSS and KLUSolve/KLUSolveX comparisons are used for benchmark
and parity validation. KLUSolve/KLUSolveX source reference:
https://github.com/dss-extensions/klusolve

JKLU 2.0.0 direct complex kernel benchmark snapshot:

The following results were refreshed on `master` with JKLU `2.0.0` and local
AMDJ/BTFJ `1.0.2-SNAPSHOT` dependencies. The weighted factor/refactor score is:

    0.2 * factor + 0.8 * refactor

Benchmark environment:

- Hardware: Apple Mac mini, Apple M4 Pro, 12 cores (8 performance, 4
  efficiency), 24 GB memory, arm64.
- OS: macOS 26.3.
- Java: Eclipse Temurin OpenJDK 21.0.10 LTS.
- Native comparison: local native KLUSolve benchmark binary using the same
  Matrix Market inputs, built from the KLUSolve/KLUSolveX source project:
  https://github.com/dss-extensions/klusolve

Complex Y-matrix target sweep:

This sweep exercises the direct Java complex numeric kernel on complex
power-system admittance matrices generated from InterPSS feeders/cases. Settings
were `warmups=10`, `iterations=50`, `repeats=3` except where noted. Rows report
the best weighted factor/refactor repeat.

| Case | Matrix size | Input nnz | Iterations | JKLU factor ms | JKLU refactor ms | Native factor ms | Native refactor ms | Weighted factor/refactor | Fill | Residual/conclusion |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Ckt7 | 3768 x 3768 | 18638 | 50 | 0.620 | 0.233 | 0.195 | 0.103 | 39.3% | 18756 match | Tight residual; tiny fragmented feeder case. |
| Ckt24 | 18177 x 18177 | 43485 | 50 | 1.320 | 0.378 | 0.432 | 0.243 | 49.5% | 44347 match | Tight residual; fragmented feeder case. |
| IEEE8500 | 14631 x 14631 | 52683 | 50 | 1.554 | 0.477 | 0.734 | 0.365 | 63.4% | 60361 match | Tight residual; fragmented feeder case. |
| ACTIVSg25k | 25000 x 25000 | 85220 | 50 | 4.963 | 1.926 | 3.601 | 1.545 | 77.2% | 186800 match | Correct fill/residual. |
| ACTIVSg70K | 69999 x 69999 | 154313 | 30 | 13.741 | 5.992 | 10.544 | 4.967 | 80.6% | 542881 match | Retained from earlier run; local matrix artifact not present in this refresh. |
| OpenEI | 78484 x 78484 | 294398 | 30 | 26.299 | 12.270 | 23.246 | 11.397 | 91.3% | 960910 match | Residual 1.907199e-16. |

Change versus the previous README snapshot:

| Case | Previous weighted factor/refactor | Refreshed weighted factor/refactor | Change |
|---|---:|---:|---:|
| Ckt24 | 47.6% | 49.5% | +1.9 percentage points |
| IEEE8500 | 52.6% | 63.4% | +10.8 percentage points |
| ACTIVSg25k | 69.6% | 77.2% | +7.6 percentage points |
| ACTIVSg70K | 80.6% | 80.6% | Retained |
| OpenEI | 87.1% | 91.3% | +4.2 percentage points |

OpenEI solve-inclusive pair from the same refresh:

| Case | Matrix size | Input nnz | JKLU factor ms | JKLU refactor ms | JKLU solve ms | Native factor ms | Native refactor ms | Native solve ms | Fill | Residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OpenEI complex Y matrix | 78484 x 78484 | 294398 | 26.299 | 12.270 | 1.684 | 23.246 | 11.397 | 1.550 | 960910 match | 1.907199e-16 |

Real-valued Newton-Raphson Jacobian benchmark:

The following comparison uses InterPSS-generated Newton-Raphson Jacobian
matrices for ACTIVSg25k and OpenEI. These are real-valued scalar Jacobian
matrices, so the results validate the direct kernel and sparse traversal on a
real-valued workload rather than the complex Y-matrix numeric path above.
Settings were `warmups=20`, `iterations=100`, `repeats=3`, symbolic reuse
enabled.

| Case | Matrix size | Input nnz | JKLU factor ms | JKLU refactor ms | JKLU solve ms | Native factor ms | Native refactor ms | Native solve ms | Weighted factor/refactor | Solve native/JKLU | Fill | Residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ACTIVSg25k best | 50000 x 50000 | 335468 | 18.268 | 9.405 | 0.934 | 15.641 | 8.672 | 0.932 | 90.05% | 99.81% | 725968 match | 7.26e-15 |
| ACTIVSg25k avg | 50000 x 50000 | 335468 | 18.441 | 9.444 | 0.948 | 15.477 | 8.574 | 0.926 | 88.55% | 97.75% | 725968 match | 7.26e-15 |
| OpenEI best | 156968 x 156968 | 1146882 | 112.716 | 77.483 | 6.191 | 99.317 | 60.681 | 5.105 | 80.93% | 82.46% | 3560808 match | 2.96e-16 |
| OpenEI avg | 156968 x 156968 | 1146882 | 113.518 | 77.714 | 6.269 | 98.622 | 60.625 | 5.263 | 80.38% | 83.96% | 3560808 match | 2.96e-16 |

Change versus the previous README snapshot:

| Case | Previous weighted factor/refactor | Refreshed weighted factor/refactor | Change |
|---|---:|---:|---:|
| ACTIVSg25k best | 88.06% | 90.05% | +1.99 percentage points |
| ACTIVSg25k avg | 87.67% | 88.55% | +0.88 percentage points |
| OpenEI best | 80.99% | 80.93% | -0.06 percentage points |
| OpenEI avg | 79.86% | 80.38% | +0.52 percentage points |

These benchmarks compare pure-Java JKLU against a local native KLUSolve
benchmark binary on the same Matrix Market inputs. They are a performance
snapshot, not a native dependency requirement.

These additions do not change JKLU's license. JKLU remains distributed under
the GNU Lesser General Public License version 2.1 or later, as described below
and in the LICENSE file.

------------------------------------------------------------------------------

JKLU is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

JKLU is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301
