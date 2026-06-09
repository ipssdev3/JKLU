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

The following results were measured on the `feature/direct-complex-kernel`
branch. The weighted factor/refactor score is:

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

This sweep exercises the new direct Java complex numeric kernel on complex
power-system admittance matrices generated from InterPSS feeders/cases plus
the TAMU ACTIVSg70K power-network matrix. Settings were `warmups=10`,
`iterations=50` except where noted.

| Case | Matrix size | Input nnz | Iterations | JKLU factor ms | JKLU refactor ms | Native factor ms | Native refactor ms | Weighted factor/refactor | Fill | Residual/conclusion |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| Ckt24 | 18177 x 18177 | 43485 | 50 | 1.534 | 0.415 | 0.467 | 0.263 | 47.6% | 44347 match | Tight residual; fragmented feeder case. |
| IEEE8500 | 14631 x 14631 | 52683 | 50 | 1.561 | 0.543 | 0.659 | 0.326 | 52.6% | 60361 match | Tight residual; fragmented feeder case. |
| ACTIVSg25k | 25000 x 25000 | 85220 | 50 | 5.666 | 1.901 | 3.495 | 1.437 | 69.6% | 186800 match | Correct fill/residual. |
| ACTIVSg70K | 69999 x 69999 | 154313 | 30 | 13.741 | 5.992 | 10.544 | 4.967 | 80.6% | 542881 match | Meets large-case 80% target in this paired run. |
| OpenEI | 78484 x 78484 | 294398 | 30 | 23.778 | 11.720 | 20.458 | 10.265 | 87.1% | 960910 match | Residual 1.907199e-16. |

Retained OpenEI solve-inclusive pair:

| Case | Matrix size | Input nnz | JKLU factor ms | JKLU refactor ms | JKLU solve ms | Native factor ms | Native refactor ms | Native solve ms | Fill | Residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OpenEI complex Y matrix | 78484 x 78484 | 294398 | 25.366 | 12.919 | 1.568 | 20.413 | 10.379 | 1.314 | 960910 match | 1.907199e-16 |

Real-valued Newton-Raphson Jacobian benchmark:

The following comparison uses InterPSS-generated Newton-Raphson Jacobian
matrices for ACTIVSg25k and OpenEI. These are real-valued scalar Jacobian
matrices, so the results validate the direct kernel and sparse traversal on a
real-valued workload rather than the complex Y-matrix numeric path above.
Settings were `warmups=20`, `iterations=100`, `repeats=3`, symbolic reuse
enabled.

| Case | Matrix size | Input nnz | JKLU factor ms | JKLU refactor ms | JKLU solve ms | Native factor ms | Native refactor ms | Native solve ms | Weighted factor/refactor | Solve native/JKLU | Fill | Residual |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| ACTIVSg25k best | 50000 x 50000 | 335468 | 17.064 | 9.133 | 0.909 | 14.622 | 8.144 | 0.905 | 88.06% | 99.62% | 725968 match | 7.26e-15 |
| ACTIVSg25k avg | 50000 x 50000 | 335468 | 17.422 | 9.159 | 0.909 | 14.660 | 8.183 | 0.899 | 87.67% | 98.93% | 725968 match | 7.26e-15 |
| OpenEI best | 156968 x 156968 | 1146882 | 111.537 | 77.089 | 6.190 | 98.468 | 60.403 | 5.264 | 80.99% | 85.04% | 3560808 match | 2.96e-16 |
| OpenEI avg | 156968 x 156968 | 1146882 | 111.653 | 77.339 | 6.174 | 96.930 | 59.819 | 5.240 | 79.86% | 84.86% | 3560808 match | 2.96e-16 |

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
