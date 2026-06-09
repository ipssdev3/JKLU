# JKLU 1.1.0 Complex Support Feature Plan

## Release Summary

JKLU `1.1.0` is planned as a Java 8-compatible release that adds pure-Java
complex sparse solve support while preserving the existing real-valued
`tdouble` API. Complex values use interleaved double arrays in the
SuiteSparse/CSparseJ style:

```text
[re0, im0, re1, im1, ...]
```

This release is the additive compatibility release. It exposes low-level
`klu_z_*` entry points and a high-level `KluComplexSparseSet` facade without
adding native SuiteSparse, Eigen, KLUSolve, or KLUSolveX dependencies.

The larger direct Java complex numeric kernel work is tracked separately in:

```text
docs/JKLU-2.0.0-Direct-Complex-Kernel-Feature-Plan.md
```

## Release Scope

- Keep `pom.xml` version at `1.1.0`.
- Keep Maven compiler release at Java 8.
- Preserve the existing real JKLU API and behavior.
- Add complex-number support as new source-compatible APIs.
- Keep the implementation pure Java.
- Do not add native SuiteSparse, KLUSolve, KLUSolveX, or Eigen dependencies.
- Keep direct-kernel performance work out of the `1.1.0` release scope.

## Completed Features

| Status | Feature | Notes |
| --- | --- | --- |
| Done | Release version | `pom.xml` version was `1.1.0` for the 1.1.0 release. |
| Done | Complex API package | Added `edu.ufl.cise.klu.tdcomplex`. |
| Done | Low-level complex factor/solve API | Added `klu_z_factor`, `klu_z_solve`, and `klu_z_refactor`. |
| Done | Complex math helpers | Added offset-based helpers for interleaved complex arrays. |
| Done | CSparseJ comparison test | Added a Maven test comparing JKLU complex solve with CSparseJ `DZcs_lusol`. |
| Done | KLUSolve native smoke comparison | Manually compared a shared complex 2x2 case against local `libklusolvex.dylib`; not added as a Maven dependency. |
| Done | Maven Central release profile | Added `central-release` profile with source jar, javadoc jar, GPG signing, and Central Portal publishing plugin. |
| Done | Maintainer metadata | Updated project metadata to `ipssdev3/JKLU` and added `ipssdev3` as maintainer while preserving `rwl` as original developer. |
| Done | `KluComplexSparseSet` facade | Added a pure-Java high-level complex sparse-set API over `klu_z_*`. |
| Done | Primitive matrix stamping | Added dense interleaved complex primitive admittance stamping with ground-node skip support. |
| Done | Same-pattern update/refactor API | Added increment, zeroise, and explicit `refactor()` support for existing compressed entries. |
| Done | Matrix extraction/export helpers | Added CSC, triplet, and Matrix Market export helpers. |
| Done | Island detection helpers | Added DFS island detection over the compressed nonzero pattern. |

## Feature Tracking

| Priority | Status | Feature | Target |
| --- | --- | --- | --- |
| P0 | Done | Low-level complex solve support | `1.1.0` |
| P0 | Done | Java 8 release compatibility | `1.1.0` |
| P0 | Done | CSparseJ validation | `1.1.0` |
| P0 | Done | Release profile and metadata | `1.1.0` |
| P1 | Done | `KluComplexSparseSet` facade | `1.1.0` |
| P1 | Done | Primitive matrix stamping | `1.1.0` |
| P1 | Done | Incremental update/refactor API | `1.1.0` |
| P2 | Done | Matrix extraction/export helpers | `1.1.0` |
| P2 | Done | Island detection helpers | `1.1.0` |

## KLUSolve-Inspired Enhancements

The following features are inspired by `dss-extensions/klusolve` and remain
pure Java in JKLU:

- Done: add `KluComplexSparseSet` as a high-level facade over low-level
  `klu_z_*` methods.
- Done: support `addMatrixElement`, `setMatrixElement`, `getMatrixElement`,
  `factor`, `refactor`, `solve`, `zero`, `getSize`, `getNNZ`, `getSparseNNZ`,
  and `getSingularCol`.
- Done: use 1-based matrix element methods and skip row/column `0` for
  OpenDSS/KLUSolve-style ground-node compatibility.
- Done: add `addPrimitiveMatrix(int[] nodes, double[] yMatrix)` for dense
  interleaved complex admittance stamping.
- Done: sum duplicate entries during CSC assembly.
- Done: add `incrementMatrixElement`, `zeroiseMatrixElement`, and `refactor`
  for unchanged sparsity patterns.
- Done: add `getCompressedMatrix`, `getTripletMatrix`, and Matrix Market export
  helpers.
- Done: add `findIslands()` using DFS over the compressed matrix nonzero
  pattern.
- Done: add `findDisconnectedSubnetwork()` based on singular-column metadata.

## Validation Matrix

| Validation | Status | Command or Evidence |
| --- | --- | --- |
| Real JKLU tests | Passing | `mvn test` |
| Complex hand-solvable systems | Passing | `Zklu_simple` tests |
| Real-as-complex compatibility | Passing | `Zklu_simple` tests |
| Multiple RHS complex solve | Passing | `Zklu_simple` tests |
| Complex refactor API | Passing | `Zklu_simple` tests |
| CSparseJ comparison | Passing | JKLU result matches CSparseJ `DZcs_lusol` within tolerance. |
| KLUSolve native smoke comparison | Passing | Local `libklusolvex.dylib` returned matching solution for shared complex 2x2 case. |
| Sparse-set facade solve | Passing | `Zklu_sparse_set` tests solve a hand-checkable complex system through `KluComplexSparseSet`. |
| Primitive stamping | Passing | `Zklu_sparse_set` tests dense interleaved complex stamping while skipping ground node `0`. |
| Same-pattern refactor | Passing | `Zklu_sparse_set` tests increment, zeroise, explicit `refactor()`, and solve. |
| Extraction/export/islands | Passing | `Zklu_sparse_set` tests CSC/triplet extraction, Matrix Market headers, and DFS islands. |
| Release artifact generation | Passing | `mvn -Pcentral-release -Dgpg.skip=true clean verify` generated main, sources, and javadoc jars. |

## Release Checklist

- Confirm `pom.xml` version is `1.1.0`.
- Confirm Maven compiler release is Java 8.
- Confirm project URL, SCM, issue tracker, and CI metadata point to
  `https://github.com/ipssdev3/JKLU`.
- Confirm `rwl` remains credited as original developer and `ipssdev3` is listed
  as maintainer.
- Run:

```bash
mvn clean test
```

- Verify release artifacts without signing:

```bash
mvn -Pcentral-release -Dgpg.skip=true clean verify
```

- Configure Central Portal credentials under server id `central` in
  `~/.m2/settings.xml`.
- Configure GPG signing for release publishing.
- Upload to Central Portal without auto-publish:

```bash
mvn -Pcentral-release clean deploy
```

- Review the validated deployment in Central Portal before publishing.

## Deferred To 2.0.0

- Direct Java complex numeric kernel.
- BTF/AMD symbolic ordering integration in complex factorization.
- Sparse stored-entry traversal for complex solve, extraction, and diagnostics.
- Public complex API hardening for malformed CSC, RHS, and extraction buffers.
- Primitive sparse workspace optimization for large matrices.
- Permanent benchmark suite against JKLU real, CSparseJ, and native
  KLUSolve/KLUSolveX where native libraries are available.
