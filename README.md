JKLU: A sparse LU factorization algorithm suited to circuit simulation.
Copyright (C) 2006-2012, Timothy A. Davis.
Copyright (C) 2011-2012, Richard Lincoln.
http://www.cise.ufl.edu/research/sparse/klu/

------------------------------------------------------------------------------

Project home: https://github.com/ipssdev3/JKLU

JKLU 1.1.0 adds Java 8-compatible pure-Java complex sparse solve support while
preserving the existing real-valued API. Complex values use interleaved double
arrays in the SuiteSparse/CSparseJ style:

    [re0, im0, re1, im1, ...]

The complex implementation provides low-level klu_z_* entry points and a
high-level edu.ufl.cise.klu.solver.KluComplexSparseSet facade for circuit-style
matrix stamping, primitive admittance stamping, same-pattern refactor updates,
matrix extraction/export helpers, and island detection.

The complex support is implemented in Java and does not bundle or require
native SuiteSparse, native KLUSolve, or KLUSolveX. Some high-level facade
features are inspired by KLUSolve/KLUSolveX behavior for sparse circuit
matrices, but no native KLUSolve source or binary dependency is included.
CSparseJ is used only as a test/comparison dependency for complex solve
validation.

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
