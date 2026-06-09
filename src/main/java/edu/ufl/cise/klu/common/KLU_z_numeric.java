/**
 * KLU: a sparse LU factorization algorithm.
 * Copyright (C) 2004-2009, Timothy A. Davis.
 * Copyright (C) 2011-2012, Richard W. Lincoln.
 * http://www.cise.ufl.edu/research/sparse/klu
 *
 * -------------------------------------------------------------------------
 *
 * KLU is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 */

package edu.ufl.cise.klu.common;

/**
 * Numeric object for complex KLU factorization.
 *
 * Complex values use interleaved double storage: [real0, imag0, real1, imag1].
 */
public class KLU_z_numeric
{

	public int n;
	public int nz;
	public int[] Ap;
	public int[] Ai;
	public int[] ApRef;
	public int[] AiRef;
	public double[] Ax;

	public int lnz;
	public int unz;
	public int nzoff;
	public int nblocks;
	public int[] Pnum;
	public int[] Pinv;
	public int[] pivots;
	public boolean PivotsIdentity;
	public int[] RowPerm;
	public int[] ColPerm;
	public int[] RowInv;
	public int[] ColInv;
	public int[] R;
	public int[] BlockOf;

	public int[] Lip;
	public int[] Uip;
	public int[] Llen;
	public int[] Ulen;
	public int[] LUsize;
	public double[][] LUbx;
	public double[] Udiag;
	public double[] UdiagInv;

	public double[] Rs;
	public int[] Offp;
	public int[] Offi;
	public double[] Offx;

	public double[] Work;
	public double[] Xwork;

	public int[][] LUrowCols;
	public double[][] LUrowValues;
	public boolean LUrowsCurrent;
	public int[][] LrowCols;
	public double[][] LrowValues;
	public int[][] UrowCols;
	public double[][] UrowValues;
	public int[][] DirectColumnRows;
	public double[][] DirectColumnReal;
	public double[][] DirectColumnImag;
	public int[] DirectColumnNext;
	public int[] DirectColumnP;
	public int[] DirectColumnI;
	public double[] DirectColumnX;
	public int DirectHasBlockDiagonalCandidates;
	public Object DirectKernelWork;
	public int[][] DirectLRows;
	public double[][] DirectLReal;
	public double[][] DirectLImag;
	public int[][] DirectURows;
	public double[][] DirectUReal;
	public double[][] DirectUImag;
	public int[] DirectLp;
	public int[] DirectLi;
	public double[] DirectLx;
	public int[] DirectUp;
	public int[] DirectUi;
	public double[] DirectUx;

}
