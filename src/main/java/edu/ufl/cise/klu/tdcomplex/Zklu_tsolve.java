/**
 * Complex transpose solve.
 */

package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

public class Zklu_tsolve extends Zklu_internal
{

	public static int klu_z_tsolve(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
			int d, int nrhs, double[] B, KLU_common Common)
	{
		return klu_z_tsolve(Symbolic, Numeric, d, nrhs, B, 0, Common);
	}

	public static int klu_z_tsolve(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
			int d, int nrhs, double[] B, int B_offset, KLU_common Common)
	{
			if (!Zklu_solve.valid(Symbolic, Numeric, d, nrhs, B, B_offset, Common))
		{
			return FALSE;
		}
		Common.status = KLU_OK;
		int n = Symbolic.n;
		KLU_z_numeric transpose = transposeNumeric(Numeric, Common);
		if (transpose == null)
		{
			return FALSE;
		}
		double[] x = new double[2 * n];
		for (int rhs = 0; rhs < nrhs; rhs++)
		{
			for (int row = 0; row < n; row++)
			{
				int source = B_offset + 2 * (rhs * d + row);
				x[2 * row] = B[source];
				x[2 * row + 1] = B[source + 1];
			}
			Zklu_solve.solveFactored(transpose, x);
			for (int row = 0; row < n; row++)
			{
				int target = B_offset + 2 * (rhs * d + row);
				B[target] = x[2 * row];
				B[target + 1] = x[2 * row + 1];
			}
		}
		return TRUE;
	}

	private static KLU_z_numeric transposeNumeric(KLU_z_numeric Numeric,
			KLU_common Common)
	{
		int n = Numeric.n;
		int nz = Numeric.nz;
		int[] counts = new int[n];
		for (int col = 0; col < n; col++)
		{
			for (int p = Numeric.Ap[col]; p < Numeric.Ap[col + 1]; p++)
			{
				counts[Numeric.Ai[p]]++;
			}
		}
		int[] Ap = new int[n + 1];
		for (int col = 0; col < n; col++)
		{
			Ap[col + 1] = Ap[col] + counts[col];
		}
		int[] next = new int[n];
		System.arraycopy(Ap, 0, next, 0, n);
		int[] Ai = new int[nz];
		double[] Ax = new double[2 * nz];
		for (int col = 0; col < n; col++)
		{
			for (int p = Numeric.Ap[col]; p < Numeric.Ap[col + 1]; p++)
			{
				int row = Numeric.Ai[p];
				int q = next[row]++;
				Ai[q] = col;
				Ax[2 * q] = Numeric.Ax[2 * p];
				Ax[2 * q + 1] = Numeric.Ax[2 * p + 1];
			}
		}
		KLU_z_numeric transpose = new KLU_z_numeric();
		transpose.n = n;
		transpose.nz = nz;
		transpose.Ap = Ap;
		transpose.Ai = Ai;
		transpose.Ax = Ax;
		transpose.Rs = null;
		return Zklu_factor.factorNumeric(transpose, Common) == TRUE ? transpose : null;
	}
}
