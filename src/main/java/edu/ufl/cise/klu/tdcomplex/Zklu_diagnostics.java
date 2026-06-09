package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

public class Zklu_diagnostics extends Zklu_internal
{
		public static int klu_z_rgrowth(int[] Ap, int[] Ai, double[] Ax,
				KLU_symbolic Symbolic, KLU_z_numeric Numeric, KLU_common Common)
		{
				if (!valid(Symbolic, Numeric, Common) ||
					!validMatrix(Ap, Ai, Ax, Numeric.n, Common))
				{
					return FALSE;
				}
			if (Zklu_factor.ensureRowStorage(Numeric, Common) == FALSE)
			{
				return FALSE;
			}
			double maxA = matrixNorm(Ap, Ai, Ax, Numeric.n);
		double maxU = upperNorm(Numeric);
		Common.rgrowth = maxA == 0.0 ? 1.0 : maxU / maxA;
		Common.status = KLU_OK;
		return TRUE;
	}

	public static int klu_z_condest(int[] Ap, double[] Ax, KLU_symbolic Symbolic,
			KLU_z_numeric Numeric, KLU_common Common)
	{
		if (!valid(Symbolic, Numeric, Common))
		{
			return FALSE;
		}
		Common.condest = 1.0 / Math.max(klu_z_rcond(Symbolic, Numeric, Common) == TRUE ?
				Common.rcond : 0.0, 1e-300);
		Common.status = KLU_OK;
		return TRUE;
	}

	public static int klu_z_flops(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
			KLU_common Common)
	{
		if (!valid(Symbolic, Numeric, Common))
		{
			return FALSE;
		}
		Common.flops = (2.0 / 3.0) * Numeric.n * Numeric.n * Numeric.n;
		Common.status = KLU_OK;
		return TRUE;
	}

	public static int klu_z_rcond(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
			KLU_common Common)
	{
		if (!valid(Symbolic, Numeric, Common))
		{
			return FALSE;
		}
		double min = Double.POSITIVE_INFINITY;
		double max = 0.0;
		for (int i = 0; i < Numeric.n; i++)
		{
			double a = Zklu_complex.abs(Numeric.Udiag, i);
			min = Math.min(min, a);
			max = Math.max(max, a);
		}
		Common.rcond = max == 0.0 ? 0.0 : min / max;
		Common.status = KLU_OK;
		return TRUE;
	}

	private static boolean valid(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
			KLU_common Common)
	{
		if (Common == null)
		{
			return false;
		}
				if (Symbolic == null || Numeric == null)
				{
					Common.status = KLU_INVALID;
					return false;
				}
				if (Symbolic.n != Numeric.n ||
					Numeric.Udiag == null || Numeric.Udiag.length < 2 * Numeric.n)
				{
					Common.status = KLU_INVALID;
					return false;
				}
					if ((Numeric.DirectLp == null || Numeric.DirectLi == null ||
							Numeric.DirectLx == null || Numeric.DirectUp == null ||
							Numeric.DirectUi == null || Numeric.DirectUx == null) &&
							(Numeric.DirectLRows == null || Numeric.DirectURows == null))
					{
						if (Numeric.LUrowCols == null || Numeric.LUrowValues == null)
						{
						Common.status = KLU_INVALID;
						return false;
					}
					for (int row = 0; row < Numeric.n; row++)
					{
						if (Numeric.LUrowCols[row] == null || Numeric.LUrowValues[row] == null ||
							Numeric.LUrowValues[row].length < 2 * Numeric.LUrowCols[row].length)
						{
							Common.status = KLU_INVALID;
							return false;
						}
					}
				}
				return true;
			}

		private static boolean validMatrix(int[] Ap, int[] Ai, double[] Ax, int n,
				KLU_common Common)
		{
			if (Ap == null || Ai == null || Ax == null || Ap.length < n + 1 ||
				Ap[0] != 0)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			for (int col = 0; col < n; col++)
			{
				if (Ap[col] < 0 || Ap[col] > Ap[col + 1])
				{
					Common.status = KLU_INVALID;
					return false;
				}
			}
			int nz = Ap[n];
			if (nz < 0 || Ai.length < nz || nz > Ax.length / 2)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			for (int p = 0; p < nz; p++)
			{
				if (Ai[p] < 0 || Ai[p] >= n)
				{
					Common.status = KLU_INVALID;
					return false;
				}
			}
			return true;
		}

	private static double matrixNorm(int[] Ap, int[] Ai, double[] Ax, int n)
	{
		if (Ap == null || Ai == null || Ax == null || Ap.length < n + 1)
		{
			return 0.0;
		}
		double max = 0.0;
		for (int col = 0; col < n; col++)
		{
			for (int p = Ap[col]; p < Ap[col + 1]; p++)
			{
				max = Math.max(max, Zklu_complex.abs(Ax, p));
			}
		}
		return max;
	}

		private static double upperNorm(KLU_z_numeric Numeric)
		{
			double max = 0.0;
			for (int row = 0; row < Numeric.n; row++)
			{
				int rowLength = Zklu_factor.rowLength(Numeric, row);
				for (int p = 0; p < rowLength; p++)
				{
					int col = Zklu_factor.rowCol(Numeric, row, p);
					if (col >= row)
					{
						max = Math.max(max, Zklu_complex.abs(
								Zklu_factor.rowReal(Numeric, row, p),
								Zklu_factor.rowImag(Numeric, row, p)));
					}
				}
			}
			return max;
		}
}
