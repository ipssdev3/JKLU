package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;

public class Zklu_scale extends Zklu_internal
{
	public static int klu_z_scale(int scale, int n, int[] Ap, int[] Ai,
			double[] Ax, double[] Rs, int[] W, KLU_common Common)
	{
		if (Common == null)
		{
			return FALSE;
		}
		if (Ap == null || Ai == null || Ax == null || Ap.length < n + 1 ||
			Ai.length < Ap[n] || Ax.length < 2 * Ap[n])
		{
			Common.status = KLU_INVALID;
			return FALSE;
		}
		Common.status = KLU_OK;
		if (Rs == null || scale < 0)
		{
			return TRUE;
		}
		for (int row = 0; row < n; row++)
		{
			Rs[row] = 0.0;
			if (W != null)
			{
				W[row] = -1;
			}
		}
		for (int col = 0; col < n; col++)
		{
			for (int p = Ap[col]; p < Ap[col + 1]; p++)
			{
				int row = Ai[p];
				if (row < 0 || row >= n)
				{
					Common.status = KLU_INVALID;
					return FALSE;
				}
				if (W != null)
				{
					if (W[row] == col)
					{
						Common.status = KLU_INVALID;
						return FALSE;
					}
					W[row] = col;
				}
				double a = Zklu_complex.abs(Ax, p);
				if (scale == 1)
				{
					Rs[row] += a;
				}
				else
				{
					Rs[row] = Math.max(Rs[row], a);
				}
			}
		}
		for (int row = 0; row < n; row++)
		{
			if (Rs[row] == 0.0)
			{
				Rs[row] = 1.0;
			}
		}
		return TRUE;
	}
}
