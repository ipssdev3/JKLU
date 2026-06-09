package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

public class Zklu_extract extends Zklu_internal
{
	public static int klu_z_extract(KLU_z_numeric Numeric, KLU_symbolic Symbolic,
			int[] Lp, int[] Li, double[] Lx, int[] Up, int[] Ui, double[] Ux,
			int[] Fp, int[] Fi, double[] Fx, int[] P, int[] Q, double[] Rs,
			int[] R, KLU_common Common)
	{
		if (Common == null)
		{
			return FALSE;
		}
			if (Numeric == null || Symbolic == null)
			{
				Common.status = KLU_INVALID;
				return FALSE;
				}
				if (Zklu_factor.ensureRowStorage(Numeric, Common) == FALSE)
				{
					return FALSE;
				}
				int n = Numeric.n;
			if (!validOutputs(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux, Fp, Fi, Fx,
					P, Q, Rs, R, Common))
			{
				return FALSE;
			}
			if (Lp != null)
			{
				int p = 0;
				for (int col = 0; col < n; col++)
				{
					Lp[col] = p;
					p = appendEntry(Li, Lx, p, col, 1.0, 0.0);
					for (int row = col + 1; row < n; row++)
					{
						int q = position(Numeric, row, col);
						if (q >= 0)
						{
							p = appendEntry(Li, Lx, p, row,
									Zklu_factor.rowReal(Numeric, row, q),
									Zklu_factor.rowImag(Numeric, row, q));
						}
					}
				}
				Lp[n] = p;
			}
			if (Up != null)
		{
			int p = 0;
				for (int col = 0; col < n; col++)
				{
					Up[col] = p;
					for (int row = 0; row <= col; row++)
					{
						int q = position(Numeric, row, col);
						if (q >= 0)
						{
							p = appendEntry(Ui, Ux, p, row,
									Zklu_factor.rowReal(Numeric, row, q),
									Zklu_factor.rowImag(Numeric, row, q));
						}
					}
				}
				Up[n] = p;
			}
			if (Fp != null)
			{
				System.arraycopy(Numeric.Offp, 0, Fp, 0, n + 1);
				if (Fi != null && Numeric.Offi != null)
				{
					System.arraycopy(Numeric.Offi, 0, Fi, 0, Numeric.nzoff);
				}
				if (Fx != null && Numeric.Offx != null)
				{
					System.arraycopy(Numeric.Offx, 0, Fx, 0, 2 * Numeric.nzoff);
				}
			}
		if (P != null && Numeric.Pnum != null)
		{
			System.arraycopy(Numeric.Pnum, 0, P, 0, n);
		}
		if (Q != null && Symbolic.Q != null)
		{
			System.arraycopy(Symbolic.Q, 0, Q, 0, n);
		}
		if (Rs != null && Numeric.Rs != null)
		{
			System.arraycopy(Numeric.Rs, 0, Rs, 0, n);
		}
		if (R != null && Symbolic.R != null)
		{
			System.arraycopy(Symbolic.R, 0, R, 0, Symbolic.nblocks + 1);
		}
			Common.status = KLU_OK;
			return TRUE;
		}

		private static boolean validOutputs(KLU_z_numeric Numeric, KLU_symbolic Symbolic,
				int[] Lp, int[] Li, double[] Lx, int[] Up, int[] Ui, double[] Ux,
				int[] Fp, int[] Fi, double[] Fx, int[] P, int[] Q, double[] Rs,
				int[] R, KLU_common Common)
		{
			int n = Numeric.n;
			if (Lp != null && Lp.length < n + 1)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			if (Up != null && Up.length < n + 1)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			if ((Li != null && Li.length < Numeric.lnz) ||
				(Lx != null && Lx.length < 2 * Numeric.lnz) ||
				(Ui != null && Ui.length < Numeric.unz) ||
				(Ux != null && Ux.length < 2 * Numeric.unz))
			{
				Common.status = KLU_INVALID;
				return false;
			}
			if (Fp != null && Fp.length < n + 1)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			if ((Fi != null && Fi.length < Numeric.nzoff) ||
				(Fx != null && Fx.length < 2 * Numeric.nzoff) ||
				(P != null && P.length < n) ||
				(Q != null && Q.length < n) ||
				(Rs != null && Rs.length < n) ||
				(R != null && Symbolic.R != null && R.length < Symbolic.nblocks + 1))
			{
				Common.status = KLU_INVALID;
				return false;
			}
			return true;
		}

		private static int appendEntry(int[] indices, double[] values, int p,
				int index, double re, double im)
		{
			if (indices != null)
			{
				indices[p] = index;
			}
			if (values != null)
			{
				values[2 * p] = re;
				values[2 * p + 1] = im;
			}
			return p + 1;
		}

		private static int position(KLU_z_numeric Numeric, int row, int col)
		{
			int rowLength = Zklu_factor.rowLength(Numeric, row);
			for (int p = 0; p < rowLength; p++)
			{
				int c = Zklu_factor.rowCol(Numeric, row, p);
				if (c == col)
				{
					return p;
				}
				if (c > col)
				{
					break;
				}
			}
			return EMPTY;
		}
	}
