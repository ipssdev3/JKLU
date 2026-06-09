/**
 * Complex KLU solve.
 */

package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

public class Zklu_solve extends Zklu_internal
{

	public static int klu_z_solve(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
			int d, int nrhs, double[] B, KLU_common Common)
	{
		return klu_z_solve(Symbolic, Numeric, d, nrhs, B, 0, Common);
	}

	public static int klu_z_solve(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
			int d, int nrhs, double[] B, int B_offset, KLU_common Common)
	{
			if (!valid(Symbolic, Numeric, d, nrhs, B, B_offset, Common))
		{
			return FALSE;
		}
			Common.status = KLU_OK;
			int n = Symbolic.n;
			double[] x = solveWorkspace(Numeric, n);
			for (int rhs = 0; rhs < nrhs; rhs++)
			{
					loadRhs(Numeric, d, B, B_offset, rhs, x);
				solveFactored(Numeric, x);
				storeRhs(Numeric, n, d, B, B_offset, rhs, x);
		}
			return TRUE;
		}

		private static double[] solveWorkspace(KLU_z_numeric Numeric, int n)
		{
			if (Numeric.Xwork == null || Numeric.Xwork.length < 2 * n)
			{
				Numeric.Xwork = new double[2 * n];
			}
			return Numeric.Xwork;
		}

			static boolean valid(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
				int d, int nrhs, double[] B, int B_offset, KLU_common Common)
	{
		if (Common == null)
		{
			return false;
		}
				if (Symbolic == null || Numeric == null || Numeric.pivots == null ||
					Symbolic.n != Numeric.n || d < Symbolic.n || nrhs < 0 ||
					B == null || B_offset < 0)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			if (nrhs > 0)
			{
				long last = B_offset + 2L * ((long) (nrhs - 1) * d + Symbolic.n - 1) + 1L;
				if (last >= B.length)
				{
					Common.status = KLU_INVALID;
					return false;
				}
			}
				boolean flatAvailable = Numeric.DirectLp != null && Numeric.DirectLi != null &&
						Numeric.DirectLx != null && Numeric.DirectUp != null &&
						Numeric.DirectUi != null && Numeric.DirectUx != null;
				boolean splitColumnAvailable = Numeric.DirectLRows != null &&
						Numeric.DirectURows != null;
				boolean splitRowAvailable = Numeric.LrowCols != null &&
						Numeric.UrowCols != null;
				boolean legacyRowAvailable = Numeric.LUrowCols != null &&
						Numeric.LUrowValues != null;
				boolean fallbackAvailable = splitColumnAvailable || splitRowAvailable || legacyRowAvailable;
				if (!flatAvailable && !fallbackAvailable)
					{
						Common.status = KLU_INVALID;
						return false;
					}
				if (Boolean.getBoolean("jklu.complex.solve.disableFlatColumns") &&
						!fallbackAvailable && !flatAvailable)
					{
						Common.status = KLU_INVALID;
						return false;
					}
				if ((!flatAvailable || (Boolean.getBoolean("jklu.complex.solve.disableFlatColumns") &&
						fallbackAvailable)) &&
						!splitColumnAvailable && !splitRowAvailable && !legacyRowAvailable)
					{
						Common.status = KLU_INVALID;
						return false;
					}
				if (legacyRowAvailable)
				{
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

		private static void loadRhs(KLU_z_numeric Numeric, int d, double[] B,
				int B_offset, int rhs, double[] x)
	{
				int n = Numeric.n;
				for (int row = 0; row < n; row++)
				{
					int oldRow = Numeric.RowPerm == null ? row : Numeric.RowPerm[row];
					int source = B_offset + 2 * (rhs * d + oldRow);
					double scale = Numeric.Rs == null ? 1.0 : Numeric.Rs[oldRow];
					x[2 * row] = B[source] / scale;
					x[2 * row + 1] = B[source + 1] / scale;
				}
	}

		private static void storeRhs(KLU_z_numeric Numeric, int n, int d, double[] B, int B_offset,
				int rhs, double[] x)
	{
			for (int row = 0; row < n; row++)
			{
				int oldCol = Numeric.ColPerm == null ? row : Numeric.ColPerm[row];
				int target = B_offset + 2 * (rhs * d + oldCol);
				B[target] = x[2 * row];
				B[target + 1] = x[2 * row + 1];
			}
	}

			static void solveFactored(KLU_z_numeric Numeric, double[] x)
			{
				boolean flatAvailable = Numeric.DirectLp != null && Numeric.DirectLi != null &&
						Numeric.DirectLx != null && Numeric.DirectUp != null &&
						Numeric.DirectUi != null && Numeric.DirectUx != null;
				boolean fallbackAvailable = (Numeric.DirectLRows != null && Numeric.DirectURows != null) ||
						(Numeric.LrowCols != null && Numeric.UrowCols != null) ||
						(Numeric.LUrowCols != null && Numeric.LUrowValues != null);
				if (flatAvailable &&
						(!fallbackAvailable ||
						!Boolean.getBoolean("jklu.complex.solve.disableFlatColumns")))
				{
						solveFactoredFlatColumns(Numeric, x);
					return;
				}
				if (Numeric.DirectLRows != null && Numeric.DirectURows != null)
				{
					solveFactoredColumns(Numeric, x);
				return;
			}
			if (Numeric.LrowCols != null && Numeric.UrowCols != null)
			{
				solveFactoredSplit(Numeric, x);
				return;
				}
				int n = Numeric.n;
				applyPivots(Numeric, x);

				for (int block = Numeric.nblocks - 1; block >= 0; block--)
			{
				int k1 = Numeric.R[block];
				int k2 = Numeric.R[block + 1];
				for (int row = k1 + 1; row < k2; row++)
				{
					int rowLength = Zklu_factor.rowLength(Numeric, row);
					for (int p = 0; p < rowLength; p++)
				{
					int col = Zklu_factor.rowCol(Numeric, row, p);
						if (col >= row)
						{
							break;
						}
						if (col < k1)
						{
							continue;
						}
						Zklu_complex.multiplySubtract(x, row,
								Zklu_factor.rowReal(Numeric, row, p),
								Zklu_factor.rowImag(Numeric, row, p),
								x[2 * col], x[2 * col + 1]);
					}
				}
				for (int row = k2 - 1; row >= k1; row--)
				{
					int rowLength = Zklu_factor.rowLength(Numeric, row);
					for (int p = 0; p < rowLength; p++)
				{
					int col = Zklu_factor.rowCol(Numeric, row, p);
						if (col <= row)
						{
							continue;
						}
						if (col >= k2)
						{
							break;
						}
						Zklu_complex.multiplySubtract(x, row,
								Zklu_factor.rowReal(Numeric, row, p),
								Zklu_factor.rowImag(Numeric, row, p),
							x[2 * col], x[2 * col + 1]);
				}
			double ar = x[2 * row];
			double ai = x[2 * row + 1];
				divideByUdiag(Numeric, x, row, ar, ai);
			}
				applyOffDiagonal(Numeric, k1, k2, x);
				}
			}

				private static void solveFactoredSplit(KLU_z_numeric Numeric, double[] x)
				{
					applyPivots(Numeric, x);
					for (int block = Numeric.nblocks - 1; block >= 0; block--)
				{
					int k1 = Numeric.R[block];
					int k2 = Numeric.R[block + 1];
					for (int row = k1 + 1; row < k2; row++)
					{
						int[] cols = Numeric.LrowCols[row];
						double[] values = Numeric.LrowValues[row];
						for (int p = 0; p < cols.length; p++)
						{
							int col = cols[p];
							if (col >= k1)
							{
								Zklu_complex.multiplySubtract(x, row,
										values[2 * p], values[2 * p + 1],
										x[2 * col], x[2 * col + 1]);
							}
						}
					}
					for (int row = k2 - 1; row >= k1; row--)
					{
						int[] cols = Numeric.UrowCols[row];
						double[] values = Numeric.UrowValues[row];
						for (int p = cols.length - 1; p >= 0; p--)
						{
							int col = cols[p];
							if (col < k2)
							{
								Zklu_complex.multiplySubtract(x, row,
										values[2 * p], values[2 * p + 1],
										x[2 * col], x[2 * col + 1]);
							}
						}
						double ar = x[2 * row];
						double ai = x[2 * row + 1];
						divideByUdiag(Numeric, x, row, ar, ai);
					}
					applyOffDiagonal(Numeric, k1, k2, x);
				}
				}

							private static void solveFactoredFlatColumns(KLU_z_numeric Numeric, double[] x)
							{
								applyPivots(Numeric, x);
								if (Numeric.nzoff == 0)
								{
									solveFactoredFlatColumnsNoOffDiagonal(Numeric, x);
									return;
								}
								int[] directLp = Numeric.DirectLp;
							int[] directLi = Numeric.DirectLi;
							double[] directLx = Numeric.DirectLx;
							int[] directUp = Numeric.DirectUp;
							int[] directUi = Numeric.DirectUi;
							double[] directUx = Numeric.DirectUx;
							double[] udiag = Numeric.Udiag;
							double[] udiagInv = Numeric.UdiagInv;
								for (int block = Numeric.nblocks - 1; block >= 0; block--)
							{
								int k1 = Numeric.R[block];
								int k2 = Numeric.R[block + 1];
								if (k2 - k1 == 1 && Numeric.nzoff == 0)
								{
									divideFlatColumn(Numeric, x, udiag, udiagInv, k1);
									continue;
								}
								for (int col = k1; col < k2; col++)
								{
								double xr = x[2 * col];
								double xi = x[2 * col + 1];
										for (int p = directLp[col]; p < directLp[col + 1]; p++)
										{
											int xp = 2 * directLi[p];
											double cr = directLx[2 * p];
											double ci = directLx[2 * p + 1];
											x[xp] -= cr * xr - ci * xi;
											x[xp + 1] -= ci * xr + cr * xi;
										}
							}
						for (int col = k2 - 1; col >= k1; col--)
						{
							double ar = x[2 * col];
							double ai = x[2 * col + 1];
								int xpivot = 2 * col;
								double xr;
								double xi;
								if (udiagInv != null && udiagInv.length > xpivot + 1)
								{
									double invr = udiagInv[xpivot];
									double invi = udiagInv[xpivot + 1];
									xr = ar * invr - ai * invi;
									xi = ar * invi + ai * invr;
									x[xpivot] = xr;
								x[xpivot + 1] = xi;
							}
								else
								{
									Zklu_complex.divide(x, col, ar, ai,
											udiag[xpivot], udiag[xpivot + 1]);
									xr = x[xpivot];
									xi = x[xpivot + 1];
								}
										for (int p = directUp[col]; p < directUp[col + 1]; p++)
										{
											int xp = 2 * directUi[p];
											double cr = directUx[2 * p];
											double ci = directUx[2 * p + 1];
											x[xp] -= cr * xr - ci * xi;
											x[xp + 1] -= ci * xr + cr * xi;
										}
						}
						applyOffDiagonal(Numeric, k1, k2, x);
						}
					}

			private static void solveFactoredFlatColumnsNoOffDiagonal(
					KLU_z_numeric Numeric, double[] x)
			{
				int[] directLp = Numeric.DirectLp;
				int[] directLi = Numeric.DirectLi;
				double[] directLx = Numeric.DirectLx;
				int[] directUp = Numeric.DirectUp;
				int[] directUi = Numeric.DirectUi;
				double[] directUx = Numeric.DirectUx;
				double[] udiag = Numeric.Udiag;
				double[] udiagInv = Numeric.UdiagInv;
				for (int block = Numeric.nblocks - 1; block >= 0; block--)
				{
					int k1 = Numeric.R[block];
					int k2 = Numeric.R[block + 1];
					if (k2 - k1 == 1)
					{
						divideFlatColumn(Numeric, x, udiag, udiagInv, k1);
						continue;
					}
					for (int col = k1; col < k2; col++)
					{
						double xr = x[2 * col];
						double xi = x[2 * col + 1];
						for (int p = directLp[col]; p < directLp[col + 1]; p++)
						{
							int xp = 2 * directLi[p];
							double cr = directLx[2 * p];
							double ci = directLx[2 * p + 1];
							x[xp] -= cr * xr - ci * xi;
							x[xp + 1] -= ci * xr + cr * xi;
						}
					}
					for (int col = k2 - 1; col >= k1; col--)
					{
						double ar = x[2 * col];
						double ai = x[2 * col + 1];
						int xpivot = 2 * col;
						double xr;
						double xi;
						if (udiagInv != null && udiagInv.length > xpivot + 1)
						{
							double invr = udiagInv[xpivot];
							double invi = udiagInv[xpivot + 1];
							xr = ar * invr - ai * invi;
							xi = ar * invi + ai * invr;
							x[xpivot] = xr;
							x[xpivot + 1] = xi;
						}
						else
						{
							Zklu_complex.divide(x, col, ar, ai,
									udiag[xpivot], udiag[xpivot + 1]);
							xr = x[xpivot];
							xi = x[xpivot + 1];
						}
						for (int p = directUp[col]; p < directUp[col + 1]; p++)
						{
							int xp = 2 * directUi[p];
							double cr = directUx[2 * p];
							double ci = directUx[2 * p + 1];
							x[xp] -= cr * xr - ci * xi;
							x[xp + 1] -= ci * xr + cr * xi;
						}
					}
				}
			}

			private static void divideFlatColumn(KLU_z_numeric Numeric, double[] x,
					double[] udiag, double[] udiagInv, int col)
			{
				int xpivot = 2 * col;
				double ar = x[xpivot];
				double ai = x[xpivot + 1];
				if (udiagInv != null && udiagInv.length > xpivot + 1)
				{
					double invr = udiagInv[xpivot];
					double invi = udiagInv[xpivot + 1];
					x[xpivot] = ar * invr - ai * invi;
					x[xpivot + 1] = ar * invi + ai * invr;
					return;
				}
				Zklu_complex.divide(x, col, ar, ai,
						udiag[xpivot], udiag[xpivot + 1]);
			}

						private static void solveFactoredColumns(KLU_z_numeric Numeric, double[] x)
						{
						applyPivots(Numeric, x);
					for (int block = Numeric.nblocks - 1; block >= 0; block--)
				{
					int k1 = Numeric.R[block];
					int k2 = Numeric.R[block + 1];
					for (int col = k1; col < k2; col++)
					{
						double xr = x[2 * col];
						double xi = x[2 * col + 1];
							int[] rows = Numeric.DirectLRows[col];
							double[] real = Numeric.DirectLReal[col];
							double[] imag = Numeric.DirectLImag[col];
							for (int p = 0; p < rows.length; p++)
						{
							int row = rows[p];
							int xp = 2 * row;
							double cr = real[p];
							double ci = imag[p];
							x[xp] -= cr * xr - ci * xi;
							x[xp + 1] -= ci * xr + cr * xi;
						}
					}
					for (int col = k2 - 1; col >= k1; col--)
					{
						double ar = x[2 * col];
						double ai = x[2 * col + 1];
						divideByUdiag(Numeric, x, col, ar, ai);
						double xr = x[2 * col];
						double xi = x[2 * col + 1];
							int[] rows = Numeric.DirectURows[col];
							double[] real = Numeric.DirectUReal[col];
							double[] imag = Numeric.DirectUImag[col];
							for (int p = 0; p < rows.length; p++)
						{
							int row = rows[p];
							int xp = 2 * row;
							double cr = real[p];
							double ci = imag[p];
							x[xp] -= cr * xr - ci * xi;
							x[xp + 1] -= ci * xr + cr * xi;
						}
					}
					applyOffDiagonal(Numeric, k1, k2, x);
			}
		}

		private static void divideByUdiag(KLU_z_numeric Numeric, double[] x,
				int row, double ar, double ai)
		{
			int p = 2 * row;
			if (Numeric.UdiagInv != null && Numeric.UdiagInv.length > p + 1)
			{
				double invr = Numeric.UdiagInv[p];
				double invi = Numeric.UdiagInv[p + 1];
				x[p] = ar * invr - ai * invi;
				x[p + 1] = ar * invi + ai * invr;
				return;
			}
			Zklu_complex.divide(x, row, ar, ai,
					Numeric.Udiag[p], Numeric.Udiag[p + 1]);
		}

			private static void applyOffDiagonal(KLU_z_numeric Numeric, int k1, int k2,
				double[] x)
		{
			if (Numeric.nzoff == 0 ||
					Numeric.Offp == null || Numeric.Offi == null || Numeric.Offx == null)
			{
				return;
			}
			for (int col = k1; col < k2; col++)
			{
				double xr = x[2 * col];
				double xi = x[2 * col + 1];
				for (int p = Numeric.Offp[col]; p < Numeric.Offp[col + 1]; p++)
				{
					int row = Numeric.Offi[p];
					Zklu_complex.multiplySubtract(x, row,
							Numeric.Offx[2 * p], Numeric.Offx[2 * p + 1], xr, xi);
				}
			}
				}

		private static void applyPivots(KLU_z_numeric Numeric, double[] x)
		{
			if (Numeric.PivotsIdentity)
			{
				return;
			}
			for (int k = 0; k < Numeric.n; k++)
			{
				int pivot = Numeric.pivots[k];
				if (pivot != k)
				{
					swap(x, k, pivot);
				}
			}
		}

		static void swap(double[] x, int a, int b)
		{
		double tr = x[2 * a];
		double ti = x[2 * a + 1];
		x[2 * a] = x[2 * b];
		x[2 * a + 1] = x[2 * b + 1];
		x[2 * b] = tr;
		x[2 * b + 1] = ti;
	}
}
