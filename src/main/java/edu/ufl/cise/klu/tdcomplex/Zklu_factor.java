/**
 * Complex KLU factorization.
 */

package edu.ufl.cise.klu.tdcomplex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

import static edu.ufl.cise.klu.tdcomplex.Zklu_scale.klu_z_scale;

public class Zklu_factor extends Zklu_internal
{

	private static final double DROP_TOL = 0.0;
	private static final ThreadLocal<ColumnKernelWork> FACTOR_WORK =
			new ThreadLocal<ColumnKernelWork>();

	static boolean validInputs(int[] Ap, int[] Ai, double[] Ax, KLU_symbolic Symbolic,
			KLU_common Common)
	{
		if (Common == null)
		{
			return false;
		}
		if (Ap == null || Ai == null || Ax == null || Symbolic == null ||
			Symbolic.n < 0 || Ap.length < Symbolic.n + 1)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			if (Ap[0] != 0)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			for (int col = 0; col < Symbolic.n; col++)
			{
				if (Ap[col] < 0 || Ap[col] > Ap[col + 1])
				{
					Common.status = KLU_INVALID;
					return false;
				}
			}
			int nz = Ap[Symbolic.n];
			if (nz < 0 || Ai.length < nz || nz > Ax.length / 2)
			{
				Common.status = KLU_INVALID;
				return false;
			}
			for (int p = 0; p < nz; p++)
			{
				if (Ai[p] < 0 || Ai[p] >= Symbolic.n)
				{
					Common.status = KLU_INVALID;
					return false;
				}
			}
			return true;
		}

	static boolean samePattern(int[] Ap, int[] Ai, KLU_z_numeric Numeric)
	{
		if (Numeric == null || Numeric.Ap == null || Numeric.Ai == null ||
			Ap == null || Ai == null || Ap.length < Numeric.n + 1)
		{
			return false;
		}
		if (Ap[Numeric.n] != Numeric.nz)
		{
			return false;
		}
		if ((Ap == Numeric.Ap && Ai == Numeric.Ai) ||
				(Ap == Numeric.ApRef && Ai == Numeric.AiRef))
		{
			return true;
		}
		for (int col = 0; col <= Numeric.n; col++)
		{
			if (Ap[col] != Numeric.Ap[col])
			{
				return false;
			}
		}
		for (int p = 0; p < Numeric.nz; p++)
		{
			if (Ai[p] != Numeric.Ai[p])
			{
				return false;
			}
		}
		return true;
	}

		static int factorNumeric(KLU_z_numeric Numeric, KLU_common Common)
	{
			int n = Numeric.n;
			ensureSymbolicState(Numeric);
			FactorProfile profile = FactorProfile.create(n);
			long directInputStartNs = profile.enabled() ? System.nanoTime() : 0L;
			DirectColumns directColumns = buildDirectColumnKernelInput(Numeric, Common);
			if (profile.enabled())
			{
				profile.directInput(System.nanoTime() - directInputStartNs,
						Numeric.DirectColumnRows != null || Numeric.DirectColumnP != null);
			}
				if (directColumns != null && hasDirectFactors(Numeric) &&
						hasBlockDiagonalCandidates(Numeric, directColumns) &&
						tryDirectRefactor(Numeric, directColumns, profile, Common))
				{
					refreshUdiagInverse(Numeric);
					return TRUE;
				}
				if (directColumns != null && hasBlockDiagonalCandidates(Numeric, directColumns) &&
						useReachKernel(n) &&
						tryReachColumnKernel(Numeric, directColumns, profile, Common))
				{
					refreshUdiagInverse(Numeric);
					return TRUE;
				}
				if (directColumns != null && hasBlockDiagonalCandidates(Numeric, directColumns) &&
						tryColumnKernel(Numeric, directColumns, profile, Common))
				{
					refreshUdiagInverse(Numeric);
					return TRUE;
				}
			long setup0 = System.nanoTime();
			SparseRow[] rows = buildSparseRows(Numeric, Common);
		long setup1 = System.nanoTime();
		if (rows == null)
		{
			return FALSE;
			}
			extractOffDiagonalBlocks(Numeric, rows);
			TreeSet<Integer>[] columnRows = buildColumnRows(rows, n);
		long setup2 = System.nanoTime();
		profile.printSetup(setup1 - setup0, setup2 - setup1, rows);
		int[] pivots = new int[n];
		int[] pnum = new int[n];
		int[] pinv = new int[n];
		for (int i = 0; i < n; i++)
		{
			pnum[i] = i;
			pinv[i] = i;
		}

		Common.status = KLU_OK;
		Common.numerical_rank = EMPTY;
		Common.singular_col = EMPTY;
		Common.noffdiag = 0;

		for (int block = 0; block < Numeric.nblocks; block++)
		{
			int k1 = Numeric.R[block];
			int k2 = Numeric.R[block + 1];
			for (int k = k1; k < k2; k++)
			{
				int pivotRow = k;
				double pivotAbs = 0.0;
				int pivotCandidates = 0;
				for (Integer rowRef : columnRows[k].subSet(Integer.valueOf(k), Integer.valueOf(k2)))
				{
					pivotCandidates++;
					int row = rowRef.intValue();
					double candidate = rows[row].abs(k);
					if (candidate > pivotAbs)
					{
						pivotAbs = candidate;
						pivotRow = row;
					}
				}
				pivots[k] = pivotRow;
				if (pivotRow != k)
				{
					swapRows(rows, columnRows, k, pivotRow);
					int t = pnum[k];
					pnum[k] = pnum[pivotRow];
					pnum[pivotRow] = t;
					Common.noffdiag++;
				}

				Complex pivot = rows[k].get(k);
				if (pivot == null || pivot.isZero())
				{
					Common.status = KLU_SINGULAR;
					if (Common.numerical_rank == EMPTY)
					{
						Common.numerical_rank = k;
						Common.singular_col = k;
					}
					profile.step(k, pivotCandidates, 0, rows[k].size(), 0, rows);
					if (Common.halt_if_singular != 0 || profile.shouldStop(k))
					{
						return FALSE;
					}
					continue;
				}

				List<Integer> rowsWithPivot = new ArrayList<Integer>(
						columnRows[k].subSet(Integer.valueOf(k + 1), Integer.valueOf(k2)));
				int updateCount = 0;
				for (Integer rowRef : rowsWithPivot)
				{
					int row = rowRef.intValue();
					Complex entry = rows[row].get(k);
					if (entry == null || entry.isZero())
					{
						continue;
					}
					Complex multiplier = entry.divide(pivot);
					rows[row].put(k, multiplier.re, multiplier.im);
					for (Map.Entry<Integer, Complex> pivotEntry : rows[k].tailMap(k + 1))
					{
						int col = pivotEntry.getKey().intValue();
						if (col >= k2)
						{
							break;
						}
						updateCount++;
						Complex u = pivotEntry.getValue();
						int change = rows[row].subtractProduct(col, multiplier, u);
						if (change > 0)
						{
							columnRows[col].add(Integer.valueOf(row));
						}
						else if (change < 0)
						{
							columnRows[col].remove(Integer.valueOf(row));
						}
					}
				}
				profile.step(k, pivotCandidates, rowsWithPivot.size(), rows[k].size(),
						updateCount, rows);
				if (profile.shouldStop(k))
				{
					Common.status = KLU_OK;
					return FALSE;
				}
			}
		}

		Numeric.pivots = pivots;
		Numeric.PivotsIdentity = pivotsIdentity(pivots);
		Numeric.Pnum = new int[n];
		Numeric.Pinv = new int[n];
		for (int i = 0; i < n; i++)
		{
			int oldRow = Numeric.RowPerm == null ? pnum[i] : Numeric.RowPerm[pnum[i]];
			Numeric.Pnum[i] = oldRow;
			Numeric.Pinv[oldRow] = i;
		}
		Numeric.Udiag = new double[2 * n];
		for (int k = 0; k < n; k++)
		{
			Complex pivot = rows[k].get(k);
			Numeric.Udiag[2 * k] = pivot == null ? 0.0 : pivot.re;
			Numeric.Udiag[2 * k + 1] = pivot == null ? 0.0 : pivot.im;
		}
		refreshUdiagInverse(Numeric);
		populateSparsePatternState(Numeric, rows);
		profile.finish(rows);
			return Common.status == KLU_OK || Common.status == KLU_SINGULAR ? TRUE : FALSE;
		}

		static int ensureRowStorage(KLU_z_numeric Numeric, KLU_common Common)
		{
			if (Numeric == null)
			{
				if (Common != null)
				{
					Common.status = KLU_INVALID;
				}
				return FALSE;
			}
			if (Numeric.LUrowsCurrent)
			{
				return TRUE;
			}
			if (((Numeric.DirectLp == null || Numeric.DirectLi == null ||
					Numeric.DirectLx == null || Numeric.DirectUp == null ||
					Numeric.DirectUi == null || Numeric.DirectUx == null) &&
					(Numeric.DirectLRows == null || Numeric.DirectURows == null)) ||
					Numeric.Udiag == null)
			{
				if (Common != null)
				{
					Common.status = KLU_INVALID;
				}
				return FALSE;
			}
			ColumnKernelWork work = columnKernelWork(Numeric, Numeric.n);
			RowBuilder[] rows = work.luRows;
			for (int row = 0; row < Numeric.n; row++)
			{
				rows[row].clear();
			}
			for (int col = 0; col < Numeric.n; col++)
			{
				if (Numeric.DirectUp != null && Numeric.DirectUi != null &&
						Numeric.DirectUx != null)
				{
					for (int p = Numeric.DirectUp[col]; p < Numeric.DirectUp[col + 1]; p++)
					{
						rows[Numeric.DirectUi[p]].add(col,
								Numeric.DirectUx[2 * p], Numeric.DirectUx[2 * p + 1]);
					}
				}
				else
					{
						int[] uRows = Numeric.DirectURows[col];
						double[] uReal = Numeric.DirectUReal[col];
						double[] uImag = Numeric.DirectUImag[col];
						for (int p = 0; p < uRows.length; p++)
						{
							rows[uRows[p]].add(col, uReal[p], uImag[p]);
						}
				}
				rows[col].add(col, Numeric.Udiag[2 * col], Numeric.Udiag[2 * col + 1]);
				if (Numeric.DirectLp != null && Numeric.DirectLi != null &&
						Numeric.DirectLx != null)
				{
					for (int p = Numeric.DirectLp[col]; p < Numeric.DirectLp[col + 1]; p++)
					{
						rows[Numeric.DirectLi[p]].add(col,
								Numeric.DirectLx[2 * p], Numeric.DirectLx[2 * p + 1]);
					}
				}
				else
					{
						int[] lRows = Numeric.DirectLRows[col];
						double[] lReal = Numeric.DirectLReal[col];
						double[] lImag = Numeric.DirectLImag[col];
						for (int p = 0; p < lRows.length; p++)
						{
							rows[lRows[p]].add(col, lReal[p], lImag[p]);
						}
				}
			}
			populateSparsePatternState(Numeric, rows);
			if (Common != null)
			{
				Common.status = KLU_OK;
			}
			return TRUE;
		}

		private static boolean tryReachColumnKernel(KLU_z_numeric Numeric,
				DirectColumns initialColumns, FactorProfile profile, KLU_common Common)
		{
			int n = Numeric.n;
			ColumnKernelWork work = factorColumnKernelWork(n);
			FlatColumnStore lStore = work.lStore;
			FlatColumnStore uStore = work.uStore;
			RowBuilder[] luRows = work.luRows;
			boolean materializeRows = Boolean.getBoolean(
					"jklu.complex.factor.materializeRows");
			boolean collectRows = materializeRows || profile.enabled();
			work.clearBuilders(collectRows);
			lStore.clear();
			uStore.clear();
				int expectedFactorEntries = expectedFactorEntries(n, Numeric.nz);
			lStore.ensureCapacity(expectedFactorEntries);
			uStore.ensureCapacity(expectedFactorEntries);
			double[] xr = work.xr;
			double[] xi = work.xi;
			int[] mark = work.mark;
			int[] lpend = work.lpend;
			IntList touched = work.touched;
			touched.ensureCapacity(n);
			int[] pivots = work.pivots;
			boolean profilePhases = profile.enabled();
			boolean reachUnsortedU = useReachUnsortedU(n);
			int reachSortInsertionLimit = reachInsertionSortLimit();
			double pivotTol2 = Common.tol * Common.tol;
			Numeric.Udiag = new double[2 * n];
			for (int i = 0; i < n; i++)
			{
				lpend[i] = EMPTY;
			}

				for (int block = 0; block < Numeric.nblocks; block++)
				{
					int k1 = Numeric.R[block];
					int k2 = Numeric.R[block + 1];
					for (int k = k1; k < k2; k++)
					{
						int inputStart = initialColumns.start(k);
						int inputEnd = initialColumns.end(k);
						if (isDirectDiagonalColumn(initialColumns, k, inputStart, inputEnd) &&
								!isZero(initialColumns.real(k, inputStart),
										initialColumns.imag(k, inputStart)))
						{
							recordDirectDiagonalColumn(Numeric, lStore, uStore, luRows,
									collectRows, profile, k, initialColumns.real(k, inputStart),
									initialColumns.imag(k, inputStart));
							continue;
						}
						int generation = work.nextGeneration();
						int reachGeneration = work.nextReachGeneration();
						int[] stack = work.reachStack;
					int[] reachMark = work.reachMark;
					int[] reachPos = work.reachPos;
					int[] reachRows = work.reachRows;
					work.reachCount = 0;
					int top = n;

						long phase0 = profilePhases ? System.nanoTime() : 0L;
							if (initialColumns.p != null)
					{
						int[] inputRows = initialColumns.i;
						double[] inputValues = initialColumns.x;
						for (int p = inputStart; p < inputEnd; p++)
						{
							int row = inputRows[p];
							if (row >= k1 && row < k2)
							{
								int xp = 2 * p;
								addReachWork(row, inputValues[xp],
										inputValues[xp + 1], xr, xi, mark,
										generation, touched);
								if (reachMark[row] != reachGeneration)
								{
									if (row < k)
									{
										top = reachDfs(row, k, lStore, reachMark,
												reachGeneration, stack, reachPos, top,
												work, lpend);
									}
									else
									{
										reachMark[row] = reachGeneration;
										reachRows[work.reachCount++] = row;
									}
								}
							}
						}
						}
							else
							{
								for (int p = inputStart; p < inputEnd; p++)
								{
									int row = initialColumns.row(k, p);
									if (row >= k1 && row < k2)
									{
										addReachWork(row, initialColumns.real(k, p),
												initialColumns.imag(k, p), xr, xi, mark,
												generation, touched);
									if (reachMark[row] != reachGeneration)
								{
								if (row < k)
								{
									top = reachDfs(row, k, lStore, reachMark,
											reachGeneration, stack, reachPos, top,
											work, lpend);
								}
								else
								{
									reachMark[row] = reachGeneration;
									reachRows[work.reachCount++] = row;
								}
								}
							}
						}
					}
					long scatterNs = profilePhases ? System.nanoTime() - phase0 : 0L;

					uStore.startColumn(k);
					int updateCount = 0;
					phase0 = profilePhases ? System.nanoTime() : 0L;
					for (int s = top; s < n; s++)
					{
						int row = stack[s];
						double ur = xr[row];
						double ui = xi[row];
							xr[row] = 0.0;
							xi[row] = 0.0;
							uStore.addFast(row, ur, ui);
						if (ur == 0.0 && ui == 0.0)
						{
							continue;
						}
								int lStart = lStore.p[row];
								int lSize = lStore.len[row];
								int[] targets = lStore.index;
								double[] lReal = lStore.real;
								double[] lImag = lStore.imag;
							if (profilePhases)
							{
								profile.updateColumnLength(lSize);
								for (int p = 0; p < lSize; p++)
								{
									int q = lStart + p;
									int target = targets[q];
									double lr = lReal[q];
									double li = lImag[q];
									profile.updateShape(ur, ui, lr, li);
									profile.updateTarget(target < k);
									updateCount++;
									xr[target] -= lr * ur - li * ui;
									xi[target] -= li * ur + lr * ui;
								}
							}
							else
							{
								int q = lStart;
								int end = lStart + lSize;
								int unrolledEnd = end - ((end - q) & 3);
								for (; q < unrolledEnd; q += 4)
								{
									int target0 = targets[q];
									double lr0 = lReal[q];
									double li0 = lImag[q];
									xr[target0] -= lr0 * ur - li0 * ui;
									xi[target0] -= li0 * ur + lr0 * ui;

									int q1 = q + 1;
									int target1 = targets[q1];
									double lr1 = lReal[q1];
									double li1 = lImag[q1];
									xr[target1] -= lr1 * ur - li1 * ui;
									xi[target1] -= li1 * ur + lr1 * ui;

									int q2 = q + 2;
									int target2 = targets[q2];
									double lr2 = lReal[q2];
									double li2 = lImag[q2];
									xr[target2] -= lr2 * ur - li2 * ui;
									xi[target2] -= li2 * ur + lr2 * ui;

									int q3 = q + 3;
									int target3 = targets[q3];
									double lr3 = lReal[q3];
									double li3 = lImag[q3];
									xr[target3] -= lr3 * ur - li3 * ui;
									xi[target3] -= li3 * ur + lr3 * ui;
								}
								for (; q < end; q++)
								{
									int target = targets[q];
									double lr = lReal[q];
									double li = lImag[q];
									xr[target] -= lr * ur - li * ui;
									xi[target] -= li * ur + lr * ui;
								}
							}
					}
						if (!reachUnsortedU)
					{
						sortFlatColumn(uStore, k, reachSortInsertionLimit);
					}
					long updateNs = profilePhases ? System.nanoTime() - phase0 : 0L;

					pivots[k] = k;
					double pivotR = xr[k];
					double pivotI = xi[k];
					Numeric.Udiag[2 * k] = pivotR;
					Numeric.Udiag[2 * k + 1] = pivotI;
					if (collectRows)
					{
						luRows[k].add(k, pivotR, pivotI);
						int uStart = uStore.p[k];
						int uLen = uStore.len[k];
						for (int p = 0; p < uLen; p++)
						{
							int q = uStart + p;
							luRows[uStore.index[q]].add(k, uStore.real[q], uStore.imag[q]);
						}
					}

					lStore.startColumn(k);
					phase0 = profilePhases ? System.nanoTime() : 0L;
					double pivotAbs2 = pivotR * pivotR + pivotI * pivotI;
					if (pivotAbs2 == 0.0)
					{
						clearReachWork(xr, xi, touched, stack, top, n, reachRows,
								work.reachCount);
						return false;
					}
					double maxAbs2 = pivotAbs2;
					double invDenom = 1.0 / pivotAbs2;
					for (int p = 0; p < work.reachCount; p++)
					{
						int row = reachRows[p];
						if (row > k && row < k2)
						{
							double candidateAbs2 = xr[row] * xr[row] + xi[row] * xi[row];
							if (candidateAbs2 > maxAbs2)
							{
								maxAbs2 = candidateAbs2;
							}
						}
						if (row <= k || row >= k2)
						{
							continue;
						}
							double lr = (xr[row] * pivotR + xi[row] * pivotI) * invDenom;
							double li = (xi[row] * pivotR - xr[row] * pivotI) * invDenom;
							lStore.addFast(row, lr, li);
							xr[row] = 0.0;
							xi[row] = 0.0;
							if (collectRows)
							{
								luRows[row].add(k, lr, li);
							}
						}
						long scaleNs = profilePhases ? System.nanoTime() - phase0 : 0L;
						phase0 = profilePhases ? System.nanoTime() : 0L;
						if (pivotAbs2 < pivotTol2 * maxAbs2)
						{
							clearReachWork(xr, xi, touched, stack, top, n, reachRows,
									work.reachCount);
							return false;
						}
						long pivotNs = profilePhases ? System.nanoTime() - phase0 : 0L;
						pruneReachColumns(lStore, uStore, k, lpend);
					if (profilePhases)
					{
						profile.step(k, work.reachCount, lStore.len[k], uStore.len[k],
								updateCount, luRows);
					}
					phase0 = profilePhases ? System.nanoTime() : 0L;
					clearReachWork(xr, xi, touched, stack, top, n, reachRows,
							work.reachCount);
					if (profilePhases)
					{
						profile.directPhases(scatterNs, updateNs, pivotNs, scaleNs,
								System.nanoTime() - phase0);
					}
					if (profilePhases && profile.shouldStop(k))
					{
						return false;
					}
				}
			}

			Numeric.pivots = pivots;
			Numeric.PivotsIdentity = true;
			Numeric.Pnum = new int[n];
			Numeric.Pinv = new int[n];
			for (int i = 0; i < n; i++)
			{
				int oldRow = Numeric.RowPerm == null ? i : Numeric.RowPerm[i];
				Numeric.Pnum[i] = oldRow;
				Numeric.Pinv[oldRow] = i;
			}
			Common.status = KLU_OK;
			long finalizeStartNs = profilePhases ? System.nanoTime() : 0L;
				if (reachUnsortedU && useReachSortCache(n))
				{
					sortFlatColumns(uStore, n);
				}
				lStore.finish(n);
				uStore.finish(n);
				cacheDirectColumnFactors(Numeric, lStore, uStore);
			if (materializeRows)
			{
				populateSparsePatternState(Numeric, luRows);
			}
			else
			{
				setDirectCounts(Numeric);
				Numeric.LUrowsCurrent = false;
			}
			if (profilePhases)
			{
				profile.directFinalize(System.nanoTime() - finalizeStartNs);
				profile.directActive(0L, 0L);
			}
			if (collectRows)
			{
				profile.finish(luRows);
			}
			return true;
		}

		private static int reachDfs(int start, int limit, FlatColumnStore lStore,
				int[] reachMark, int generation, int[] stack, int[] reachPos,
				int top, ColumnKernelWork work, int[] lpend)
		{
			int head = 0;
			stack[0] = start;
			while (head >= 0)
			{
				int row = stack[head];
				if (reachMark[row] != generation)
				{
					reachMark[row] = generation;
					reachPos[head] = lpend[row] == EMPTY ? lStore.len[row] : lpend[row];
				}
				int startPos = lStore.p[row];
				int pos;
				for (pos = --reachPos[head]; pos >= 0; pos--)
				{
					int target = lStore.index[startPos + pos];
					if (reachMark[target] != generation)
					{
						if (target < limit)
						{
							reachPos[head] = pos;
							stack[++head] = target;
							break;
						}
						reachMark[target] = generation;
						work.reachRows[work.reachCount++] = target;
					}
				}
				if (pos == -1)
				{
					head--;
					stack[--top] = row;
				}
			}
			return top;
		}

		private static void clearReachWork(double[] xr, double[] xi, IntList touched,
				int[] stack, int top, int n, int[] reachRows, int reachCount)
		{
			for (int p = 0; p < touched.size; p++)
			{
				int row = touched.values[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			for (int p = top; p < n; p++)
			{
				int row = stack[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			for (int p = 0; p < reachCount; p++)
			{
				int row = reachRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			touched.clear();
		}

		private static void pruneReachColumns(FlatColumnStore lStore,
				FlatColumnStore uStore, int k, int[] lpend)
		{
			int uStart = uStore.p[k];
			int uEnd = uStart + uStore.len[k];
			for (int up = uStart; up < uEnd; up++)
			{
				int col = uStore.index[up];
				if (col >= k || lpend[col] != EMPTY)
				{
					continue;
				}
				int lStart = lStore.p[col];
				int lEnd = lStart + lStore.len[col];
				boolean containsPivot = false;
				for (int p = lStart; p < lEnd; p++)
				{
					if (lStore.index[p] == k)
					{
						containsPivot = true;
						break;
					}
				}
				if (!containsPivot)
				{
					continue;
				}
				int head = lStart;
				int tail = lEnd;
				while (head < tail)
				{
					int row = lStore.index[head];
					if (row <= k)
					{
						head++;
					}
					else
					{
						tail--;
						swapFlatEntries(lStore, head, tail);
					}
				}
				lpend[col] = tail - lStart;
			}
		}

		private static void sortFlatColumn(FlatColumnStore store, int col)
		{
			sortFlatColumn(store, col, reachInsertionSortLimit());
		}

		private static int reachInsertionSortLimit()
		{
			return Integer.getInteger("jklu.complex.factor.reachInsertionSortLimit", 16)
					.intValue();
		}

		private static void sortFlatColumn(FlatColumnStore store, int col,
				int insertionLimit)
		{
			int start = store.p[col];
			int end = start + store.len[col];
			if (end - start <= 1)
			{
				return;
			}
			if (end - start > insertionLimit)
			{
				quickSortFlatColumn(store, start, end - 1);
				return;
			}
			for (int p = start + 1; p < end; p++)
			{
				int row = store.index[p];
				double re = store.real[p];
				double im = store.imag[p];
				int q = p - 1;
				while (q >= start && store.index[q] > row)
				{
					store.index[q + 1] = store.index[q];
					store.real[q + 1] = store.real[q];
					store.imag[q + 1] = store.imag[q];
					q--;
				}
				store.index[q + 1] = row;
				store.real[q + 1] = re;
				store.imag[q + 1] = im;
			}
		}

		private static void sortFlatColumns(FlatColumnStore store, int n)
		{
			int insertionLimit = reachInsertionSortLimit();
			for (int col = 0; col < n; col++)
			{
				sortFlatColumn(store, col, insertionLimit);
			}
		}

		private static void quickSortFlatColumn(FlatColumnStore store, int left,
				int right)
		{
			int i = left;
			int j = right;
			int pivot = store.index[(left + right) >>> 1];
			while (i <= j)
			{
				while (store.index[i] < pivot)
				{
					i++;
				}
				while (store.index[j] > pivot)
				{
					j--;
				}
				if (i <= j)
				{
					swapFlatEntries(store, i, j);
					i++;
					j--;
				}
			}
			if (left < j)
			{
				quickSortFlatColumn(store, left, j);
			}
			if (i < right)
			{
				quickSortFlatColumn(store, i, right);
			}
		}

		private static void swapFlatEntries(FlatColumnStore store, int a, int b)
		{
			if (a == b)
			{
				return;
			}
			int row = store.index[a];
			double re = store.real[a];
			double im = store.imag[a];
			store.index[a] = store.index[b];
			store.real[a] = store.real[b];
			store.imag[a] = store.imag[b];
			store.index[b] = row;
			store.real[b] = re;
			store.imag[b] = im;
		}

		private static boolean tryColumnKernel(KLU_z_numeric Numeric, DirectColumns initialColumns,
				FactorProfile profile, KLU_common Common)
		{
			int n = Numeric.n;
			ColumnKernelWork work = factorColumnKernelWork(n);
			FlatColumnStore lStore = work.lStore;
			FlatColumnStore uStore = work.uStore;
			RowBuilder[] luRows = work.luRows;
			boolean materializeRows = Boolean.getBoolean(
					"jklu.complex.factor.materializeRows");
			boolean collectRows = materializeRows || profile.enabled();
			work.clearBuilders(collectRows);
			lStore.clear();
			uStore.clear();
			int expectedFactorEntries = expectedFactorEntries(n, Numeric.nz);
			lStore.ensureCapacity(expectedFactorEntries);
			uStore.ensureCapacity(expectedFactorEntries);
			double[] xr = work.xr;
			double[] xi = work.xi;
			int[] mark = work.mark;
			IntList touched = work.touched;
			ActiveList active = work.active;
			touched.ensureCapacity(n);
			active.ensureCapacity(n);
			active.resetStats();
			int[] pivots = work.pivots;
			boolean profilePhases = profile.enabled();
			Numeric.Udiag = new double[2 * n];

			for (int block = 0; block < Numeric.nblocks; block++)
			{
				int k1 = Numeric.R[block];
				int k2 = Numeric.R[block + 1];
					if (k2 - k1 == 1)
					{
						int k = k1;
						lStore.startColumn(k);
						uStore.startColumn(k);
					int inputStart = initialColumns.start(k);
					int inputEnd = initialColumns.end(k);
					double pr = 0.0;
					double pi = 0.0;
					if (isDirectDiagonalColumn(initialColumns, k, inputStart, inputEnd))
					{
						pr = initialColumns.real(k, inputStart);
						pi = initialColumns.imag(k, inputStart);
					}
					else
					{
						for (int p = inputStart; p < inputEnd; p++)
						{
							if (initialColumns.row(k, p) == k)
							{
								pr += initialColumns.real(k, p);
								pi += initialColumns.imag(k, p);
							}
						}
					}
						if (pr == 0.0 && pi == 0.0)
						{
							Common.status = KLU_SINGULAR;
						if (Common.numerical_rank == EMPTY)
						{
							Common.numerical_rank = k;
							Common.singular_col = k;
						}
						if (Common.halt_if_singular != 0)
						{
							return false;
						}
					}
					pivots[k] = k;
					Numeric.Udiag[2 * k] = pr;
					Numeric.Udiag[2 * k + 1] = pi;
					if (collectRows)
					{
						luRows[k].add(k, pr, pi);
					}
					if (profilePhases)
					{
						profile.step(k, inputEnd - inputStart, 0, 0, 0, luRows);
					}
					continue;
					}
					for (int k = k1; k < k2; k++)
					{
						int inputStart = initialColumns.start(k);
						int inputEnd = initialColumns.end(k);
						if (isDirectDiagonalColumn(initialColumns, k, inputStart, inputEnd) &&
								!isZero(initialColumns.real(k, inputStart),
										initialColumns.imag(k, inputStart)))
						{
							recordDirectDiagonalColumn(Numeric, lStore, uStore, luRows,
									collectRows, profile, k, initialColumns.real(k, inputStart),
									initialColumns.imag(k, inputStart));
							continue;
						}
						int generation = work.nextGeneration();
						long phase0 = profilePhases ? System.nanoTime() : 0L;
						for (int p = inputStart; p < inputEnd; p++)
					{
						int row = initialColumns.row(k, p);
							if (row >= k1 && row < k2)
							{
								addWork(row, initialColumns.real(k, p),
										initialColumns.imag(k, p), xr, xi, mark,
										generation, touched, active, k);
							}
					}
					long scatterNs = profilePhases ? System.nanoTime() - phase0 : 0L;

					uStore.startColumn(k);
					int updateCount = 0;
					phase0 = profilePhases ? System.nanoTime() : 0L;
					while (true)
					{
						int row = active.pollMinBelow(k);
						if (row == EMPTY)
						{
							break;
						}
						double ur = xr[row];
						double ui = xi[row];
						xr[row] = 0.0;
						xi[row] = 0.0;
						uStore.add(row, ur, ui);
						if (ur == 0.0 && ui == 0.0)
						{
							continue;
						}
							int lStart = lStore.p[row];
							int lSize = lStore.len[row];
							int[] targets = lStore.index;
							double[] lReal = lStore.real;
							double[] lImag = lStore.imag;
							if (profilePhases)
							{
								for (int p = 0; p < lSize; p++)
								{
									int q = lStart + p;
									int target = targets[q];
									double lr = lReal[q];
									double li = lImag[q];
								profile.updateShape(ur, ui, lr, li);
								profile.updateTarget(target < k);
									boolean prior = target < k;
							if (mark[target] != generation)
							{
								mark[target] = generation;
								touched.add(target);
								if (prior)
								{
									active.add(target);
									}
									profile.newTouch(prior);
								}
								xr[target] -= lr * ur - li * ui;
								xi[target] -= li * ur + lr * ui;
								updateCount++;
							}
							}
							else
							{
								for (int p = 0; p < lSize; p++)
								{
									int q = lStart + p;
									int target = targets[q];
									double lr = lReal[q];
									double li = lImag[q];
							if (mark[target] != generation)
							{
								mark[target] = generation;
								touched.add(target);
									if (target < k)
									{
										active.add(target);
										}
									}
									xr[target] -= lr * ur - li * ui;
									xi[target] -= li * ur + lr * ui;
								}
								updateCount += lSize;
											}
										}
								long updateNs = profilePhases ? System.nanoTime() - phase0 : 0L;

					phase0 = profilePhases ? System.nanoTime() : 0L;
					double pivotAbs2 = absSquared(xr[k], xi[k]);
					double maxAbs2 = pivotAbs2;
					for (int p = 0; p < touched.size; p++)
					{
							int row = touched.values[p];
							if (row >= k && row < k2)
							{
								double candidateAbs2 = xr[row] * xr[row] + xi[row] * xi[row];
								if (candidateAbs2 > maxAbs2)
								{
									maxAbs2 = candidateAbs2;
								}
							}
						}
						if (pivotAbs2 == 0.0 ||
								pivotAbs2 < Common.tol * Common.tol * maxAbs2)
						{
							clearWork(xr, xi, touched, active);
								return false;
						}
					long pivotNs = profilePhases ? System.nanoTime() - phase0 : 0L;

					pivots[k] = k;
					double pivotR = xr[k];
					double pivotI = xi[k];
					Numeric.Udiag[2 * k] = pivotR;
					Numeric.Udiag[2 * k + 1] = pivotI;
						if (collectRows)
						{
							luRows[k].add(k, pivotR, pivotI);
							int uStart = uStore.p[k];
							int uLen = uStore.len[k];
							for (int p = 0; p < uLen; p++)
							{
								int q = uStart + p;
								luRows[uStore.index[q]].add(k, uStore.real[q], uStore.imag[q]);
							}
						}

					lStore.startColumn(k);
					phase0 = profilePhases ? System.nanoTime() : 0L;
					double denom = pivotR * pivotR + pivotI * pivotI;
					double invDenom = 1.0 / denom;
						for (int p = 0; p < touched.size; p++)
						{
							int row = touched.values[p];
							if (row <= k || row >= k2)
						{
							continue;
						}
							double lr = (xr[row] * pivotR + xi[row] * pivotI) * invDenom;
							double li = (xi[row] * pivotR - xr[row] * pivotI) * invDenom;
							lStore.add(row, lr, li);
						if (collectRows)
						{
							luRows[row].add(k, lr, li);
							}
						}
						long scaleNs = profilePhases ? System.nanoTime() - phase0 : 0L;
					if (profilePhases)
					{
						profile.step(k, touched.size, lStore.len[k], uStore.len[k],
								updateCount, luRows);
					}
					phase0 = profilePhases ? System.nanoTime() : 0L;
					clearWork(xr, xi, touched, active);
					if (profilePhases)
					{
						profile.directPhases(scatterNs, updateNs, pivotNs, scaleNs,
								System.nanoTime() - phase0);
					}
						if (profilePhases && profile.shouldStop(k))
						{
								return false;
						}
				}
			}

			Numeric.pivots = pivots;
				Numeric.PivotsIdentity = true;
				Numeric.Pnum = new int[n];
				Numeric.Pinv = new int[n];
			for (int i = 0; i < n; i++)
			{
				int oldRow = Numeric.RowPerm == null ? i : Numeric.RowPerm[i];
				Numeric.Pnum[i] = oldRow;
				Numeric.Pinv[oldRow] = i;
				}
				Common.status = KLU_OK;
				long finalizeStartNs = profilePhases ? System.nanoTime() : 0L;
				lStore.finish(n);
				uStore.finish(n);
				cacheDirectColumnFactors(Numeric, lStore, uStore);
			if (materializeRows)
			{
				populateSparsePatternState(Numeric, luRows);
			}
			else
			{
				setDirectCounts(Numeric);
				Numeric.LUrowsCurrent = false;
			}
			if (profilePhases)
			{
				profile.directFinalize(System.nanoTime() - finalizeStartNs);
				profile.directActive(active.addCount, active.pollCount);
			}
			if (collectRows)
			{
				profile.finish(luRows);
			}
			return true;
		}

		private static boolean isDirectDiagonalColumn(DirectColumns columns, int k,
				int inputStart, int inputEnd)
		{
			return inputEnd == inputStart + 1 && columns.row(k, inputStart) == k;
		}

		private static void recordDirectDiagonalColumn(KLU_z_numeric Numeric,
				FlatColumnStore lStore, FlatColumnStore uStore, RowBuilder[] luRows,
				boolean collectRows, FactorProfile profile, int k, double pr, double pi)
		{
			lStore.startColumn(k);
			uStore.startColumn(k);
			Numeric.Udiag[2 * k] = pr;
			Numeric.Udiag[2 * k + 1] = pi;
			if (collectRows)
			{
				luRows[k].add(k, pr, pi);
			}
			if (profile.enabled())
			{
				profile.step(k, 1, 0, 0, 0, luRows);
			}
		}

			private static boolean tryDirectRefactor(KLU_z_numeric Numeric,
				DirectColumns initialColumns, FactorProfile profile, KLU_common Common)
		{
			int n = Numeric.n;
			if (Numeric.DirectLp != null && Numeric.DirectLi != null &&
					Numeric.DirectLx != null && Numeric.DirectUp != null &&
					Numeric.DirectUi != null && Numeric.DirectUx != null)
			{
				return tryDirectRefactorFlat(Numeric, initialColumns, profile, Common);
			}
			if (Numeric.DirectLRows == null || Numeric.DirectURows == null)
			{
				return false;
			}
			ColumnKernelWork work = columnKernelWork(Numeric, n);
			RowBuilder[] luRows = work.luRows;
			boolean materializeRows = Boolean.getBoolean(
					"jklu.complex.refactor.materializeRows");
			if (materializeRows)
			{
				for (int row = 0; row < n; row++)
				{
					luRows[row].clear();
				}
			}
			double[] xr = work.xr;
			double[] xi = work.xi;
				boolean profilePhases = profile.enabled();
				Numeric.Udiag = Numeric.Udiag != null && Numeric.Udiag.length == 2 * n
						? Numeric.Udiag : new double[2 * n];
					for (int block = 0; block < Numeric.nblocks; block++)
			{
					int k1 = Numeric.R[block];
					int k2 = Numeric.R[block + 1];
					if (k2 - k1 == 1)
					{
						int k = k1;
						int inputStart = initialColumns.start(k);
						int inputEnd = initialColumns.end(k);
						double pr = 0.0;
						double pi = 0.0;
						for (int p = inputStart; p < inputEnd; p++)
						{
							if (initialColumns.row(k, p) == k)
							{
								pr += initialColumns.real(k, p);
								pi += initialColumns.imag(k, p);
							}
						}
						if (pr == 0.0 && pi == 0.0)
						{
							return false;
						}
						Numeric.Udiag[2 * k] = pr;
						Numeric.Udiag[2 * k + 1] = pi;
						if (materializeRows)
						{
							luRows[k].add(k, pr, pi);
						}
						if (profilePhases)
						{
							profile.step(k, inputEnd - inputStart, 0, 0, 0, luRows);
						}
						continue;
					}
					for (int k = k1; k < k2; k++)
					{
						long phase0 = profilePhases ? System.nanoTime() : 0L;
							int inputStart = initialColumns.start(k);
							int inputEnd = initialColumns.end(k);
							if (initialColumns.p != null && !initialColumns.cscView)
							{
								int[] inputRows = initialColumns.i;
								double[] inputValues = initialColumns.x;
								for (int p = inputStart; p < inputEnd; p++)
								{
									int row = inputRows[p];
									int xp = 2 * p;
									xr[row] += inputValues[xp];
									xi[row] += inputValues[xp + 1];
								}
							}
							else
							{
								for (int p = inputStart; p < inputEnd; p++)
								{
									int row = initialColumns.row(k, p);
									xr[row] += initialColumns.real(k, p);
									xi[row] += initialColumns.imag(k, p);
								}
							}
						long scatterNs = profilePhases ? System.nanoTime() - phase0 : 0L;

							int[] uRows = Numeric.DirectURows[k];
							double[] uReal = Numeric.DirectUReal[k];
							double[] uImag = Numeric.DirectUImag[k];
							int uLen = uRows.length;
						int updateCount = 0;
						phase0 = profilePhases ? System.nanoTime() : 0L;
							for (int up = 0; up < uLen; up++)
							{
								int row = uRows[up];
							double ur = xr[row];
							double ui = xi[row];
							uReal[up] = ur;
							uImag[up] = ui;
							xr[row] = 0.0;
							xi[row] = 0.0;
								if (materializeRows)
								{
									luRows[row].add(k, ur, ui);
								}
								if (ur == 0.0 && ui == 0.0)
								{
									continue;
								}

							int[] lRows = Numeric.DirectLRows[row];
							double[] lReal = Numeric.DirectLReal[row];
							double[] lImag = Numeric.DirectLImag[row];
							int lUpdateLen = lRows.length;
						for (int lp = 0; lp < lUpdateLen; lp++)
						{
							int target = lRows[lp];
								double lr = lReal[lp];
								double li = lImag[lp];
								xr[target] -= lr * ur - li * ui;
								xi[target] -= li * ur + lr * ui;
								updateCount++;
							}
						}
						long updateNs = profilePhases ? System.nanoTime() - phase0 : 0L;

						phase0 = profilePhases ? System.nanoTime() : 0L;
						if (absSquared(xr[k], xi[k]) == 0.0)
						{
									initialColumns.clear(xr, xi, k);
									clearRefactorColumn(xr, xi, uRows, k,
											Numeric.DirectLRows[k], uLen,
											Numeric.DirectLRows[k].length);
							return false;
					}
					Numeric.Udiag[2 * k] = xr[k];
					Numeric.Udiag[2 * k + 1] = xi[k];
						if (materializeRows)
						{
							luRows[k].add(k, xr[k], xi[k]);
						}
						long pivotNs = profilePhases ? System.nanoTime() - phase0 : 0L;

							int[] lRows = Numeric.DirectLRows[k];
							double[] lReal = Numeric.DirectLReal[k];
							double[] lImag = Numeric.DirectLImag[k];
							int lLen = lRows.length;
						double denom = xr[k] * xr[k] + xi[k] * xi[k];
						double invDenom = 1.0 / denom;
						phase0 = profilePhases ? System.nanoTime() : 0L;
						for (int lp = 0; lp < lLen; lp++)
						{
							int row = lRows[lp];
						double lr = (xr[row] * xr[k] + xi[row] * xi[k]) * invDenom;
						double li = (xi[row] * xr[k] - xr[row] * xi[k]) * invDenom;
						lReal[lp] = lr;
						lImag[lp] = li;
						if (materializeRows)
						{
								luRows[row].add(k, lr, li);
							}
						}
						long scaleNs = profilePhases ? System.nanoTime() - phase0 : 0L;
						if (profilePhases)
						{
							profile.step(k, inputEnd - inputStart, lLen, uLen,
									updateCount, luRows);
						}
						phase0 = profilePhases ? System.nanoTime() : 0L;
						initialColumns.clear(xr, xi, k);
						clearRefactorColumn(xr, xi, uRows, k, lRows, uLen, lLen);
						if (profilePhases)
						{
							profile.directPhases(scatterNs, updateNs, pivotNs, scaleNs,
									System.nanoTime() - phase0);
						}
						if (profilePhases && profile.shouldStop(k))
						{
							return false;
						}
				}
			}
			Common.status = KLU_OK;
			if (materializeRows)
			{
				populateSparsePatternState(Numeric, luRows);
				profile.finish(luRows);
			}
			else
			{
				Numeric.LUrowsCurrent = false;
			}
			return true;
		}

		private static boolean tryDirectRefactorFlat(KLU_z_numeric Numeric,
				DirectColumns initialColumns, FactorProfile profile, KLU_common Common)
		{
			int n = Numeric.n;
			ColumnKernelWork work = factorColumnKernelWork(n);
			RowBuilder[] luRows = work.luRows;
			boolean materializeRows = Boolean.getBoolean(
					"jklu.complex.refactor.materializeRows");
			if (materializeRows)
			{
				for (int row = 0; row < n; row++)
				{
					luRows[row].clear();
				}
			}
			double[] xr = work.xr;
			double[] xi = work.xi;
			boolean profilePhases = profile.enabled();
			Numeric.Udiag = Numeric.Udiag != null && Numeric.Udiag.length == 2 * n
					? Numeric.Udiag : new double[2 * n];

			for (int block = 0; block < Numeric.nblocks; block++)
			{
				int k1 = Numeric.R[block];
					int k2 = Numeric.R[block + 1];
					if (k2 - k1 == 1)
					{
						int k = k1;
						int inputStart = initialColumns.start(k);
						int inputEnd = initialColumns.end(k);
						double pr = 0.0;
						double pi = 0.0;
						if (isDirectDiagonalColumn(initialColumns, k, inputStart, inputEnd))
						{
							pr = initialColumns.real(k, inputStart);
							pi = initialColumns.imag(k, inputStart);
						}
						else
						{
							for (int p = inputStart; p < inputEnd; p++)
							{
								if (initialColumns.row(k, p) == k)
								{
									pr += initialColumns.real(k, p);
									pi += initialColumns.imag(k, p);
								}
							}
						}
							if (pr == 0.0 && pi == 0.0)
							{
								return false;
							}
						Numeric.Udiag[2 * k] = pr;
						Numeric.Udiag[2 * k + 1] = pi;
						if (materializeRows)
						{
							luRows[k].add(k, pr, pi);
						}
						if (profilePhases)
						{
							profile.step(k, inputEnd - inputStart, 0, 0, 0, luRows);
						}
						continue;
					}
					for (int k = k1; k < k2; k++)
						{
						long phase0 = profilePhases ? System.nanoTime() : 0L;
							int inputStart = initialColumns.start(k);
							int inputEnd = initialColumns.end(k);
							if (initialColumns.p != null)
							{
								int[] inputRows = initialColumns.i;
								double[] inputValues = initialColumns.x;
								for (int p = inputStart; p < inputEnd; p++)
								{
									int row = inputRows[p];
									int xp = 2 * p;
									xr[row] += inputValues[xp];
									xi[row] += inputValues[xp + 1];
								}
							}
							else
							{
								for (int p = inputStart; p < inputEnd; p++)
								{
									int row = initialColumns.row(k, p);
									xr[row] += initialColumns.real(k, p);
									xi[row] += initialColumns.imag(k, p);
								}
							}
						long scatterNs = profilePhases ? System.nanoTime() - phase0 : 0L;

						int updateCount = 0;
						phase0 = profilePhases ? System.nanoTime() : 0L;
							for (int up = Numeric.DirectUp[k]; up < Numeric.DirectUp[k + 1]; up++)
							{
								int row = Numeric.DirectUi[up];
							double ur = xr[row];
							double ui = xi[row];
							Numeric.DirectUx[2 * up] = ur;
							Numeric.DirectUx[2 * up + 1] = ui;
							xr[row] = 0.0;
							xi[row] = 0.0;
							if (materializeRows)
							{
								luRows[row].add(k, ur, ui);
							}
								if (ur == 0.0 && ui == 0.0)
								{
									continue;
								}

									for (int lp = Numeric.DirectLp[row]; lp < Numeric.DirectLp[row + 1]; lp++)
									{
										int target = Numeric.DirectLi[lp];
										double lr = Numeric.DirectLx[2 * lp];
										double li = Numeric.DirectLx[2 * lp + 1];
										xr[target] -= lr * ur - li * ui;
										xi[target] -= li * ur + lr * ui;
										updateCount++;
									}
							}
						long updateNs = profilePhases ? System.nanoTime() - phase0 : 0L;

						phase0 = profilePhases ? System.nanoTime() : 0L;
							if (absSquared(xr[k], xi[k]) == 0.0)
							{
								initialColumns.clear(xr, xi, k);
								clearRefactorColumn(xr, xi, Numeric.DirectUi,
										Numeric.DirectUp[k], Numeric.DirectUp[k + 1], k,
								Numeric.DirectLi, Numeric.DirectLp[k], Numeric.DirectLp[k + 1]);
						return false;
					}
						Numeric.Udiag[2 * k] = xr[k];
						Numeric.Udiag[2 * k + 1] = xi[k];
						double pivotR = xr[k];
						double pivotI = xi[k];
							if (materializeRows)
							{
								luRows[k].add(k, pivotR, pivotI);
							}
							long pivotNs = profilePhases ? System.nanoTime() - phase0 : 0L;

							double denom = pivotR * pivotR + pivotI * pivotI;
							double invDenom = 1.0 / denom;
							phase0 = profilePhases ? System.nanoTime() : 0L;
							for (int lp = Numeric.DirectLp[k]; lp < Numeric.DirectLp[k + 1]; lp++)
							{
								int row = Numeric.DirectLi[lp];
							double lr = (xr[row] * pivotR + xi[row] * pivotI) * invDenom;
							double li = (xi[row] * pivotR - xr[row] * pivotI) * invDenom;
								Numeric.DirectLx[2 * lp] = lr;
								Numeric.DirectLx[2 * lp + 1] = li;
							xr[row] = 0.0;
							xi[row] = 0.0;
								if (materializeRows)
							{
									luRows[row].add(k, lr, li);
							}
						}
						long scaleNs = profilePhases ? System.nanoTime() - phase0 : 0L;
						if (profilePhases)
						{
							profile.step(k, inputEnd - inputStart,
									Numeric.DirectLp[k + 1] - Numeric.DirectLp[k],
									Numeric.DirectUp[k + 1] - Numeric.DirectUp[k],
									updateCount, luRows);
						}
						phase0 = profilePhases ? System.nanoTime() : 0L;
						if (Numeric.nzoff != 0)
						{
							initialColumns.clear(xr, xi, k);
						}
						xr[k] = 0.0;
						xi[k] = 0.0;
						if (profilePhases)
						{
							profile.directPhases(scatterNs, updateNs, pivotNs, scaleNs,
									System.nanoTime() - phase0);
						}
						if (profilePhases && profile.shouldStop(k))
						{
						return false;
					}
				}
			}
			Common.status = KLU_OK;
			if (materializeRows)
			{
				populateSparsePatternState(Numeric, luRows);
				profile.finish(luRows);
			}
			else
			{
				Numeric.LUrowsCurrent = false;
			}
			return true;
		}

		private static void cacheDirectColumnFactors(KLU_z_numeric Numeric,
				FlatColumnStore lStore, FlatColumnStore uStore)
		{
			int n = Numeric.n;
			boolean cacheFlat = cacheFlatColumns(n);
			boolean keepColumnArrays = !cacheFlat ||
					Boolean.getBoolean("jklu.complex.factor.keepColumnArrays");
			if (keepColumnArrays)
			{
				Numeric.DirectLRows = new int[n][];
				Numeric.DirectLReal = new double[n][];
				Numeric.DirectLImag = new double[n][];
				Numeric.DirectURows = new int[n][];
				Numeric.DirectUReal = new double[n][];
				Numeric.DirectUImag = new double[n][];
			}
			else
			{
				Numeric.DirectLRows = null;
				Numeric.DirectLReal = null;
				Numeric.DirectLImag = null;
				Numeric.DirectURows = null;
				Numeric.DirectUReal = null;
				Numeric.DirectUImag = null;
			}
			if (cacheFlat)
			{
				lStore.p[n] = lStore.size;
				uStore.p[n] = uStore.size;
				Numeric.DirectLp = Arrays.copyOf(lStore.p, n + 1);
				Numeric.DirectUp = Arrays.copyOf(uStore.p, n + 1);
				Numeric.DirectLi = Arrays.copyOf(lStore.index, lStore.size);
				Numeric.DirectUi = Arrays.copyOf(uStore.index, uStore.size);
				Numeric.DirectLx = interleavedCopy(lStore.real, lStore.imag,
						lStore.size);
				Numeric.DirectUx = interleavedCopy(uStore.real, uStore.imag,
						uStore.size);
			}
			else
			{
				Numeric.DirectLp = null;
				Numeric.DirectLi = null;
				Numeric.DirectLx = null;
				Numeric.DirectUp = null;
				Numeric.DirectUi = null;
				Numeric.DirectUx = null;
			}
			for (int k = 0; k < n; k++)
			{
				if (keepColumnArrays)
				{
					Numeric.DirectLRows[k] = copy(lStore.index, lStore.p[k], lStore.len[k]);
					Numeric.DirectLReal[k] = copy(lStore.real, lStore.p[k], lStore.len[k]);
					Numeric.DirectLImag[k] = copy(lStore.imag, lStore.p[k], lStore.len[k]);
					Numeric.DirectURows[k] = copy(uStore.index, uStore.p[k], uStore.len[k]);
					Numeric.DirectUReal[k] = copy(uStore.real, uStore.p[k], uStore.len[k]);
					Numeric.DirectUImag[k] = copy(uStore.imag, uStore.p[k], uStore.len[k]);
				}
			}
		}

		private static double[] interleavedCopy(double[] real, double[] imag, int size)
		{
			double[] values = new double[2 * size];
			for (int p = 0; p < size; p++)
			{
				int q = 2 * p;
				values[q] = real[p];
				values[q + 1] = imag[p];
			}
			return values;
		}

		private static boolean cacheFlatColumns(int n)
		{
			if (Boolean.getBoolean("jklu.complex.factor.cacheFlatColumns"))
			{
				return true;
			}
			int maxN = Integer.getInteger("jklu.complex.factor.flatColumnMaxN", 100000)
					.intValue();
			return n <= maxN;
		}

			private static boolean useReachKernel(int n)
			{
				Boolean configured = configuredBoolean("jklu.complex.factor.reachKernel");
				if (configured != null)
				{
				return configured.booleanValue();
			}
			int minN = Integer.getInteger("jklu.complex.factor.reachKernelMinN",
						0).intValue();
				return n >= minN;
			}

			private static int expectedFactorEntries(int n, int nz)
			{
				int multiplier = Integer.getInteger(
						"jklu.complex.factor.storeInitialNnzMultiplier",
						2).intValue();
				if (multiplier < 1)
				{
					multiplier = 1;
				}
				long expected = Math.max((long) n, (long) nz * multiplier);
				return (int) Math.min(Integer.MAX_VALUE,
						Math.max(1024L, expected));
			}

			private static boolean useReachUnsortedU(int n)
		{
			Boolean configured = configuredBoolean("jklu.complex.factor.reachUnsortedU");
			if (configured != null)
			{
				return configured.booleanValue();
			}
					return true;
				}

		private static boolean useReachSortCache(int n)
		{
			Boolean configured = configuredBoolean("jklu.complex.factor.reachSortCache");
			if (configured != null)
			{
				return configured.booleanValue();
			}
			return false;
		}

		private static Boolean configuredBoolean(String name)
		{
			String value = System.getProperty(name);
			if (value == null)
			{
				return null;
			}
			return Boolean.valueOf(value);
		}

		private static void copyColumnToFlat(ColumnBuilder column, int offset,
				int[] indices, double[] values)
		{
			System.arraycopy(column.index, 0, indices, offset, column.size);
			for (int p = 0; p < column.size; p++)
			{
				int q = offset + p;
				values[2 * q] = column.real[p];
				values[2 * q + 1] = column.imag[p];
			}
		}

		private static void copyColumnToFlat(FlatColumnStore store, int col, int offset,
				int[] indices, double[] values)
		{
			int start = store.p[col];
			int size = store.len[col];
			System.arraycopy(store.index, start, indices, offset, size);
			for (int p = 0; p < size; p++)
			{
				int source = start + p;
				int target = offset + p;
				values[2 * target] = store.real[source];
				values[2 * target + 1] = store.imag[source];
			}
		}

		private static boolean pivotsIdentity(int[] pivots)
		{
			for (int i = 0; i < pivots.length; i++)
			{
				if (pivots[i] != i)
				{
					return false;
				}
			}
			return true;
		}

		private static double absSquared(double re, double im)
		{
			return re * re + im * im;
		}

		static void refreshUdiagInverse(KLU_z_numeric Numeric)
		{
			if (!Boolean.getBoolean("jklu.complex.factor.precomputeUdiagInv"))
			{
				Numeric.UdiagInv = null;
				return;
			}
			int n = Numeric.n;
			if (Numeric.Udiag == null || Numeric.Udiag.length < 2 * n)
			{
				Numeric.UdiagInv = null;
				return;
			}
			if (Numeric.UdiagInv == null || Numeric.UdiagInv.length < 2 * n)
			{
				Numeric.UdiagInv = new double[2 * n];
			}
			for (int k = 0; k < n; k++)
			{
				int p = 2 * k;
				double br = Numeric.Udiag[p];
				double bi = Numeric.Udiag[p + 1];
				double denom = br * br + bi * bi;
				if (denom == 0.0)
				{
					Numeric.UdiagInv[p] = 0.0;
					Numeric.UdiagInv[p + 1] = 0.0;
				}
				else
				{
					Numeric.UdiagInv[p] = br / denom;
					Numeric.UdiagInv[p + 1] = -bi / denom;
				}
			}
		}

		private static void setDirectCounts(KLU_z_numeric Numeric)
		{
			int lnz = Numeric.n;
			int unz = Numeric.n;
			if (Numeric.DirectLp != null)
			{
				lnz += Numeric.DirectLp[Numeric.n];
			}
				else if (Numeric.DirectLRows != null)
				{
					for (int col = 0; col < Numeric.DirectLRows.length; col++)
					{
						lnz += Numeric.DirectLRows[col].length;
					}
				}
			if (Numeric.DirectUp != null)
			{
				unz += Numeric.DirectUp[Numeric.n];
			}
				else if (Numeric.DirectURows != null)
				{
					for (int col = 0; col < Numeric.DirectURows.length; col++)
					{
						unz += Numeric.DirectURows[col].length;
					}
				}
			Numeric.lnz = lnz;
			Numeric.unz = unz;
			Numeric.Work = new double[8 * Math.max(1, Numeric.n)];
			Numeric.Xwork = Numeric.Work;
		}

		private static boolean hasDirectFactors(KLU_z_numeric Numeric)
		{
			return (Numeric.DirectLp != null && Numeric.DirectLi != null &&
					Numeric.DirectLx != null && Numeric.DirectUp != null &&
					Numeric.DirectUi != null && Numeric.DirectUx != null) ||
					(Numeric.DirectLRows != null && Numeric.DirectURows != null);
		}

		private static ColumnKernelWork factorColumnKernelWork(int n)
		{
			ColumnKernelWork work = FACTOR_WORK.get();
			if (work == null || work.n != n)
			{
				work = new ColumnKernelWork(n);
				FACTOR_WORK.set(work);
			}
			return work;
		}

			private static ColumnKernelWork columnKernelWork(KLU_z_numeric Numeric, int n)
		{
			if (Numeric.DirectKernelWork instanceof ColumnKernelWork)
			{
				ColumnKernelWork work = (ColumnKernelWork) Numeric.DirectKernelWork;
				if (work.n == n)
				{
					return work;
				}
			}
			ColumnKernelWork work = new ColumnKernelWork(n);
			Numeric.DirectKernelWork = work;
			return work;
		}

		private static DirectColumns buildDirectColumnKernelInput(KLU_z_numeric Numeric,
				KLU_common Common)
		{
			int n = Numeric.n;
			if (Boolean.getBoolean("jklu.complex.factor.disableColumnKernel"))
			{
				return null;
			}
			if (!allowDirectColumnKernelForBtf(Numeric))
			{
				return null;
			}
			int minN = Integer.getInteger("jklu.complex.factor.columnKernelMinN",
					10000).intValue();
				if (n < minN && Numeric.nblocks <= 1)
				{
					return null;
				}
				if (Numeric.Rs != null)
			{
				if (klu_z_scale(Common.scale, n, Numeric.Ap, Numeric.Ai, Numeric.Ax,
						Numeric.Rs, new int[n], Common) == FALSE)
				{
						return null;
					}
						}
						if (!hasDirectFactors(Numeric) && streamFirstFactorInput(n))
						{
							return new DirectColumns(Numeric);
						}
					boolean builtPattern = false;
					if (Numeric.DirectColumnRows == null && Numeric.DirectColumnP == null)
					{
						buildDirectColumnPattern(Numeric, Common);
						if (Common.status < KLU_OK)
					{
						return null;
					}
						builtPattern = true;
				}
				if (!builtPattern)
				{
					refillDirectColumnValues(Numeric, Common);
					if (Common.status < KLU_OK)
					{
						return null;
					}
				}
			if (Numeric.DirectColumnP != null)
			{
				return new DirectColumns(Numeric.DirectColumnP, Numeric.DirectColumnI,
						Numeric.DirectColumnX);
			}
			return new DirectColumns(Numeric.DirectColumnRows, Numeric.DirectColumnReal,
					Numeric.DirectColumnImag);
		}

		private static boolean allowDirectColumnKernelForBtf(KLU_z_numeric Numeric)
		{
			return !Boolean.getBoolean(
					"jklu.complex.factor.disableFragmentedBtfColumnKernel");
		}

		private static boolean hasBlockDiagonalCandidates(KLU_z_numeric Numeric,
				DirectColumns columns)
		{
			if (Numeric.DirectHasBlockDiagonalCandidates != 0)
			{
				return Numeric.DirectHasBlockDiagonalCandidates > 0;
			}
			for (int block = 0; block < Numeric.nblocks; block++)
			{
				int k1 = Numeric.R[block];
				int k2 = Numeric.R[block + 1];
				for (int col = k1; col < k2; col++)
				{
					boolean found = false;
					for (int p = columns.start(col); p < columns.end(col); p++)
					{
						if (columns.row(col, p) == col)
						{
							found = true;
							break;
						}
					}
					if (!found)
					{
						Numeric.DirectHasBlockDiagonalCandidates = -1;
						return false;
					}
				}
			}
			Numeric.DirectHasBlockDiagonalCandidates = 1;
			return true;
		}

		private static void buildDirectColumnPattern(KLU_z_numeric Numeric,
				KLU_common Common)
		{
			int n = Numeric.n;
			boolean packedInput = packedDirectInput(n);
			if (packedInput && Numeric.nzoff == 0)
			{
				buildPackedNoOffDiagonalPattern(Numeric);
				Common.status = KLU_OK;
				return;
			}
					int[] counts = new int[n];
					int[] offCounts = new int[n];
					boolean[] diagonalFound = new boolean[n];
					for (int oldCol = 0; oldCol < n; oldCol++)
					{
						int col = Numeric.ColInv == null ? oldCol : Numeric.ColInv[oldCol];
					for (int p = Numeric.Ap[oldCol]; p < Numeric.Ap[oldCol + 1]; p++)
				{
					int oldRow = Numeric.Ai[p];
					if (oldRow < 0 || oldRow >= n)
					{
						Common.status = KLU_INVALID;
						return;
					}
					int row = Numeric.RowInv == null ? oldRow : Numeric.RowInv[oldRow];
						if (Numeric.BlockOf[row] == Numeric.BlockOf[col])
						{
							counts[col]++;
							if (row == col)
							{
								diagonalFound[col] = true;
							}
						}
						else
						{
						offCounts[col]++;
				}
					}
				}
				Numeric.DirectHasBlockDiagonalCandidates = allTrue(diagonalFound) ? 1 : -1;
				Numeric.DirectColumnRows = new int[n][];
			Numeric.DirectColumnReal = new double[n][];
			Numeric.DirectColumnImag = new double[n][];
			if (packedInput)
			{
				Numeric.DirectColumnRows = null;
				Numeric.DirectColumnReal = null;
				Numeric.DirectColumnImag = null;
				Numeric.DirectColumnP = new int[n + 1];
				for (int col = 0; col < n; col++)
				{
					Numeric.DirectColumnP[col + 1] = Numeric.DirectColumnP[col] + counts[col];
				}
				Numeric.DirectColumnI = new int[Numeric.DirectColumnP[n]];
				Numeric.DirectColumnX = new double[2 * Numeric.DirectColumnP[n]];
			}
			else
			{
			for (int col = 0; col < n; col++)
			{
				Numeric.DirectColumnRows[col] = new int[counts[col]];
				Numeric.DirectColumnReal[col] = new double[counts[col]];
				Numeric.DirectColumnImag[col] = new double[counts[col]];
			}
			}
			Numeric.Offp = new int[n + 1];
				for (int col = 0; col < n; col++)
				{
					Numeric.Offp[col + 1] = Numeric.Offp[col] + offCounts[col];
				}
				Numeric.nzoff = Numeric.Offp[n];
			Numeric.Offi = new int[Numeric.nzoff];
			Numeric.Offx = new double[2 * Numeric.nzoff];
			Numeric.DirectColumnNext = packedInput ? copy(Numeric.DirectColumnP, n) : new int[n];
			int[] offNext = Numeric.nzoff == 0 ? null : copy(Numeric.Offp, n);
		for (int oldCol = 0; oldCol < n; oldCol++)
		{
			int col = Numeric.ColInv == null ? oldCol : Numeric.ColInv[oldCol];
			for (int p = Numeric.Ap[oldCol]; p < Numeric.Ap[oldCol + 1]; p++)
				{
					int oldRow = Numeric.Ai[p];
					int row = Numeric.RowInv == null ? oldRow : Numeric.RowInv[oldRow];
					double scale = Numeric.Rs == null ? 1.0 : Numeric.Rs[oldRow];
					double re = Numeric.Ax[2 * p] / scale;
					double im = Numeric.Ax[2 * p + 1] / scale;
					if (Numeric.BlockOf[row] == Numeric.BlockOf[col])
					{
						if (packedInput)
						{
							int q = Numeric.DirectColumnNext[col]++;
							Numeric.DirectColumnI[q] = row;
							Numeric.DirectColumnX[2 * q] = re;
							Numeric.DirectColumnX[2 * q + 1] = im;
						}
						else
						{
							int q = Numeric.DirectColumnNext[col]++;
							Numeric.DirectColumnRows[col][q] = row;
							Numeric.DirectColumnReal[col][q] = re;
							Numeric.DirectColumnImag[col][q] = im;
						}
					}
					else
					{
						if (offNext != null)
						{
							int q = offNext[col]++;
							Numeric.Offi[q] = row;
							Numeric.Offx[2 * q] = re;
							Numeric.Offx[2 * q + 1] = im;
						}
					}
				}
		}
	}

		private static void buildPackedNoOffDiagonalPattern(KLU_z_numeric Numeric)
		{
			int n = Numeric.n;
			Numeric.DirectColumnRows = null;
			Numeric.DirectColumnReal = null;
			Numeric.DirectColumnImag = null;
					Numeric.DirectColumnP = new int[n + 1];
					int[] counts = new int[n];
					boolean[] diagonalFound = new boolean[n];
					for (int oldCol = 0; oldCol < n; oldCol++)
					{
						int col = Numeric.ColInv == null ? oldCol : Numeric.ColInv[oldCol];
						counts[col] = Numeric.Ap[oldCol + 1] - Numeric.Ap[oldCol];
					}
			for (int col = 0; col < n; col++)
			{
				Numeric.DirectColumnP[col + 1] = Numeric.DirectColumnP[col] + counts[col];
			}
			Numeric.DirectColumnI = new int[Numeric.DirectColumnP[n]];
			Numeric.DirectColumnX = new double[2 * Numeric.DirectColumnP[n]];
			Numeric.Offp = new int[n + 1];
			Numeric.Offi = new int[0];
			Numeric.Offx = new double[0];
			Numeric.DirectColumnNext = copy(Numeric.DirectColumnP, n);
			for (int oldCol = 0; oldCol < n; oldCol++)
			{
				int col = Numeric.ColInv == null ? oldCol : Numeric.ColInv[oldCol];
					for (int p = Numeric.Ap[oldCol]; p < Numeric.Ap[oldCol + 1]; p++)
					{
							int oldRow = Numeric.Ai[p];
							int row = Numeric.RowInv == null ? oldRow : Numeric.RowInv[oldRow];
								int q = Numeric.DirectColumnNext[col]++;
								double scale = Numeric.Rs == null ? 1.0 : Numeric.Rs[oldRow];
								Numeric.DirectColumnI[q] = row;
								Numeric.DirectColumnX[2 * q] = Numeric.Ax[2 * p] / scale;
								Numeric.DirectColumnX[2 * q + 1] =
										Numeric.Ax[2 * p + 1] / scale;
								if (row == col)
								{
									diagonalFound[col] = true;
								}
							}
						}
				Numeric.DirectHasBlockDiagonalCandidates = allTrue(diagonalFound) ? 1 : -1;
			}

		private static void refillDirectColumnValues(KLU_z_numeric Numeric,
				KLU_common Common)
		{
			int n = Numeric.n;
			int[] next = Numeric.DirectColumnNext;
			if (Numeric.DirectColumnP != null && Numeric.nzoff == 0)
			{
				refillPackedNoOffDiagonalValues(Numeric);
				Common.status = KLU_OK;
				return;
			}
			if (Numeric.DirectColumnP != null)
			{
				System.arraycopy(Numeric.DirectColumnP, 0, next, 0, n);
			}
			else
			{
				for (int col = 0; col < n; col++)
				{
					next[col] = 0;
					}
			}
			int[] offNext = Numeric.nzoff == 0 ? null : copy(Numeric.Offp, n);
			for (int oldCol = 0; oldCol < n; oldCol++)
			{
				int col = Numeric.ColInv == null ? oldCol : Numeric.ColInv[oldCol];
				for (int p = Numeric.Ap[oldCol]; p < Numeric.Ap[oldCol + 1]; p++)
				{
					int oldRow = Numeric.Ai[p];
					if (oldRow < 0 || oldRow >= n)
					{
						Common.status = KLU_INVALID;
						return;
					}
					int row = Numeric.RowInv == null ? oldRow : Numeric.RowInv[oldRow];
					double scale = Numeric.Rs == null ? 1.0 : Numeric.Rs[oldRow];
					double re = Numeric.Ax[2 * p] / scale;
					double im = Numeric.Ax[2 * p + 1] / scale;
					if (Numeric.BlockOf[row] == Numeric.BlockOf[col])
					{
						int q = next[col]++;
						if (Numeric.DirectColumnP != null)
						{
							Numeric.DirectColumnX[2 * q] = re;
							Numeric.DirectColumnX[2 * q + 1] = im;
						}
						else
						{
							Numeric.DirectColumnReal[col][q] = re;
							Numeric.DirectColumnImag[col][q] = im;
						}
					}
					else
					{
						if (offNext != null)
						{
							int q = offNext[col]++;
							Numeric.Offx[2 * q] = re;
							Numeric.Offx[2 * q + 1] = im;
						}
					}
				}
			}
			Common.status = KLU_OK;
		}

		private static void refillPackedNoOffDiagonalValues(KLU_z_numeric Numeric)
		{
			int n = Numeric.n;
			int[] next = Numeric.DirectColumnNext;
			System.arraycopy(Numeric.DirectColumnP, 0, next, 0, n);
			for (int oldCol = 0; oldCol < n; oldCol++)
			{
				int col = Numeric.ColInv == null ? oldCol : Numeric.ColInv[oldCol];
				for (int p = Numeric.Ap[oldCol]; p < Numeric.Ap[oldCol + 1]; p++)
				{
					int q = next[col]++;
					double scale = Numeric.Rs == null ? 1.0 : Numeric.Rs[Numeric.Ai[p]];
					Numeric.DirectColumnX[2 * q] = Numeric.Ax[2 * p] / scale;
					Numeric.DirectColumnX[2 * q + 1] = Numeric.Ax[2 * p + 1] / scale;
				}
			}
		}

			private static boolean packedDirectInput(int n)
			{
				Boolean configured = configuredBoolean("jklu.complex.factor.packedInput");
			if (configured != null)
			{
				return configured.booleanValue();
			}
				return useReachKernel(n);
			}

			private static boolean streamFirstFactorInput(int n)
			{
				Boolean configured = configuredBoolean(
						"jklu.complex.factor.streamFirstInput");
				if (configured != null)
				{
					return configured.booleanValue();
				}
				return false;
			}

			private static Column[] columnsFromRows(SparseRow[] rows, int n)
		{
			ColumnBuilder[] builders = new ColumnBuilder[n];
			for (int col = 0; col < n; col++)
			{
				builders[col] = new ColumnBuilder();
			}
			for (int row = 0; row < n; row++)
			{
				for (Map.Entry<Integer, Complex> entry : rows[row].entrySet())
				{
					int col = entry.getKey().intValue();
					Complex value = entry.getValue();
					builders[col].add(row, value.re, value.im);
				}
			}
			Column[] columns = new Column[n];
			for (int col = 0; col < n; col++)
			{
				columns[col] = builders[col].toColumn();
			}
			return columns;
		}

		private static void addWork(int row, double re, double im, double[] xr,
				double[] xi, int[] mark, int generation, IntList touched,
				ActiveList active, int orderedLimit)
		{
			if (mark[row] != generation)
			{
				mark[row] = generation;
				touched.add(row);
				if (row < orderedLimit)
				{
					active.add(row);
				}
			}
			xr[row] += re;
			xi[row] += im;
		}

		private static void addReachWork(int row, double re, double im, double[] xr,
				double[] xi, int[] mark, int generation, IntList touched)
		{
			if (mark[row] != generation)
			{
				mark[row] = generation;
				touched.add(row);
			}
			xr[row] += re;
			xi[row] += im;
		}

		private static void clearRefactorColumn(double[] xr, double[] xi,
				int[] inputRows, int[] uRows, int diag, int[] lRows)
		{
			for (int p = 0; p < inputRows.length; p++)
			{
				int row = inputRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			for (int p = 0; p < uRows.length; p++)
			{
				int row = uRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			xr[diag] = 0.0;
			xi[diag] = 0.0;
			for (int p = 0; p < lRows.length; p++)
			{
				int row = lRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
		}

		private static void clearRefactorColumn(double[] xr, double[] xi,
				int[] uRows, int diag, int[] lRows)
		{
			for (int p = 0; p < uRows.length; p++)
			{
				int row = uRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			xr[diag] = 0.0;
			xi[diag] = 0.0;
			for (int p = 0; p < lRows.length; p++)
			{
				int row = lRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
		}

		private static void clearRefactorColumn(double[] xr, double[] xi,
				int[] uRows, int diag, int[] lRows, int uLen, int lLen)
		{
			for (int p = 0; p < uLen; p++)
			{
				int row = uRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			xr[diag] = 0.0;
			xi[diag] = 0.0;
			for (int p = 0; p < lLen; p++)
			{
				int row = lRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
		}

		private static void clearRefactorColumn(double[] xr, double[] xi,
				int[] inputRows, int[] uRows, int uStart, int uEnd, int diag,
				int[] lRows, int lStart, int lEnd)
		{
			for (int p = 0; p < inputRows.length; p++)
			{
				int row = inputRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			for (int p = uStart; p < uEnd; p++)
			{
				int row = uRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			xr[diag] = 0.0;
			xi[diag] = 0.0;
			for (int p = lStart; p < lEnd; p++)
			{
				int row = lRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
		}

		private static void clearRefactorColumn(double[] xr, double[] xi,
				int[] uRows, int uStart, int uEnd, int diag,
				int[] lRows, int lStart, int lEnd)
		{
			for (int p = uStart; p < uEnd; p++)
			{
				int row = uRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			xr[diag] = 0.0;
			xi[diag] = 0.0;
			for (int p = lStart; p < lEnd; p++)
			{
				int row = lRows[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
		}

		private static void clearWork(double[] xr, double[] xi, IntList touched,
				ActiveList active)
		{
			for (int p = 0; p < touched.size; p++)
			{
				int row = touched.values[p];
				xr[row] = 0.0;
				xi[row] = 0.0;
			}
			touched.clear();
			active.clear();
		}

		private static boolean isZero(double re, double im)
		{
			return Math.abs(re) <= DROP_TOL && Math.abs(im) <= DROP_TOL;
		}

		private static boolean allTrue(boolean[] values)
		{
			for (int i = 0; i < values.length; i++)
			{
				if (!values[i])
				{
					return false;
				}
			}
			return true;
		}

	private static void ensureSymbolicState(KLU_z_numeric Numeric)
	{
		int n = Numeric.n;
		if (Numeric.RowPerm == null)
		{
			Numeric.RowPerm = identity(n);
			Numeric.RowInv = identity(n);
		}
		if (Numeric.ColPerm == null)
		{
			Numeric.ColPerm = identity(n);
			Numeric.ColInv = identity(n);
		}
		if (Numeric.RowInv == null)
		{
			Numeric.RowInv = inverse(Numeric.RowPerm);
		}
		if (Numeric.ColInv == null)
		{
			Numeric.ColInv = inverse(Numeric.ColPerm);
		}
		if (Numeric.R == null)
		{
			Numeric.nblocks = 1;
			Numeric.R = naturalR(n);
		}
		if (Numeric.nblocks <= 0)
		{
			Numeric.nblocks = Numeric.R.length - 1;
		}
		if (Numeric.BlockOf == null)
		{
			Numeric.BlockOf = blockOf(n, Numeric.R, Numeric.nblocks);
		}
	}

	private static SparseRow[] buildSparseRows(KLU_z_numeric Numeric, KLU_common Common)
	{
		int n = Numeric.n;
		SparseRow[] rows = new SparseRow[n];
		for (int row = 0; row < n; row++)
		{
			rows[row] = new SparseRow();
		}
		int[] W = new int[n];
		if (Numeric.Rs != null)
		{
			if (klu_z_scale(Common.scale, n, Numeric.Ap, Numeric.Ai, Numeric.Ax,
					Numeric.Rs, W, Common) == FALSE)
			{
				return null;
			}
		}
		for (int col = 0; col < n; col++)
		{
			for (int p = Numeric.Ap[col]; p < Numeric.Ap[col + 1]; p++)
			{
					int oldRow = Numeric.Ai[p];
					if (oldRow < 0 || oldRow >= n)
					{
						Common.status = KLU_INVALID;
						return null;
					}
					int row = Numeric.RowInv == null ? oldRow : Numeric.RowInv[oldRow];
					int newCol = Numeric.ColInv == null ? col : Numeric.ColInv[col];
					double scale = Numeric.Rs == null ? 1.0 : Numeric.Rs[oldRow];
					rows[row].add(newCol, Numeric.Ax[2 * p] / scale,
							Numeric.Ax[2 * p + 1] / scale);
				}
			}
		return rows;
	}

	@SuppressWarnings("unchecked")
	private static TreeSet<Integer>[] buildColumnRows(SparseRow[] rows, int n)
	{
		TreeSet<Integer>[] columnRows = new TreeSet[n];
		for (int col = 0; col < n; col++)
		{
			columnRows[col] = new TreeSet<Integer>();
		}
		for (int row = 0; row < n; row++)
		{
			for (Integer col : rows[row].keySet())
			{
				columnRows[col.intValue()].add(Integer.valueOf(row));
			}
		}
		return columnRows;
	}

	private static void extractOffDiagonalBlocks(KLU_z_numeric Numeric, SparseRow[] rows)
	{
		int n = Numeric.n;
		int[] offCounts = new int[n];
		int nzoff = 0;
		for (int row = 0; row < n; row++)
		{
			List<Integer> offCols = new ArrayList<Integer>();
			for (Integer colRef : rows[row].keySet())
			{
				int col = colRef.intValue();
				if (Numeric.BlockOf[row] != Numeric.BlockOf[col])
				{
					offCols.add(colRef);
					offCounts[col]++;
					nzoff++;
				}
			}
			for (Integer colRef : offCols)
			{
				rows[row].moveToOff(colRef.intValue());
			}
		}

		Numeric.nzoff = nzoff;
		Numeric.Offp = new int[n + 1];
		for (int col = 0; col < n; col++)
		{
			Numeric.Offp[col + 1] = Numeric.Offp[col] + offCounts[col];
		}
		Numeric.Offi = new int[nzoff];
		Numeric.Offx = new double[2 * nzoff];
		int[] next = copy(Numeric.Offp, n);
		for (int col = 0; col < n; col++)
		{
			next[col] = Numeric.Offp[col];
		}

		for (int row = 0; row < n; row++)
		{
			for (Map.Entry<Integer, Complex> entry : rows[row].offEntries())
			{
				int col = entry.getKey().intValue();
				Complex value = entry.getValue();
				int p = next[col]++;
				Numeric.Offi[p] = row;
				Numeric.Offx[2 * p] = value.re;
				Numeric.Offx[2 * p + 1] = value.im;
			}
			rows[row].clearOffEntries();
		}
	}

	private static void swapRows(SparseRow[] rows, TreeSet<Integer>[] columnRows,
			int row1, int row2)
	{
		if (row1 == row2)
		{
			return;
		}
		SparseRow first = rows[row1];
		SparseRow second = rows[row2];
		TreeSet<Integer> cols = new TreeSet<Integer>();
		for (Integer col : first.keySet())
		{
			cols.add(col);
		}
		for (Integer col : second.keySet())
		{
			cols.add(col);
		}
		for (Integer colRef : cols)
		{
			int col = colRef.intValue();
			boolean firstHas = first.get(col) != null;
			boolean secondHas = second.get(col) != null;
			columnRows[col].remove(Integer.valueOf(row1));
			columnRows[col].remove(Integer.valueOf(row2));
			if (secondHas)
			{
				columnRows[col].add(Integer.valueOf(row1));
			}
			if (firstHas)
			{
				columnRows[col].add(Integer.valueOf(row2));
			}
		}
		rows[row1] = second;
		rows[row2] = first;
	}

	private static void populateSparsePatternState(KLU_z_numeric Numeric, SparseRow[] rows)
	{
			int n = Numeric.n;
			Numeric.Lip = new int[n];
		Numeric.Uip = new int[n];
		Numeric.Llen = new int[n];
		Numeric.Ulen = new int[n];
		Numeric.LUsize = new int[] {0};
		Numeric.LUbx = new double[][] {new double[0]};

			int lnz = 0;
			int unz = 0;
			Numeric.LUrowCols = new int[n][];
			Numeric.LUrowValues = new double[n][];
			Numeric.LrowCols = new int[n][];
			Numeric.LrowValues = new double[n][];
			Numeric.UrowCols = new int[n][];
			Numeric.UrowValues = new double[n][];
			for (int row = 0; row < n; row++)
			{
				SparseRow sparseRow = rows[row];
				int size = sparseRow.size();
				int[] cols = new int[size];
				double[] values = new double[2 * size];
				int lower = 0;
				int upper = 0;
				for (Map.Entry<Integer, Complex> entry : sparseRow.entrySet())
				{
					int col = entry.getKey().intValue();
					if (col < row)
					{
						lower++;
					}
					else if (col > row)
					{
						upper++;
					}
				}
				int[] lcols = new int[lower];
				double[] lvalues = new double[2 * lower];
				int[] ucols = new int[upper];
				double[] uvalues = new double[2 * upper];
				int p = 0;
				int lp = 0;
				int up = 0;
				for (Map.Entry<Integer, Complex> entry : sparseRow.entrySet())
				{
					int col = entry.getKey().intValue();
					Complex value = entry.getValue();
					cols[p] = col;
				values[2 * p] = value.re;
				values[2 * p + 1] = value.im;
					if (col < row)
					{
						lnz++;
						lcols[lp] = col;
						lvalues[2 * lp] = value.re;
						lvalues[2 * lp + 1] = value.im;
						lp++;
					}
					else if (col > row)
					{
						unz++;
						ucols[up] = col;
						uvalues[2 * up] = value.re;
						uvalues[2 * up + 1] = value.im;
						up++;
					}
					p++;
				}
				Numeric.LUrowCols[row] = cols;
				Numeric.LUrowValues[row] = values;
				Numeric.LrowCols[row] = lcols;
				Numeric.LrowValues[row] = lvalues;
				Numeric.UrowCols[row] = ucols;
				Numeric.UrowValues[row] = uvalues;
			}
			Numeric.lnz = lnz + n;
			Numeric.unz = unz + n;
				Numeric.Work = new double[8 * Math.max(1, n)];
				Numeric.Xwork = Numeric.Work;
				Numeric.LUrowsCurrent = true;
			}

		private static void populateSparsePatternState(KLU_z_numeric Numeric, RowBuilder[] rows)
		{
			int n = Numeric.n;
			Numeric.Lip = new int[n];
			Numeric.Uip = new int[n];
			Numeric.Llen = new int[n];
			Numeric.Ulen = new int[n];
			Numeric.LUsize = new int[] {0};
			Numeric.LUbx = new double[][] {new double[0]};

			int lnz = 0;
			int unz = 0;
			boolean splitRows = Boolean.getBoolean("jklu.complex.solve.splitRows");
			Numeric.LUrowCols = new int[n][];
			Numeric.LUrowValues = new double[n][];
			Numeric.LrowCols = splitRows ? new int[n][] : null;
			Numeric.LrowValues = splitRows ? new double[n][] : null;
			Numeric.UrowCols = splitRows ? new int[n][] : null;
			Numeric.UrowValues = splitRows ? new double[n][] : null;
			for (int row = 0; row < n; row++)
			{
				RowBuilder builder = rows[row];
				int[] cols = new int[builder.size];
				double[] values = new double[2 * builder.size];
				System.arraycopy(builder.cols, 0, cols, 0, builder.size);
				System.arraycopy(builder.values, 0, values, 0, 2 * builder.size);
				for (int p = 0; p < builder.size; p++)
				{
					int col = cols[p];
					if (col < row)
					{
						lnz++;
					}
					else if (col > row)
					{
						unz++;
					}
				}
				if (splitRows)
				{
					int lower = 0;
					int upper = 0;
					for (int p = 0; p < builder.size; p++)
					{
						int col = cols[p];
						if (col < row)
						{
							lower++;
						}
						else if (col > row)
						{
							upper++;
						}
					}
					int[] lcols = new int[lower];
					double[] lvalues = new double[2 * lower];
					int[] ucols = new int[upper];
					double[] uvalues = new double[2 * upper];
					int lp = 0;
					int up = 0;
					for (int p = 0; p < builder.size; p++)
					{
						int col = cols[p];
						if (col < row)
						{
							lcols[lp] = col;
							lvalues[2 * lp] = values[2 * p];
							lvalues[2 * lp + 1] = values[2 * p + 1];
							lp++;
						}
						else if (col > row)
						{
							ucols[up] = col;
							uvalues[2 * up] = values[2 * p];
							uvalues[2 * up + 1] = values[2 * p + 1];
							up++;
						}
					}
					Numeric.LrowCols[row] = lcols;
					Numeric.LrowValues[row] = lvalues;
					Numeric.UrowCols[row] = ucols;
					Numeric.UrowValues[row] = uvalues;
				}
				Numeric.LUrowCols[row] = cols;
				Numeric.LUrowValues[row] = values;
			}
				Numeric.lnz = lnz + n;
				Numeric.unz = unz + n;
				Numeric.Work = new double[8 * Math.max(1, n)];
				Numeric.Xwork = Numeric.Work;
				Numeric.LUrowsCurrent = true;
			}

	/**
	 * Factor a complex sparse matrix in compressed-column form.
	 *
	 * Ax stores interleaved complex values: [real0, imag0, real1, imag1, ...].
	 */
	public static KLU_z_numeric klu_z_factor(int[] Ap, int[] Ai, double[] Ax,
			KLU_symbolic Symbolic, KLU_common Common)
	{
		if (!validInputs(Ap, Ai, Ax, Symbolic, Common))
		{
			return null;
		}

		KLU_z_numeric Numeric = new KLU_z_numeric();
			Numeric.n = Symbolic.n;
			Numeric.nz = Ap[Symbolic.n];
		Numeric.Ap = copy(Ap, Symbolic.n + 1);
			Numeric.ApRef = Ap;
			Numeric.Ai = copy(Ai, Numeric.nz);
			Numeric.AiRef = Ai;
			Numeric.Ax = Ax;
			Numeric.Rs = Common.scale > 0 ? new double[Numeric.n] : null;
			Numeric.RowPerm = copy(Symbolic.P, Symbolic.n);
			Numeric.ColPerm = copy(Symbolic.Q, Symbolic.n);
			Numeric.RowInv = inverse(Numeric.RowPerm);
			Numeric.ColInv = inverse(Numeric.ColPerm);
			Numeric.nblocks = Math.max(1, Symbolic.nblocks);
			Numeric.R = Symbolic.R == null ? naturalR(Symbolic.n) : copy(Symbolic.R, Numeric.nblocks + 1);
			Numeric.BlockOf = blockOf(Symbolic.n, Numeric.R, Numeric.nblocks);
			Numeric.nzoff = Math.max(0, Symbolic.nzoff);

			return factorNumeric(Numeric, Common) == TRUE ? Numeric : null;
		}


	static double real(KLU_z_numeric Numeric, int row, int col)
	{
		int p = position(Numeric, row, col);
		return p < 0 ? 0.0 : Numeric.LUrowValues[row][2 * p];
	}

	static double imag(KLU_z_numeric Numeric, int row, int col)
	{
		int p = position(Numeric, row, col);
		return p < 0 ? 0.0 : Numeric.LUrowValues[row][2 * p + 1];
	}

		static double abs(KLU_z_numeric Numeric, int row, int col)
		{
			return Zklu_complex.abs(real(Numeric, row, col), imag(Numeric, row, col));
		}

		static int rowLength(KLU_z_numeric Numeric, int row)
		{
			return Numeric.LUrowCols[row].length;
		}

		static int rowCol(KLU_z_numeric Numeric, int row, int p)
		{
			return Numeric.LUrowCols[row][p];
		}

		static double rowReal(KLU_z_numeric Numeric, int row, int p)
		{
			return Numeric.LUrowValues[row][2 * p];
		}

		static double rowImag(KLU_z_numeric Numeric, int row, int p)
		{
			return Numeric.LUrowValues[row][2 * p + 1];
		}

	private static int position(KLU_z_numeric Numeric, int row, int col)
	{
		int[] cols = Numeric.LUrowCols[row];
		int lo = 0;
		int hi = cols.length - 1;
		while (lo <= hi)
		{
			int mid = (lo + hi) >>> 1;
			int c = cols[mid];
			if (c == col)
			{
				return mid;
			}
			if (c < col)
			{
				lo = mid + 1;
			}
			else
			{
				hi = mid - 1;
			}
			}
			return -1;
		}

		private static final class Column
		{
			final int size;
			final int[] index;
			final double[] real;
			final double[] imag;

			Column(int size, int[] index, double[] real, double[] imag)
			{
				this.size = size;
				this.index = index;
				this.real = real;
				this.imag = imag;
			}
		}

		private static final class DirectColumns
		{
				final int[][] rows;
				final double[][] real;
				final double[][] imag;
				final int[] p;
				final int[] i;
				final double[] x;
				final KLU_z_numeric numeric;
				final boolean cscView;

				DirectColumns(int[][] rows, double[][] real, double[][] imag)
				{
					this.rows = rows;
					this.real = real;
					this.imag = imag;
					this.p = null;
					this.i = null;
					this.x = null;
					this.numeric = null;
					this.cscView = false;
				}

				DirectColumns(int[] p, int[] i, double[] x)
				{
					this.rows = null;
					this.real = null;
					this.imag = null;
					this.p = p;
					this.i = i;
					this.x = x;
					this.numeric = null;
					this.cscView = false;
				}

				DirectColumns(KLU_z_numeric numeric)
				{
					this.rows = null;
					this.real = null;
					this.imag = null;
					this.p = numeric.Ap;
					this.i = numeric.Ai;
					this.x = numeric.Ax;
					this.numeric = numeric;
					this.cscView = true;
				}

				int start(int col)
				{
					if (cscView)
					{
						int oldCol = numeric.ColPerm == null ? col : numeric.ColPerm[col];
						return p[oldCol];
					}
					return p == null ? 0 : p[col];
				}

				int end(int col)
				{
					if (cscView)
					{
						int oldCol = numeric.ColPerm == null ? col : numeric.ColPerm[col];
						return p[oldCol + 1];
					}
					return p == null ? rows[col].length : p[col + 1];
				}

			int count(int col)
			{
				return end(col) - start(col);
			}

				int row(int col, int pos)
				{
					if (cscView)
					{
						int oldRow = i[pos];
						return numeric.RowInv == null ? oldRow : numeric.RowInv[oldRow];
					}
					return p == null ? rows[col][pos] : i[pos];
				}

				double real(int col, int pos)
				{
					if (cscView)
					{
						double scale = numeric.Rs == null ? 1.0 : numeric.Rs[i[pos]];
						return x[2 * pos] / scale;
					}
					return p == null ? real[col][pos] : x[2 * pos];
				}

				double imag(int col, int pos)
				{
					if (cscView)
					{
						double scale = numeric.Rs == null ? 1.0 : numeric.Rs[i[pos]];
						return x[2 * pos + 1] / scale;
					}
					return p == null ? imag[col][pos] : x[2 * pos + 1];
				}

					void clear(double[] xr, double[] xi, int col)
			{
				for (int pos = start(col); pos < end(col); pos++)
				{
					int row = row(col, pos);
					xr[row] = 0.0;
					xi[row] = 0.0;
				}
			}
		}

			private static final class ColumnKernelWork
			{
				final int n;
				final FlatColumnStore lStore;
				final FlatColumnStore uStore;
				final RowBuilder[] luRows;
			final double[] xr;
			final double[] xi;
			final int[] mark;
			final int[] reachMark;
			final int[] reachStack;
			final int[] reachPos;
			final int[] reachRows;
			final int[] lpend;
			final IntList touched;
			final ActiveList active;
			final int[] pivots;
			int generation;
			int reachGeneration;
			int reachCount;

			ColumnKernelWork(int n)
			{
					this.n = n;
					this.lStore = new FlatColumnStore(n);
					this.uStore = new FlatColumnStore(n);
					this.luRows = new RowBuilder[n];
					for (int i = 0; i < n; i++)
					{
						luRows[i] = new RowBuilder();
					}
					this.xr = new double[n];
					this.xi = new double[n];
					this.mark = new int[n];
					this.reachMark = new int[n];
					this.reachStack = new int[n];
					this.reachPos = new int[n];
					this.reachRows = new int[n];
					this.lpend = new int[n];
					this.touched = new IntList();
					this.active = new ActiveList();
					this.pivots = new int[n];
				}

			int nextGeneration()
			{
				generation++;
				if (generation == 0)
				{
					for (int i = 0; i < mark.length; i++)
					{
						mark[i] = 0;
					}
					generation = 1;
				}
				return generation;
			}

			int nextReachGeneration()
			{
				reachGeneration++;
				if (reachGeneration == 0)
				{
					for (int i = 0; i < reachMark.length; i++)
					{
						reachMark[i] = 0;
					}
					reachGeneration = 1;
				}
				return reachGeneration;
			}

			void clearBuilders(boolean clearRows)
			{
					for (int i = 0; i < n; i++)
					{
						if (clearRows)
					{
						luRows[i].clear();
					}
				}
			}
		}

			private static final class ColumnBuilder
			{
				int size;
				int[] index = new int[4];
				double[] real = new double[4];
			double[] imag = new double[4];

			void add(int row, double re, double im)
			{
				if (isZero(re, im))
				{
					return;
				}
				addNonZero(row, re, im);
			}

			void addNonZero(int row, double re, double im)
			{
				ensure(size + 1);
				index[size] = row;
				real[size] = re;
				imag[size] = im;
				size++;
			}

				Column toColumn()
				{
					int[] columnIndex = new int[size];
				double[] columnReal = new double[size];
				double[] columnImag = new double[size];
				System.arraycopy(index, 0, columnIndex, 0, size);
				System.arraycopy(real, 0, columnReal, 0, size);
				System.arraycopy(imag, 0, columnImag, 0, size);
					return new Column(size, columnIndex, columnReal, columnImag);
				}

				void clear()
				{
					size = 0;
				}

				private void ensure(int capacity)
			{
				if (capacity <= index.length)
				{
					return;
				}
				int newCapacity = Math.max(capacity, index.length * 2);
				int[] newIndex = new int[newCapacity];
				double[] newReal = new double[newCapacity];
				double[] newImag = new double[newCapacity];
				System.arraycopy(index, 0, newIndex, 0, size);
				System.arraycopy(real, 0, newReal, 0, size);
				System.arraycopy(imag, 0, newImag, 0, size);
				index = newIndex;
				real = newReal;
				imag = newImag;
				}
			}

				private static final class FlatColumnStore
				{
					final int[] p;
					final int[] len;
				int size;
				int currentCol;
				int[] index = new int[1024];
				double[] real = new double[1024];
				double[] imag = new double[1024];

				FlatColumnStore(int n)
					{
						this.p = new int[n + 1];
						this.len = new int[n];
					}

				void clear()
				{
					size = 0;
						for (int i = 0; i < len.length; i++)
						{
							len[i] = 0;
						}
				}

					void startColumn(int col)
						{
							currentCol = col;
							p[col] = size;
							len[col] = 0;
						}

					void finish(int n)
					{
						p[n] = size;
					}

					void add(int row, double re, double im)
				{
					ensure(size + 1);
					index[size] = row;
					real[size] = re;
					imag[size] = im;
						size++;
						len[currentCol]++;
					}

				void addFast(int row, double re, double im)
				{
					if (size == index.length)
					{
						ensureCapacity(size + 1);
					}
					index[size] = row;
					real[size] = re;
					imag[size] = im;
					size++;
					len[currentCol]++;
				}

				void ensureCapacity(int capacity)
				{
					if (capacity <= index.length)
					{
						return;
					}
					int newCapacity = Math.max(capacity, index.length * 2);
					int[] newIndex = new int[newCapacity];
					double[] newReal = new double[newCapacity];
					double[] newImag = new double[newCapacity];
					System.arraycopy(index, 0, newIndex, 0, size);
					System.arraycopy(real, 0, newReal, 0, size);
					System.arraycopy(imag, 0, newImag, 0, size);
					index = newIndex;
					real = newReal;
					imag = newImag;
				}

					private void ensure(int capacity)
					{
						ensureCapacity(capacity);
					}

				}

			private static final class RowBuilder
		{
			int size;
			int[] cols = new int[4];
			double[] values = new double[8];

			void add(int col, double re, double im)
			{
				if (isZero(re, im))
				{
					return;
				}
				ensure(size + 1);
				cols[size] = col;
				values[2 * size] = re;
				values[2 * size + 1] = im;
				size++;
			}

			void clear()
			{
				size = 0;
			}

			private void ensure(int capacity)
			{
				if (capacity <= cols.length)
				{
					return;
				}
				int newCapacity = Math.max(capacity, cols.length * 2);
				int[] newCols = new int[newCapacity];
				double[] newValues = new double[2 * newCapacity];
				System.arraycopy(cols, 0, newCols, 0, size);
				System.arraycopy(values, 0, newValues, 0, 2 * size);
				cols = newCols;
				values = newValues;
			}
		}

		private static final class IntList
		{
			int size;
			int[] values = new int[16];

				void add(int value)
				{
					ensureCapacity(size + 1);
					values[size++] = value;
				}

				void ensureCapacity(int capacity)
				{
					if (capacity <= values.length)
					{
						return;
					}
					int[] next = new int[Math.max(capacity, values.length * 2)];
					System.arraycopy(values, 0, next, 0, size);
					values = next;
				}

			void clear()
			{
				size = 0;
			}
		}

		private static final class ActiveList
		{
			int size;
			int[] values = new int[16];
			long addCount;
			long pollCount;

				void add(int value)
				{
					addCount++;
					ensureCapacity(size + 1);
					int pos = size++;
				while (pos > 0)
				{
					int parent = (pos - 1) >>> 1;
					int parentValue = values[parent];
					if (parentValue <= value)
					{
						break;
					}
					values[pos] = parentValue;
					pos = parent;
				}
				values[pos] = value;
			}

			int pollMinBelow(int limit)
			{
				if (size == 0)
				{
					return EMPTY;
				}
				int best = values[0];
				if (best >= limit)
				{
					return EMPTY;
				}
				int replacement = values[--size];
				if (size > 0)
				{
					int pos = 0;
					while (true)
					{
						int child = (pos << 1) + 1;
						if (child >= size)
						{
							break;
						}
						int right = child + 1;
						if (right < size && values[right] < values[child])
						{
							child = right;
						}
						if (values[child] >= replacement)
						{
							break;
						}
						values[pos] = values[child];
						pos = child;
					}
					values[pos] = replacement;
				}
				pollCount++;
				return best;
			}

			void resetStats()
			{
				addCount = 0L;
				pollCount = 0L;
			}

			void clear()
			{
				size = 0;
			}

				void ensureCapacity(int capacity)
				{
					if (capacity <= values.length)
					{
					return;
				}
				int[] next = new int[Math.max(capacity, values.length * 2)];
				System.arraycopy(values, 0, next, 0, size);
				values = next;
			}
		}

				static int[] copy(int[] source, int length)
				{
					int[] copy = new int[length];
					System.arraycopy(source, 0, copy, 0, length);
					return copy;
				}

				static int[] copy(int[] source, int offset, int length)
				{
					int[] copy = new int[length];
					System.arraycopy(source, offset, copy, 0, length);
					return copy;
				}

				static double[] copy(double[] source, int length)
				{
					double[] copy = new double[length];
					System.arraycopy(source, 0, copy, 0, length);
					return copy;
				}

				static double[] copy(double[] source, int offset, int length)
				{
					double[] copy = new double[length];
					System.arraycopy(source, offset, copy, 0, length);
					return copy;
				}

			private static int[] inverse(int[] permutation)
		{
			int[] inverse = new int[permutation.length];
			for (int i = 0; i < permutation.length; i++)
			{
				inverse[permutation[i]] = i;
			}
			return inverse;
		}

		private static int[] naturalR(int n)
		{
			return new int[] {0, n};
		}

		private static int[] identity(int n)
		{
			int[] identity = new int[n];
			for (int i = 0; i < n; i++)
			{
				identity[i] = i;
			}
			return identity;
		}

		private static int[] blockOf(int n, int[] R, int nblocks)
		{
			int[] blockOf = new int[n];
			for (int block = 0; block < nblocks; block++)
			{
				for (int k = R[block]; k < R[block + 1]; k++)
				{
					blockOf[k] = block;
				}
			}
			return blockOf;
		}

	private static final class SparseRow
	{
		private final TreeMap<Integer, Complex> values = new TreeMap<Integer, Complex>();
		private final TreeMap<Integer, Complex> offValues = new TreeMap<Integer, Complex>();

		Complex get(int col)
		{
			return values.get(Integer.valueOf(col));
		}

		double abs(int col)
		{
			Complex value = get(col);
			return value == null ? 0.0 : value.abs();
		}

		void add(int col, double re, double im)
		{
			Complex value = get(col);
			if (value == null)
			{
				if (re != 0.0 || im != 0.0)
				{
					values.put(Integer.valueOf(col), new Complex(re, im));
				}
				return;
			}
			put(col, value.re + re, value.im + im);
		}

		boolean put(int col, double re, double im)
		{
			if (Math.abs(re) <= DROP_TOL && Math.abs(im) <= DROP_TOL)
			{
				values.remove(Integer.valueOf(col));
				return false;
			}
			else
			{
				values.put(Integer.valueOf(col), new Complex(re, im));
				return true;
			}
		}

		int subtractProduct(int col, Complex a, Complex b)
		{
			Complex value = get(col);
			boolean hadValue = value != null;
			double re = value == null ? 0.0 : value.re;
			double im = value == null ? 0.0 : value.im;
			boolean hasValue = put(col, re - (a.re * b.re - a.im * b.im),
					im - (a.im * b.re + a.re * b.im));
			if (!hadValue && hasValue)
			{
				return 1;
			}
			if (hadValue && !hasValue)
			{
				return -1;
			}
			return 0;
		}

		int size()
		{
			return values.size();
		}

		Iterable<Map.Entry<Integer, Complex>> entrySet()
		{
			return values.entrySet();
		}

		Iterable<Integer> keySet()
		{
			return values.keySet();
		}

			Iterable<Map.Entry<Integer, Complex>> tailMap(int col)
			{
				return values.tailMap(Integer.valueOf(col)).entrySet();
			}

			void moveToOff(int col)
			{
				Complex value = values.remove(Integer.valueOf(col));
				if (value != null)
				{
					offValues.put(Integer.valueOf(col), value);
				}
			}

			Iterable<Map.Entry<Integer, Complex>> offEntries()
			{
				return offValues.entrySet();
			}

			void clearOffEntries()
			{
				offValues.clear();
			}
		}

	private static final class FactorProfile
	{
		private final boolean enabled;
		private final int interval;
		private final int maxPivots;
		private final int n;
		private final long startNs;
		private long candidateCount;
		private long rowsWithPivotCount;
		private long updateCount;
		private long directScatterNs;
		private long directUpdateNs;
		private long directPivotNs;
		private long directScaleNs;
		private long directClearNs;
		private long directFinalizeNs;
		private long directInputNs;
		private long directPriorUpdateTargets;
		private long directFutureUpdateTargets;
		private long directPriorNewTouches;
			private long directFutureNewTouches;
			private long directActiveAdds;
			private long directActivePolls;
			private long directRealUUpdates;
			private long directRealLUpdates;
			private long directPureRealUpdates;
			private long directUpdateColumns;
			private long directTinyUpdateColumns;
			private long directSmallUpdateColumns;
			private long directMediumUpdateColumns;
			private long directLargeUpdateColumns;
			private long directTinyUpdateEntries;
			private long directSmallUpdateEntries;
			private long directMediumUpdateEntries;
			private long directLargeUpdateEntries;
			private int maxRowNnz;

		private FactorProfile(boolean enabled, int interval, int maxPivots, int n)
		{
			this.enabled = enabled;
			this.interval = interval <= 0 ? 1000 : interval;
			this.maxPivots = maxPivots;
			this.n = n;
			this.startNs = System.nanoTime();
		}

		static FactorProfile create(int n)
		{
			boolean enabled = Boolean.getBoolean("jklu.profile.factor");
			int interval = intProperty("jklu.profile.factor.interval", 1000);
			int maxPivots = intProperty("jklu.profile.factor.maxPivots", -1);
			return new FactorProfile(enabled, interval, maxPivots, n);
		}

		boolean enabled()
		{
			return enabled;
		}

		void directPhases(long scatterNs, long updateNs, long pivotNs,
				long scaleNs, long clearNs)
		{
			if (!enabled)
			{
				return;
			}
			directScatterNs += scatterNs;
			directUpdateNs += updateNs;
			directPivotNs += pivotNs;
			directScaleNs += scaleNs;
			directClearNs += clearNs;
		}

		void directFinalize(long finalizeNs)
		{
			if (!enabled)
			{
				return;
			}
			directFinalizeNs += finalizeNs;
		}

		void directInput(long inputNs, boolean directInputAvailable)
		{
			if (!enabled)
			{
				return;
			}
			directInputNs += inputNs;
			System.out.println("factorProfileDirectInput n=" + n +
					",available=" + directInputAvailable +
					",directInputMs=" + millis(inputNs));
			System.out.flush();
		}

		void directActive(long addCount, long pollCount)
		{
			if (!enabled)
			{
				return;
			}
			directActiveAdds += addCount;
			directActivePolls += pollCount;
		}

		void updateTarget(boolean prior)
		{
			if (!enabled)
			{
				return;
			}
			if (prior)
			{
				directPriorUpdateTargets++;
			}
			else
			{
				directFutureUpdateTargets++;
			}
		}

			void newTouch(boolean prior)
			{
				if (!enabled)
				{
					return;
			}
			if (prior)
			{
				directPriorNewTouches++;
			}
			else
			{
				directFutureNewTouches++;
			}
		}

		void printSetup(long buildRowsNs, long buildColumnRowsNs, SparseRow[] rows)
		{
			if (!enabled)
			{
				return;
			}
			maxRowNnz = maxRowNnz(rows);
			System.out.println("factorProfileSetup n=" + n +
					",buildRowsMs=" + millis(buildRowsNs) +
					",buildColumnRowsMs=" + millis(buildColumnRowsNs) +
					",initialMaxRowNnz=" + maxRowNnz);
			System.out.flush();
		}

			void step(int k, int pivotCandidates, int rowsWithPivot, int pivotRowNnz,
					int updates, SparseRow[] rows)
		{
			if (!enabled)
			{
				return;
			}
			candidateCount += pivotCandidates;
			rowsWithPivotCount += rowsWithPivot;
			updateCount += updates;
			if ((k + 1) % interval == 0 || k + 1 == n || shouldStop(k))
			{
				maxRowNnz = Math.max(maxRowNnz, maxRowNnz(rows));
				System.out.println("factorProfileStep pivot=" + (k + 1) + "/" + n +
						",elapsedMs=" + millis(System.nanoTime() - startNs) +
						",pivotCandidates=" + pivotCandidates +
						",rowsWithPivot=" + rowsWithPivot +
						",pivotRowNnz=" + pivotRowNnz +
						",updates=" + updates +
						",cumCandidates=" + candidateCount +
						",cumRowsWithPivot=" + rowsWithPivotCount +
						",cumUpdates=" + updateCount +
						",maxRowNnz=" + maxRowNnz +
						directPhaseSummary());
				System.out.flush();
				}
			}

			void updateShape(double ur, double ui, double lr, double li)
			{
				if (!enabled)
				{
					return;
				}
				boolean realU = ui == 0.0;
				boolean realL = li == 0.0;
				if (realU)
				{
					directRealUUpdates++;
				}
				if (realL)
				{
					directRealLUpdates++;
				}
				if (realU && realL)
				{
					directPureRealUpdates++;
				}
			}

			void updateColumnLength(int length)
			{
				if (!enabled)
				{
					return;
				}
				directUpdateColumns++;
				if (length <= 4)
				{
					directTinyUpdateColumns++;
					directTinyUpdateEntries += length;
				}
				else if (length <= 16)
				{
					directSmallUpdateColumns++;
					directSmallUpdateEntries += length;
				}
				else if (length <= 64)
				{
					directMediumUpdateColumns++;
					directMediumUpdateEntries += length;
				}
				else
				{
					directLargeUpdateColumns++;
					directLargeUpdateEntries += length;
				}
			}

			void step(int k, int pivotCandidates, int rowsWithPivot, int pivotRowNnz,
					int updates, RowBuilder[] rows)
			{
				if (!enabled)
				{
					return;
				}
				candidateCount += pivotCandidates;
				rowsWithPivotCount += rowsWithPivot;
				updateCount += updates;
				if ((k + 1) % interval == 0 || k + 1 == n || shouldStop(k))
				{
					maxRowNnz = Math.max(maxRowNnz, maxRowNnz(rows));
					System.out.println("factorProfileStep pivot=" + (k + 1) + "/" + n +
							",elapsedMs=" + millis(System.nanoTime() - startNs) +
							",pivotCandidates=" + pivotCandidates +
							",rowsWithPivot=" + rowsWithPivot +
							",pivotRowNnz=" + pivotRowNnz +
							",updates=" + updates +
							",cumCandidates=" + candidateCount +
							",cumRowsWithPivot=" + rowsWithPivotCount +
							",cumUpdates=" + updateCount +
								",maxRowNnz=" + maxRowNnz +
								directPhaseSummary());
					System.out.flush();
				}
			}

		boolean shouldStop(int k)
		{
			return enabled && maxPivots > 0 && k + 1 >= maxPivots;
		}

			void finish(SparseRow[] rows)
		{
			if (!enabled)
			{
				return;
			}
			System.out.println("factorProfileFinish elapsedMs=" +
					millis(System.nanoTime() - startNs) +
					",cumCandidates=" + candidateCount +
					",cumRowsWithPivot=" + rowsWithPivotCount +
					",cumUpdates=" + updateCount +
					",finalMaxRowNnz=" + maxRowNnz(rows) +
					directPhaseSummary());
				System.out.flush();
			}

			void finish(RowBuilder[] rows)
			{
				if (!enabled)
				{
					return;
				}
				System.out.println("factorProfileFinish elapsedMs=" +
						millis(System.nanoTime() - startNs) +
						",cumCandidates=" + candidateCount +
						",cumRowsWithPivot=" + rowsWithPivotCount +
						",cumUpdates=" + updateCount +
							",finalMaxRowNnz=" + maxRowNnz(rows) +
							directPhaseSummary());
				System.out.flush();
			}

		private static int intProperty(String name, int defaultValue)
		{
			Integer value = Integer.getInteger(name);
			return value == null ? defaultValue : value.intValue();
		}

			private static int maxRowNnz(SparseRow[] rows)
		{
			int max = 0;
			for (SparseRow row : rows)
			{
				max = Math.max(max, row.size());
			}
				return max;
			}

			private static int maxRowNnz(RowBuilder[] rows)
			{
				int max = 0;
				for (RowBuilder row : rows)
				{
					max = Math.max(max, row.size);
				}
				return max;
			}

		private static double millis(long nanos)
		{
			return nanos / 1.0e6;
		}

		private String directPhaseSummary()
		{
			if (directInputNs == 0L && directScatterNs == 0L && directUpdateNs == 0L &&
					directPivotNs == 0L && directScaleNs == 0L &&
					directClearNs == 0L && directFinalizeNs == 0L)
			{
				return "";
			}
			return ",directInputMs=" + millis(directInputNs) +
					",directScatterMs=" + millis(directScatterNs) +
					",directUpdateMs=" + millis(directUpdateNs) +
					",directPivotMs=" + millis(directPivotNs) +
					",directScaleMs=" + millis(directScaleNs) +
					",directClearMs=" + millis(directClearNs) +
					",directFinalizeMs=" + millis(directFinalizeNs) +
					",directPriorUpdateTargets=" + directPriorUpdateTargets +
					",directFutureUpdateTargets=" + directFutureUpdateTargets +
						",directPriorNewTouches=" + directPriorNewTouches +
						",directFutureNewTouches=" + directFutureNewTouches +
						",directActiveAdds=" + directActiveAdds +
						",directActivePolls=" + directActivePolls +
						",directRealUUpdates=" + directRealUUpdates +
						",directRealLUpdates=" + directRealLUpdates +
						",directPureRealUpdates=" + directPureRealUpdates +
						",directUpdateColumns=" + directUpdateColumns +
						",directTinyUpdateColumns=" + directTinyUpdateColumns +
						",directSmallUpdateColumns=" + directSmallUpdateColumns +
						",directMediumUpdateColumns=" + directMediumUpdateColumns +
						",directLargeUpdateColumns=" + directLargeUpdateColumns +
						",directTinyUpdateEntries=" + directTinyUpdateEntries +
						",directSmallUpdateEntries=" + directSmallUpdateEntries +
						",directMediumUpdateEntries=" + directMediumUpdateEntries +
						",directLargeUpdateEntries=" + directLargeUpdateEntries;
			}
	}

	private static final class Complex
	{
		final double re;
		final double im;

		Complex(double re, double im)
		{
			this.re = re;
			this.im = im;
		}

		double abs()
		{
			return Zklu_complex.abs(re, im);
		}

		boolean isZero()
		{
			return re == 0.0 && im == 0.0;
		}

		Complex divide(Complex other)
		{
			double scalar;
			if (Math.abs(other.re) >= Math.abs(other.im))
			{
				double ratio = other.im / other.re;
				scalar = 1.0 / (other.re + other.im * ratio);
				return new Complex(scalar * (re + im * ratio),
						scalar * (im - re * ratio));
			}
			double ratio = other.re / other.im;
			scalar = 1.0 / (other.re * ratio + other.im);
			return new Complex(scalar * (re * ratio + im),
					scalar * (im * ratio - re));
		}
	}
}
