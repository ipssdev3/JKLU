package edu.ufl.cise.klu.solver;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdouble.Dklu_version.KLU_OK;
import static edu.ufl.cise.klu.tdouble.Dklu_version.KLU_SINGULAR;
import static edu.ufl.cise.klu.tdcomplex.Zklu_factor.klu_z_factor;
import static edu.ufl.cise.klu.tdcomplex.Zklu_refactor.klu_z_refactor;
import static edu.ufl.cise.klu.tdcomplex.Zklu_solve.klu_z_solve;

/**
 * High-level complex sparse system facade inspired by KLUSolveX.
 *
 * Matrix element APIs use 1-based node numbers. Node 0 is treated as ground and
 * skipped when stamping primitive matrices.
 */
public class KluComplexSparseSet
{
	public static final int OK = 1;
	public static final int SINGULAR = 2;
	public static final int ERROR = 0;

	private final int size;
	private final Map<Long, ComplexValue> entries;
	private KLU_common common;
	private KLU_symbolic symbolic;
	private KLU_z_numeric numeric;
	private int[] Ap;
	private int[] Ai;
	private double[] Ax;
	private Map<Long, Integer> cscPosition;
	private boolean compressedDirty;
	private boolean factored;
	private int sparseNNZ;
	private int singularCol;

	public KluComplexSparseSet(int size)
	{
		if (size < 0)
		{
			throw new IllegalArgumentException("size must be non-negative");
		}
		this.size = size;
		this.entries = new HashMap<Long, ComplexValue>();
		this.common = new KLU_common();
		klu_defaults(this.common);
		this.common.halt_if_singular = 0;
		this.compressedDirty = true;
	}

	public int getSize()
	{
		return size;
	}

	public int getNNZ()
	{
		if (Ap != null && !compressedDirty)
		{
			return Ap[size];
		}
		int count = 0;
		for (ComplexValue value : entries.values())
		{
			if (!value.isZero())
			{
				count++;
			}
		}
		return count;
	}

	public int getSparseNNZ()
	{
		return sparseNNZ;
	}

	public int getSingularCol()
	{
		return singularCol;
	}

	public KLU_common getCommon()
	{
		return common;
	}

	public void zero()
	{
		entries.clear();
		Ap = null;
		Ai = null;
		Ax = null;
		cscPosition = null;
		symbolic = null;
		numeric = null;
		sparseNNZ = 0;
		singularCol = 0;
		compressedDirty = true;
		factored = false;
		common = new KLU_common();
		klu_defaults(common);
		common.halt_if_singular = 0;
	}

	public boolean addMatrixElement(int row, int col, double re, double im)
	{
		if (!validNode(row) || !validNode(col))
		{
			return false;
		}
		if (row == 0 || col == 0 || (re == 0.0 && im == 0.0))
		{
			return true;
		}
		addDirected(row, col, re, im);
		if (row != col)
		{
			addDirected(col, row, re, im);
		}
		return true;
	}

	public boolean setMatrixElement(int row, int col, double re, double im)
	{
		if (!validNode(row) || !validNode(col))
		{
			return false;
		}
		if (row == 0 || col == 0)
		{
			return true;
		}
		long key = key(row - 1, col - 1);
		ComplexValue value = entries.get(key);
		if (value == null)
		{
			if (re == 0.0 && im == 0.0)
			{
				return true;
			}
			entries.put(key, new ComplexValue(re, im));
			compressedDirty = true;
			factored = false;
			return true;
		}
		value.re = re;
		value.im = im;
		updateCompressedValue(key, re, im);
		return true;
	}

	public double[] getMatrixElement(int row, int col)
	{
		if (!validNode(row) || !validNode(col) || row == 0 || col == 0)
		{
			return new double[] {0.0, 0.0};
		}
		ComplexValue value = entries.get(key(row - 1, col - 1));
		if (value == null)
		{
			return new double[] {0.0, 0.0};
		}
		return new double[] {value.re, value.im};
	}

	public boolean addPrimitiveMatrix(int[] nodes, double[] yMatrix)
	{
		if (nodes == null || yMatrix == null || yMatrix.length < 2 * nodes.length * nodes.length)
		{
			return false;
		}
		int order = nodes.length;
		for (int i = 0; i < order; i++)
		{
			if (!validNode(nodes[i]))
			{
				return false;
			}
		}
		for (int i = 0; i < order; i++)
		{
			if (nodes[i] == 0)
			{
				continue;
			}
			for (int j = 0; j < order; j++)
			{
				int valueIndex = 2 * (i * order + j);
				if (nodes[j] == 0)
				{
					continue;
				}
				double re = yMatrix[valueIndex];
				double im = yMatrix[valueIndex + 1];
				if (re != 0.0 || im != 0.0)
				{
					addDirected(nodes[i], nodes[j], re, im);
				}
			}
		}
		return true;
	}

	public int factor()
	{
		if (compressedDirty || Ap == null)
		{
			compress();
			symbolic = klu_analyze(size, Ap, Ai, common);
			numeric = symbolic == null ? null : klu_z_factor(Ap, Ai, Ax, symbolic, common);
		}
		else if (numeric != null)
		{
			return refactor();
		}
		else
		{
			symbolic = klu_analyze(size, Ap, Ai, common);
			numeric = symbolic == null ? null : klu_z_factor(Ap, Ai, Ax, symbolic, common);
		}
		factored = numeric != null && common.status == KLU_OK;
		updateFactorStats();
		if (factored)
		{
			return OK;
		}
		if (common.status == KLU_SINGULAR)
		{
			return SINGULAR;
		}
		return ERROR;
	}

	public int refactor()
	{
		if (compressedDirty || Ap == null || symbolic == null || numeric == null)
		{
			return factor();
		}
		int ok = klu_z_refactor(Ap, Ai, Ax, symbolic, numeric, common);
		factored = ok == OK && common.status == KLU_OK;
		if (!factored)
		{
			numeric = null;
		}
		updateFactorStats();
		if (factored)
		{
			return OK;
		}
		if (common.status == KLU_SINGULAR)
		{
			return SINGULAR;
		}
		return ERROR;
	}

	public int solve(double[] rhs, double[] solution)
	{
		if (rhs == null || solution == null || rhs.length < 2 * size || solution.length < 2 * size)
		{
			return ERROR;
		}
		System.arraycopy(rhs, 0, solution, 0, 2 * size);
		return solveInPlace(solution);
	}

	public int solveInPlace(double[] rhsAndSolution)
	{
		if (rhsAndSolution == null || rhsAndSolution.length < 2 * size)
		{
			return ERROR;
		}
		if (!factored)
		{
			int factorStatus = factor();
			if (factorStatus != OK)
			{
				return factorStatus;
			}
		}
		return klu_z_solve(symbolic, numeric, size, 1, rhsAndSolution, common) == OK ? OK : ERROR;
	}

	public boolean incrementMatrixElement(int row, int col, double re, double im)
	{
		if (!validExistingCompressedEntry(row, col))
		{
			return false;
		}
		long key = key(row - 1, col - 1);
		ComplexValue value = entries.get(key);
		value.re += re;
		value.im += im;
		updateCompressedValue(key, value.re, value.im);
		return true;
	}

	public boolean zeroiseMatrixElement(int row, int col)
	{
		if (!validExistingCompressedEntry(row, col))
		{
			return false;
		}
		long key = key(row - 1, col - 1);
		ComplexValue value = entries.get(key);
		value.re = 0.0;
		value.im = 0.0;
		updateCompressedValue(key, 0.0, 0.0);
		return true;
	}

	public CscMatrix getCompressedMatrix()
	{
		if (compressedDirty || Ap == null)
		{
			compress();
		}
		return new CscMatrix(size, copy(Ap), copy(Ai), copy(Ax));
	}

	public TripletMatrix getTripletMatrix()
	{
		if (compressedDirty || Ap == null)
		{
			compress();
		}
		int nnz = Ap[size];
		int[] rows = new int[nnz];
		int[] cols = new int[nnz];
		double[] values = new double[2 * nnz];
		int p2 = 0;
		for (int col = 0; col < size; col++)
		{
			for (int p = Ap[col]; p < Ap[col + 1]; p++)
			{
				rows[p2] = Ai[p];
				cols[p2] = col;
				values[2 * p2] = Ax[2 * p];
				values[2 * p2 + 1] = Ax[2 * p + 1];
				p2++;
			}
		}
		return new TripletMatrix(size, rows, cols, values);
	}

	public int[] findIslands()
	{
		if (compressedDirty || Ap == null)
		{
			compress();
		}
		List<List<Integer>> adjacency = new ArrayList<List<Integer>>(size);
		for (int i = 0; i < size; i++)
		{
			adjacency.add(new ArrayList<Integer>());
		}
		for (int col = 0; col < size; col++)
		{
			for (int p = Ap[col]; p < Ap[col + 1]; p++)
			{
				int row = Ai[p];
				adjacency.get(col).add(row);
				if (row != col)
				{
					adjacency.get(row).add(col);
				}
			}
		}
		int[] islands = new int[size];
		int island = 0;
		List<Integer> stack = new ArrayList<Integer>();
		for (int i = 0; i < size; i++)
		{
			if (islands[i] != 0)
			{
				continue;
			}
			island++;
			stack.add(i);
			islands[i] = island;
			while (!stack.isEmpty())
			{
				int node = stack.remove(stack.size() - 1);
				for (Integer next : adjacency.get(node))
				{
					if (islands[next.intValue()] == 0)
					{
						islands[next.intValue()] = island;
						stack.add(next);
					}
				}
			}
		}
		return islands;
	}

	public int findDisconnectedSubnetwork()
	{
		if (!factored)
		{
			factor();
		}
		return singularCol;
	}

	public boolean saveAsMatrixMarket(String matrixFile, double[] b, String vectorFile)
			throws IOException
	{
		CscMatrix matrix = getCompressedMatrix();
		PrintWriter matrixOut = new PrintWriter(new FileWriter(matrixFile));
		try
		{
			matrixOut.println("%%MatrixMarket matrix coordinate complex general");
			matrixOut.println(size + " " + size + " " + matrix.values.length / 2);
			for (int col = 0; col < size; col++)
			{
				for (int p = matrix.colPointers[col]; p < matrix.colPointers[col + 1]; p++)
				{
					matrixOut.println((matrix.rowIndices[p] + 1) + " " + (col + 1) + " " +
						matrix.values[2 * p] + " " + matrix.values[2 * p + 1]);
				}
			}
		}
		finally
		{
			matrixOut.close();
		}
		if (b != null && vectorFile != null)
		{
			PrintWriter vectorOut = new PrintWriter(new FileWriter(vectorFile));
			try
			{
				vectorOut.println("%%MatrixMarket matrix array complex general");
				vectorOut.println(size + " 1");
				for (int i = 0; i < size; i++)
				{
					vectorOut.println(b[2 * i] + " " + b[2 * i + 1]);
				}
			}
			finally
			{
				vectorOut.close();
			}
		}
		return true;
	}

	private boolean validNode(int node)
	{
		return node >= 0 && node <= size;
	}

	private void addDirected(int row, int col, double re, double im)
	{
		long key = key(row - 1, col - 1);
		ComplexValue value = entries.get(key);
		if (value == null)
		{
			entries.put(key, new ComplexValue(re, im));
			compressedDirty = true;
			factored = false;
			return;
		}
		value.re += re;
		value.im += im;
		updateCompressedValue(key, value.re, value.im);
	}

	private boolean validExistingCompressedEntry(int row, int col)
	{
		if (!validNode(row) || !validNode(col) || row == 0 || col == 0 ||
			compressedDirty || cscPosition == null)
		{
			return false;
		}
		return cscPosition.containsKey(key(row - 1, col - 1));
	}

	private void updateCompressedValue(long key, double re, double im)
	{
		Integer position = cscPosition == null ? null : cscPosition.get(key);
		if (position == null)
		{
			compressedDirty = true;
		}
		else
		{
			Ax[2 * position.intValue()] = re;
			Ax[2 * position.intValue() + 1] = im;
		}
		factored = false;
	}

	private void compress()
	{
		List<Entry> nonzeros = new ArrayList<Entry>();
		for (Map.Entry<Long, ComplexValue> entry : entries.entrySet())
		{
			ComplexValue value = entry.getValue();
			if (!value.isZero())
			{
				int row = row(entry.getKey());
				int col = col(entry.getKey());
				nonzeros.add(new Entry(row, col, value.re, value.im));
			}
		}
		Collections.sort(nonzeros, new Comparator<Entry>()
		{
			public int compare(Entry a, Entry b)
			{
				if (a.col != b.col)
				{
					return a.col < b.col ? -1 : 1;
				}
				if (a.row == b.row)
				{
					return 0;
				}
				return a.row < b.row ? -1 : 1;
			}
		});
		Ap = new int[size + 1];
		Ai = new int[nonzeros.size()];
		Ax = new double[2 * nonzeros.size()];
		cscPosition = new HashMap<Long, Integer>();
		int p = 0;
		for (int col = 0; col < size; col++)
		{
			Ap[col] = p;
			while (p < nonzeros.size() && nonzeros.get(p).col == col)
			{
				Entry entry = nonzeros.get(p);
				Ai[p] = entry.row;
				Ax[2 * p] = entry.re;
				Ax[2 * p + 1] = entry.im;
				cscPosition.put(key(entry.row, entry.col), Integer.valueOf(p));
				p++;
			}
		}
		Ap[size] = p;
		compressedDirty = false;
		factored = false;
		sparseNNZ = 0;
	}

	private void updateFactorStats()
	{
		if (numeric != null)
		{
			sparseNNZ = numeric.lnz + numeric.unz - numeric.n +
				(numeric.Offp != null ? numeric.Offp[numeric.n] : 0);
		}
		else
		{
			sparseNNZ = 0;
		}
		if (common.singular_col >= 0 && common.singular_col < 2 * size)
		{
			singularCol = (common.singular_col % size) + 1;
		}
		else
		{
			singularCol = 0;
		}
	}

	private long key(int row, int col)
	{
		return ((long) col) * size + row;
	}

	private int row(long key)
	{
		return (int) (key % size);
	}

	private int col(long key)
	{
		return (int) (key / size);
	}

	private static int[] copy(int[] source)
	{
		int[] copy = new int[source.length];
		System.arraycopy(source, 0, copy, 0, source.length);
		return copy;
	}

	private static double[] copy(double[] source)
	{
		double[] copy = new double[source.length];
		System.arraycopy(source, 0, copy, 0, source.length);
		return copy;
	}

	private static final class ComplexValue
	{
		double re;
		double im;

		ComplexValue(double re, double im)
		{
			this.re = re;
			this.im = im;
		}

		boolean isZero()
		{
			return re == 0.0 && im == 0.0;
		}
	}

	private static final class Entry
	{
		final int row;
		final int col;
		final double re;
		final double im;

		Entry(int row, int col, double re, double im)
		{
			this.row = row;
			this.col = col;
			this.re = re;
			this.im = im;
		}
	}

	public static final class CscMatrix
	{
		public final int n;
		public final int[] colPointers;
		public final int[] rowIndices;
		public final double[] values;

		CscMatrix(int n, int[] colPointers, int[] rowIndices, double[] values)
		{
			this.n = n;
			this.colPointers = colPointers;
			this.rowIndices = rowIndices;
			this.values = values;
		}
	}

	public static final class TripletMatrix
	{
		public final int n;
		public final int[] rows;
		public final int[] cols;
		public final double[] values;

		TripletMatrix(int n, int[] rows, int[] cols, double[] values)
		{
			this.n = n;
			this.rows = rows;
			this.cols = cols;
			this.values = values;
		}
	}
}
