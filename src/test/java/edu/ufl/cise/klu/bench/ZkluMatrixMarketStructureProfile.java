package edu.ufl.cise.klu.bench;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import edu.ufl.cise.klu.bench.ZkluMatrixMarketBenchmark.Matrix;
import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;

import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;

/**
 * Structural profiler for large Matrix Market replay cases.
 */
public class ZkluMatrixMarketStructureProfile
{
	public static void main(String[] args) throws Exception
	{
		if (args.length < 1)
		{
			System.out.println("usage: ZkluMatrixMarketStructureProfile <matrix.mtx>");
			return;
		}

		long t0 = System.nanoTime();
		Matrix matrix = ZkluMatrixMarketBenchmark.readMatrix(args[0]);
		long t1 = System.nanoTime();

		Stats stats = matrixStats(matrix);
		Components components = components(matrix);

		KLU_common common = new KLU_common();
		klu_defaults(common);
		long t2 = System.nanoTime();
		KLU_symbolic symbolic = klu_analyze(matrix.n, matrix.Ap, matrix.Ai, common);
		long t3 = System.nanoTime();

		System.out.println("file=" + args[0]);
		System.out.println("readMs=" + millis(t1 - t0));
		System.out.println("n=" + matrix.n + ",nnz=" + matrix.Ap[matrix.n] +
				",diagNnz=" + stats.diagNnz +
				",rowNnzMax=" + stats.rowMax +
				",colNnzMax=" + stats.colMax +
				",rowNnzAvg=" + stats.rowAvg +
				",colNnzAvg=" + stats.colAvg);
		System.out.println("undirectedIslands=" + components.count +
				",largestIsland=" + components.largest +
				",singletonIslands=" + components.singletons +
				",topIslandSizes=" + components.topSizes);
		if (symbolic == null)
		{
			System.out.println("kluAnalyzeStatus=" + common.status);
			return;
		}
		System.out.println("analyzeMs=" + millis(t3 - t2) +
				",status=" + common.status +
				",btf=" + symbolic.do_btf +
				",ordering=" + symbolic.ordering +
				",structuralRank=" + symbolic.structural_rank);
		System.out.println("btfBlocks=" + symbolic.nblocks +
				",maxBlock=" + symbolic.maxblock +
				",nzoff=" + symbolic.nzoff +
				",symmetry=" + symbolic.symmetry +
				",estLnz=" + symbolic.lnz +
				",estUnz=" + symbolic.unz +
				",estFlops=" + symbolic.est_flops);
		System.out.println("symbolicPermutation pIdentity=" + identity(symbolic.P) +
				",qIdentity=" + identity(symbolic.Q) +
				",pEqualsQ=" + same(symbolic.P, symbolic.Q) +
				",naturalBlockCompatible=" + naturalBlockCompatible(matrix, symbolic));
		System.out.println("btfBlockSizeSummary=" + blockSummary(symbolic));
		System.out.println("currentDirectKernelEffectiveBlocks=" + symbolic.nblocks +
				",currentDirectKernelEffectiveMaxBlock=" + symbolic.maxblock +
				",note=tdcomplex factorization applies Symbolic.P/Q/R boundaries but still lacks full SuiteSparse-style per-block numeric storage");
	}

	private static Stats matrixStats(Matrix matrix)
	{
		int n = matrix.n;
		int[] rowCount = new int[n];
		int diag = 0;
		int colMax = 0;
		for (int col = 0; col < n; col++)
		{
			int colNnz = matrix.Ap[col + 1] - matrix.Ap[col];
			colMax = Math.max(colMax, colNnz);
			for (int p = matrix.Ap[col]; p < matrix.Ap[col + 1]; p++)
			{
				int row = matrix.Ai[p];
				rowCount[row]++;
				if (row == col)
				{
					diag++;
				}
			}
		}
		int rowMax = 0;
		for (int count : rowCount)
		{
			rowMax = Math.max(rowMax, count);
		}
		Stats stats = new Stats();
		stats.diagNnz = diag;
		stats.rowMax = rowMax;
		stats.colMax = colMax;
		stats.rowAvg = ((double) matrix.Ap[n]) / n;
		stats.colAvg = stats.rowAvg;
		return stats;
	}

	private static Components components(Matrix matrix)
	{
		int n = matrix.n;
		int[] degree = new int[n];
		for (int col = 0; col < n; col++)
		{
			for (int p = matrix.Ap[col]; p < matrix.Ap[col + 1]; p++)
			{
				int row = matrix.Ai[p];
				if (row != col)
				{
					degree[row]++;
					degree[col]++;
				}
			}
		}
		int[][] graph = new int[n][];
		for (int i = 0; i < n; i++)
		{
			graph[i] = new int[degree[i]];
		}
		int[] next = new int[n];
		for (int col = 0; col < n; col++)
		{
			for (int p = matrix.Ap[col]; p < matrix.Ap[col + 1]; p++)
			{
				int row = matrix.Ai[p];
				if (row != col)
				{
					graph[row][next[row]++] = col;
					graph[col][next[col]++] = row;
				}
			}
		}

		boolean[] seen = new boolean[n];
		ArrayDeque<Integer> queue = new ArrayDeque<Integer>();
		List<Integer> sizes = new ArrayList<Integer>();
		for (int seed = 0; seed < n; seed++)
		{
			if (seen[seed])
			{
				continue;
			}
			seen[seed] = true;
			queue.add(Integer.valueOf(seed));
			int size = 0;
			while (!queue.isEmpty())
			{
				int node = queue.removeFirst().intValue();
				size++;
				for (int nbr : graph[node])
				{
					if (!seen[nbr])
					{
						seen[nbr] = true;
						queue.add(Integer.valueOf(nbr));
					}
				}
			}
			sizes.add(Integer.valueOf(size));
		}
		Collections.sort(sizes, Collections.reverseOrder());
		Components components = new Components();
		components.count = sizes.size();
		components.largest = sizes.isEmpty() ? 0 : sizes.get(0).intValue();
		components.singletons = 0;
		for (Integer size : sizes)
		{
			if (size.intValue() == 1)
			{
				components.singletons++;
			}
		}
		components.topSizes = sizes.subList(0, Math.min(10, sizes.size())).toString();
		return components;
	}

	private static String blockSummary(KLU_symbolic symbolic)
	{
		List<Integer> sizes = new ArrayList<Integer>();
		int tiny = 0;
		for (int block = 0; block < symbolic.nblocks; block++)
		{
			int size = symbolic.R[block + 1] - symbolic.R[block];
			sizes.add(Integer.valueOf(size));
			if (size <= 3)
			{
				tiny++;
			}
		}
		Collections.sort(sizes, Collections.reverseOrder());
		int topCount = Math.min(15, sizes.size());
		List<Integer> top = sizes.subList(0, topCount);
		double avg = 0.0;
		for (Integer size : sizes)
		{
			avg += size.intValue();
		}
		avg = sizes.isEmpty() ? 0.0 : avg / sizes.size();
		return "avg=" + avg + ",tinyBlocksLe3=" + tiny + ",topSizes=" + top;
	}

	private static boolean identity(int[] permutation)
	{
		for (int i = 0; i < permutation.length; i++)
		{
			if (permutation[i] != i)
			{
				return false;
			}
		}
		return true;
	}

	private static boolean same(int[] a, int[] b)
	{
		if (a.length != b.length)
		{
			return false;
		}
		for (int i = 0; i < a.length; i++)
		{
			if (a[i] != b[i])
			{
				return false;
			}
		}
		return true;
	}

	private static boolean naturalBlockCompatible(Matrix matrix, KLU_symbolic symbolic)
	{
		int[] blockOf = new int[matrix.n];
		for (int block = 0; block < symbolic.nblocks; block++)
		{
			for (int k = symbolic.R[block]; k < symbolic.R[block + 1]; k++)
			{
				blockOf[k] = block;
			}
		}
		for (int col = 0; col < matrix.n; col++)
		{
			for (int p = matrix.Ap[col]; p < matrix.Ap[col + 1]; p++)
			{
				if (blockOf[matrix.Ai[p]] != blockOf[col])
				{
					return false;
				}
			}
		}
		return true;
	}

	private static double millis(long nanos)
	{
		return nanos / 1.0e6;
	}

	private static final class Stats
	{
		int diagNnz;
		int rowMax;
		int colMax;
		double rowAvg;
		double colAvg;
	}

	private static final class Components
	{
		int count;
		int largest;
		int singletons;
		String topSizes;
	}
}
