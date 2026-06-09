package edu.ufl.cise.klu.bench;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.StringTokenizer;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdcomplex.Zklu_factor.klu_z_factor;
import static edu.ufl.cise.klu.tdcomplex.Zklu_refactor.klu_z_refactor;
import static edu.ufl.cise.klu.tdcomplex.Zklu_solve.klu_z_solve;

/**
 * Matrix Market replay benchmark for InterPSS/KLUSolve-exported matrices.
 *
 * Matrix file format: Matrix Market coordinate, complex/real/integer/pattern.
 * RHS file format: Matrix Market array, complex/real/integer.
 */
public class ZkluMatrixMarketBenchmark
{
	private static final int DEFAULT_WARMUPS = 1;
	private static final int DEFAULT_ITERATIONS = 3;

	public static void main(String[] args) throws Exception
	{
		if (args.length < 1)
		{
			System.out.println("usage: ZkluMatrixMarketBenchmark <matrix.mtx> [rhs.mtx] [warmups] [iterations]");
			return;
		}
		String matrixFile = args[0];
		String rhsFile = args.length > 1 ? args[1] : null;
		int warmups = intArg(args, 2, DEFAULT_WARMUPS);
		int iterations = intArg(args, 3, DEFAULT_ITERATIONS);

		Matrix matrix = readMatrix(matrixFile);
		double[] rhs = rhsFile == null || "-".equals(rhsFile) ?
				onesRhs(matrix) : readVector(rhsFile, matrix.n);
		Result result = run(matrix, rhs, warmups, iterations);
		System.out.println("case=matrix-market,file=" + matrixFile +
				",n=" + matrix.n + ",nnz=" + matrix.Ap[matrix.n] +
				",iterations=" + iterations +
				",analyzeMs=" + result.analyzeMs +
				",factorMs=" + result.factorMs +
				",refactorMs=" + result.refactorMs +
				",solveMs=" + result.solveMs +
					",luEntries=" + result.luEntries +
					",fillRatio=" + result.fillRatio +
					",maxResidual=" + result.maxResidual +
					",relativeResidual=" + result.relativeResidual);
	}

	public static Result run(Matrix matrix, double[] rhs, int warmups, int iterations)
	{
		if (Boolean.getBoolean("jklu.benchmark.reuseSymbolic"))
		{
			return runReuseSymbolic(matrix, rhs, warmups, iterations);
		}
		for (int i = 0; i < warmups; i++)
		{
			runOnce(matrix, rhs);
		}

		long analyzeNs = 0;
		long factorNs = 0;
		long refactorNs = 0;
		long solveNs = 0;
		double residual = 0.0;
		int luEntries = 0;
		for (int i = 0; i < iterations; i++)
		{
			Iteration iteration = runOnce(matrix, rhs);
			analyzeNs += iteration.analyzeNs;
			factorNs += iteration.factorNs;
			refactorNs += iteration.refactorNs;
			solveNs += iteration.solveNs;
			residual = Math.max(residual, iteration.residual);
			luEntries = iteration.luEntries;
		}

		Result result = new Result();
		result.analyzeMs = millis(analyzeNs, iterations);
		result.factorMs = millis(factorNs, iterations);
		result.refactorMs = millis(refactorNs, iterations);
		result.solveMs = millis(solveNs, iterations);
			result.maxResidual = residual;
			result.relativeResidual = residual / Math.max(rhsNorm(rhs), 1.0);
			result.luEntries = luEntries;
			result.fillRatio = (double) luEntries / matrix.Ap[matrix.n];
			return result;
		}

	private static Result runReuseSymbolic(Matrix matrix, double[] rhs,
			int warmups, int iterations)
	{
			KLU_common common = new KLU_common();
			klu_defaults(common);
			common.scale = 0;
			common.ordering = Integer.getInteger("jklu.benchmark.ordering",
					common.ordering).intValue();
			common.btf = Integer.getInteger("jklu.benchmark.btf",
					common.btf).intValue();

		long t0 = System.nanoTime();
		KLU_symbolic symbolic = klu_analyze(matrix.n, matrix.Ap, matrix.Ai, common);
		long t1 = System.nanoTime();
		if (symbolic == null)
		{
			throw new IllegalStateException("analyze failed, status=" + common.status);
		}

		for (int i = 0; i < warmups; i++)
		{
			runNumericOnceTimed(matrix, rhs, symbolic, common);
		}

		long factorNs = 0;
		long refactorNs = 0;
		long solveNs = 0;
		double residual = 0.0;
		int luEntries = 0;
		for (int i = 0; i < iterations; i++)
		{
			Iteration iteration = runNumericOnceTimed(matrix, rhs, symbolic, common);
			factorNs += iteration.factorNs;
			refactorNs += iteration.refactorNs;
			solveNs += iteration.solveNs;
			residual = Math.max(residual, iteration.residual);
			luEntries = iteration.luEntries;
		}

		Result result = new Result();
		result.analyzeMs = millis(t1 - t0, 1);
		result.factorMs = millis(factorNs, iterations);
		result.refactorMs = millis(refactorNs, iterations);
		result.solveMs = millis(solveNs, iterations);
		result.maxResidual = residual;
		result.relativeResidual = residual / Math.max(rhsNorm(rhs), 1.0);
		result.luEntries = luEntries;
		result.fillRatio = (double) luEntries / matrix.Ap[matrix.n];
		return result;
	}

	private static Iteration runNumericOnceTimed(Matrix matrix, double[] rhs,
			KLU_symbolic symbolic, KLU_common common)
	{
		long t0 = System.nanoTime();
		KLU_z_numeric numeric = klu_z_factor(matrix.Ap, matrix.Ai, matrix.Ax, symbolic, common);
		long t1 = System.nanoTime();
		if (numeric == null)
		{
			throw new IllegalStateException("factorization failed, status=" + common.status);
		}

		if (klu_z_refactor(matrix.Ap, matrix.Ai, matrix.Ax, symbolic, numeric, common) == 0)
		{
			throw new IllegalStateException("refactorization failed, status=" + common.status);
		}
		long t2 = System.nanoTime();

		double[] b = rhs.clone();
		if (klu_z_solve(symbolic, numeric, matrix.n, 1, b, common) == 0)
		{
			throw new IllegalStateException("solve failed, status=" + common.status);
		}
		long t3 = System.nanoTime();

		Iteration result = new Iteration();
		result.factorNs = t1 - t0;
		result.refactorNs = t2 - t1;
		result.solveNs = t3 - t2;
		result.residual = residual(matrix, b, rhs);
		result.luEntries = numeric.lnz + numeric.unz - numeric.n;
		return result;
	}

	private static Iteration runOnce(Matrix matrix, double[] rhs)
	{
			KLU_common common = new KLU_common();
			klu_defaults(common);
			common.scale = 0;
			common.ordering = Integer.getInteger("jklu.benchmark.ordering",
					common.ordering).intValue();
			common.btf = Integer.getInteger("jklu.benchmark.btf",
					common.btf).intValue();

		long t0 = System.nanoTime();
		KLU_symbolic symbolic = klu_analyze(matrix.n, matrix.Ap, matrix.Ai, common);
		long t1 = System.nanoTime();
		KLU_z_numeric numeric = klu_z_factor(matrix.Ap, matrix.Ai, matrix.Ax, symbolic, common);
		long t2 = System.nanoTime();
		if (numeric == null)
		{
			throw new IllegalStateException("factorization failed, status=" + common.status);
		}

		if (klu_z_refactor(matrix.Ap, matrix.Ai, matrix.Ax, symbolic, numeric, common) == 0)
		{
			throw new IllegalStateException("refactorization failed, status=" + common.status);
		}
		long t3 = System.nanoTime();

		double[] b = rhs.clone();
		if (klu_z_solve(symbolic, numeric, matrix.n, 1, b, common) == 0)
		{
			throw new IllegalStateException("solve failed, status=" + common.status);
		}
		long t4 = System.nanoTime();

		Iteration result = new Iteration();
		result.analyzeNs = t1 - t0;
		result.factorNs = t2 - t1;
		result.refactorNs = t3 - t2;
		result.solveNs = t4 - t3;
		result.residual = residual(matrix, b, rhs);
		result.luEntries = numeric.lnz + numeric.unz - numeric.n;
		return result;
	}

	public static Matrix readMatrix(String file) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		try
		{
			String header = reader.readLine();
			if (header == null)
			{
				throw new IOException("unsupported Matrix Market matrix header: " + header);
			}
			String[] headerTokens = header.toLowerCase().split("\\s+");
			if (headerTokens.length < 5 ||
				!"%%matrixmarket".equals(headerTokens[0]) ||
				!"matrix".equals(headerTokens[1]) ||
				!"coordinate".equals(headerTokens[2]))
			{
				throw new IOException("unsupported Matrix Market matrix header: " + header);
			}
			String field = headerTokens[3];
			String symmetry = headerTokens[4];
			String line = nextDataLine(reader);
			StringTokenizer size = new StringTokenizer(line);
			int rows = Integer.parseInt(size.nextToken());
			int cols = Integer.parseInt(size.nextToken());
			int nz = Integer.parseInt(size.nextToken());
			if (rows != cols)
			{
				throw new IOException("JKLU benchmark requires a square matrix: " + rows + "x" + cols);
			}
			List<Entry> entries = new ArrayList<Entry>(nz);
			for (int i = 0; i < nz; i++)
			{
				line = nextDataLine(reader);
				StringTokenizer tokens = new StringTokenizer(line);
				Entry entry = new Entry();
				entry.row = Integer.parseInt(tokens.nextToken()) - 1;
				entry.col = Integer.parseInt(tokens.nextToken()) - 1;
				readMatrixValue(tokens, field, entry);
				entries.add(entry);
				addSymmetricEntry(entries, entry, symmetry);
			}
			return compress(rows, entries);
		}
		finally
		{
			reader.close();
		}
	}

	public static double[] readVector(String file, int n) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		try
		{
			String header = reader.readLine();
			if (header == null)
			{
				throw new IOException("unsupported Matrix Market vector header: " + header);
			}
			String[] headerTokens = header.toLowerCase().split("\\s+");
			if (headerTokens.length < 5 ||
				!"%%matrixmarket".equals(headerTokens[0]) ||
				!"matrix".equals(headerTokens[1]) ||
				!"array".equals(headerTokens[2]))
			{
				throw new IOException("unsupported Matrix Market vector header: " + header);
			}
			String field = headerTokens[3];
			String line = nextDataLine(reader);
			StringTokenizer size = new StringTokenizer(line);
			int rows = Integer.parseInt(size.nextToken());
			int cols = Integer.parseInt(size.nextToken());
			if (rows != n || cols != 1)
			{
				throw new IOException("RHS vector size mismatch: " + rows + "x" + cols + ", n=" + n);
			}
			double[] values = new double[2 * n];
			for (int row = 0; row < n; row++)
			{
				line = nextDataLine(reader);
				StringTokenizer tokens = new StringTokenizer(line);
				values[2 * row] = Double.parseDouble(tokens.nextToken());
				values[2 * row + 1] = "complex".equals(field) ?
						Double.parseDouble(tokens.nextToken()) : 0.0;
			}
			return values;
		}
		finally
		{
			reader.close();
		}
	}

	private static void readMatrixValue(StringTokenizer tokens, String field, Entry entry)
			throws IOException
	{
		if ("pattern".equals(field))
		{
			entry.re = 1.0;
			entry.im = 0.0;
		}
		else if ("complex".equals(field))
		{
			entry.re = Double.parseDouble(tokens.nextToken());
			entry.im = Double.parseDouble(tokens.nextToken());
		}
		else if ("real".equals(field) || "integer".equals(field))
		{
			entry.re = Double.parseDouble(tokens.nextToken());
			entry.im = 0.0;
		}
		else
		{
			throw new IOException("unsupported Matrix Market matrix field: " + field);
		}
	}

	private static void addSymmetricEntry(List<Entry> entries, Entry entry, String symmetry)
			throws IOException
	{
		if ("general".equals(symmetry) || entry.row == entry.col)
		{
			return;
		}
		Entry mirror = new Entry();
		mirror.row = entry.col;
		mirror.col = entry.row;
		if ("symmetric".equals(symmetry))
		{
			mirror.re = entry.re;
			mirror.im = entry.im;
		}
		else if ("skew-symmetric".equals(symmetry))
		{
			mirror.re = -entry.re;
			mirror.im = -entry.im;
		}
		else if ("hermitian".equals(symmetry))
		{
			mirror.re = entry.re;
			mirror.im = -entry.im;
		}
		else
		{
			throw new IOException("unsupported Matrix Market matrix symmetry: " + symmetry);
		}
		entries.add(mirror);
	}

	private static Matrix compress(int n, List<Entry> entries)
	{
		Collections.sort(entries, new Comparator<Entry>()
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
		int[] Ap = new int[n + 1];
		int[] Ai = new int[entries.size()];
		double[] Ax = new double[2 * entries.size()];
		int p = 0;
		for (int col = 0; col < n; col++)
		{
			Ap[col] = p;
			while (p < entries.size() && entries.get(p).col == col)
			{
				Entry entry = entries.get(p);
				Ai[p] = entry.row;
				Ax[2 * p] = entry.re;
				Ax[2 * p + 1] = entry.im;
				p++;
			}
		}
		Ap[n] = p;

		Matrix matrix = new Matrix();
		matrix.n = n;
		matrix.Ap = Ap;
		matrix.Ai = Ai;
		matrix.Ax = Ax;
		return matrix;
	}

	private static double[] onesRhs(Matrix matrix)
	{
		double[] rhs = new double[2 * matrix.n];
		for (int col = 0; col < matrix.n; col++)
		{
			for (int p = matrix.Ap[col]; p < matrix.Ap[col + 1]; p++)
			{
				int row = matrix.Ai[p];
				rhs[2 * row] += matrix.Ax[2 * p];
				rhs[2 * row + 1] += matrix.Ax[2 * p + 1];
			}
		}
		return rhs;
	}

		private static double residual(Matrix matrix, double[] x, double[] b)
	{
		double[] ax = new double[2 * matrix.n];
		for (int col = 0; col < matrix.n; col++)
		{
			double xr = x[2 * col];
			double xi = x[2 * col + 1];
			for (int p = matrix.Ap[col]; p < matrix.Ap[col + 1]; p++)
			{
				int row = matrix.Ai[p];
				double ar = matrix.Ax[2 * p];
				double ai = matrix.Ax[2 * p + 1];
				ax[2 * row] += ar * xr - ai * xi;
				ax[2 * row + 1] += ai * xr + ar * xi;
			}
		}
		double max = 0.0;
		for (int row = 0; row < matrix.n; row++)
		{
			double dr = ax[2 * row] - b[2 * row];
			double di = ax[2 * row + 1] - b[2 * row + 1];
			max = Math.max(max, Math.sqrt(dr * dr + di * di));
		}
			return max;
		}

		private static double rhsNorm(double[] b)
		{
			double max = 0.0;
			for (int row = 0; row < b.length / 2; row++)
			{
				double re = b[2 * row];
				double im = b[2 * row + 1];
				max = Math.max(max, Math.sqrt(re * re + im * im));
			}
			return max;
		}

	private static String nextDataLine(BufferedReader reader) throws IOException
	{
		String line;
		while ((line = reader.readLine()) != null)
		{
			line = line.trim();
			if (line.length() != 0 && !line.startsWith("%"))
			{
				return line;
			}
		}
		throw new IOException("unexpected end of Matrix Market file");
	}

	private static int intArg(String[] args, int index, int defaultValue)
	{
		if (args.length <= index)
		{
			return defaultValue;
		}
		return Integer.parseInt(args[index]);
	}

	private static double millis(long nanos, int iterations)
	{
		return nanos / 1.0e6 / iterations;
	}

	public static final class Matrix
	{
		public int n;
		public int[] Ap;
		public int[] Ai;
		public double[] Ax;
	}

		public static final class Result
		{
			public double analyzeMs;
			public double factorMs;
			public double refactorMs;
			public double solveMs;
			public double maxResidual;
			public double relativeResidual;
			public double fillRatio;
			public int luEntries;
		}

	private static final class Iteration
	{
		long analyzeNs;
		long factorNs;
		long refactorNs;
		long solveNs;
		double residual;
		int luEntries;
	}

	private static final class Entry
	{
		int row;
		int col;
		double re;
		double im;
	}
}
