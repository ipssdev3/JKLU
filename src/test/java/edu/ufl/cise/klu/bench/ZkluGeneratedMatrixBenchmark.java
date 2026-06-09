package edu.ufl.cise.klu.bench;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdcomplex.Zklu_factor.klu_z_factor;
import static edu.ufl.cise.klu.tdcomplex.Zklu_refactor.klu_z_refactor;
import static edu.ufl.cise.klu.tdcomplex.Zklu_solve.klu_z_solve;

/**
 * Lightweight generated-matrix benchmark for the direct complex kernel.
 *
 * This class intentionally lives in test sources so it does not become part of
 * the published JKLU API. Run it after test compilation with a classpath that
 * includes target/classes and target/test-classes.
 */
public class ZkluGeneratedMatrixBenchmark
{
	private static final int DEFAULT_N = 1000;
	private static final int DEFAULT_HALFBAND = 3;
	private static final int DEFAULT_NRHS = 1;
	private static final int DEFAULT_WARMUPS = 2;
	private static final int DEFAULT_ITERATIONS = 5;

	public static void main(String[] args)
	{
		int n = intArg(args, 0, DEFAULT_N);
		int halfBand = intArg(args, 1, DEFAULT_HALFBAND);
		int nrhs = intArg(args, 2, DEFAULT_NRHS);
		int warmups = intArg(args, 3, DEFAULT_WARMUPS);
		int iterations = intArg(args, 4, DEFAULT_ITERATIONS);

		Matrix matrix = bandedMatrix(n, halfBand);
		System.out.println("case=banded-complex,n=" + n + ",nnz=" + matrix.Ap[n] +
				",halfBand=" + halfBand + ",nrhs=" + nrhs);

		for (int i = 0; i < warmups; i++)
		{
			runOnce(matrix, nrhs, false);
		}

		long analyzeNs = 0;
		long factorNs = 0;
		long refactorNs = 0;
		long solveNs = 0;
		double residual = 0.0;
		int luEntries = 0;
		for (int i = 0; i < iterations; i++)
		{
			Result result = runOnce(matrix, nrhs, true);
			analyzeNs += result.analyzeNs;
			factorNs += result.factorNs;
			refactorNs += result.refactorNs;
			solveNs += result.solveNs;
			residual = Math.max(residual, result.residual);
			luEntries = result.luEntries;
		}

		System.out.println("iterations=" + iterations +
				",analyzeMs=" + millis(analyzeNs, iterations) +
				",factorMs=" + millis(factorNs, iterations) +
				",refactorMs=" + millis(refactorNs, iterations) +
				",solveMs=" + millis(solveNs, iterations) +
				",luEntries=" + luEntries +
				",fillRatio=" + ((double) luEntries / matrix.Ap[n]) +
				",maxResidual=" + residual);
	}

	private static Result runOnce(Matrix matrix, int nrhs, boolean measureRefactor)
	{
		KLU_common common = new KLU_common();
		klu_defaults(common);
		common.scale = 0;

		long t0 = System.nanoTime();
		KLU_symbolic symbolic = klu_analyze(matrix.n, matrix.Ap, matrix.Ai, common);
		long t1 = System.nanoTime();
		KLU_z_numeric numeric = klu_z_factor(matrix.Ap, matrix.Ai, matrix.Ax, symbolic, common);
		long t2 = System.nanoTime();
		if (numeric == null)
		{
			throw new IllegalStateException("factorization failed, status=" + common.status);
		}

		long t3 = t2;
		if (measureRefactor)
		{
			if (klu_z_refactor(matrix.Ap, matrix.Ai, matrix.Ax, symbolic, numeric, common) == 0)
			{
				throw new IllegalStateException("refactorization failed, status=" + common.status);
			}
			t3 = System.nanoTime();
		}

		double[] rhs = rhsForOnesSolution(matrix, nrhs);
		if (klu_z_solve(symbolic, numeric, matrix.n, nrhs, rhs, common) == 0)
		{
			throw new IllegalStateException("solve failed, status=" + common.status);
		}
		long t4 = System.nanoTime();

		Result result = new Result();
		result.analyzeNs = t1 - t0;
		result.factorNs = t2 - t1;
		result.refactorNs = measureRefactor ? t3 - t2 : 0;
		result.solveNs = t4 - t3;
		result.residual = maxSolutionError(rhs, matrix.n, nrhs);
		result.luEntries = numeric.lnz + numeric.unz - numeric.n;
		return result;
	}

	private static Matrix bandedMatrix(int n, int halfBand)
	{
		int nnz = 0;
		for (int col = 0; col < n; col++)
		{
			nnz += Math.min(n - 1, col + halfBand) - Math.max(0, col - halfBand) + 1;
		}
		int[] Ap = new int[n + 1];
		int[] Ai = new int[nnz];
		double[] Ax = new double[2 * nnz];
		int p = 0;
		for (int col = 0; col < n; col++)
		{
			Ap[col] = p;
			int first = Math.max(0, col - halfBand);
			int last = Math.min(n - 1, col + halfBand);
			for (int row = first; row <= last; row++)
			{
				Ai[p] = row;
				if (row == col)
				{
					Ax[2 * p] = 10.0 + halfBand;
					Ax[2 * p + 1] = 0.25;
				}
				else
				{
					int distance = Math.abs(row - col);
					Ax[2 * p] = -1.0 / distance;
					Ax[2 * p + 1] = 0.05 * (row < col ? -1.0 : 1.0);
				}
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

	private static double[] rhsForOnesSolution(Matrix matrix, int nrhs)
	{
		double[] rhs = new double[2 * matrix.n * nrhs];
		for (int rhsIndex = 0; rhsIndex < nrhs; rhsIndex++)
		{
			for (int col = 0; col < matrix.n; col++)
			{
				double xr = 1.0 + rhsIndex;
				double xi = -0.25 * rhsIndex;
				for (int p = matrix.Ap[col]; p < matrix.Ap[col + 1]; p++)
				{
					int row = matrix.Ai[p];
					int b = 2 * (rhsIndex * matrix.n + row);
					double ar = matrix.Ax[2 * p];
					double ai = matrix.Ax[2 * p + 1];
					rhs[b] += ar * xr - ai * xi;
					rhs[b + 1] += ai * xr + ar * xi;
				}
			}
		}
		return rhs;
	}

	private static double maxSolutionError(double[] solution, int n, int nrhs)
	{
		double max = 0.0;
		for (int rhsIndex = 0; rhsIndex < nrhs; rhsIndex++)
		{
			double expectedRe = 1.0 + rhsIndex;
			double expectedIm = -0.25 * rhsIndex;
			for (int row = 0; row < n; row++)
			{
				int p = 2 * (rhsIndex * n + row);
				double dr = solution[p] - expectedRe;
				double di = solution[p + 1] - expectedIm;
				max = Math.max(max, Math.sqrt(dr * dr + di * di));
			}
		}
		return max;
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

	private static final class Matrix
	{
		int n;
		int[] Ap;
		int[] Ai;
		double[] Ax;
	}

	private static final class Result
	{
		long analyzeNs;
		long factorNs;
		long refactorNs;
		long solveNs;
		double residual;
		int luEntries;
	}
}
