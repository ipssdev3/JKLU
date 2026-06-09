package edu.ufl.cise.klu.test;

import junit.framework.TestCase;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcs;
import edu.emory.mathcs.csparsej.tdcomplex.DZcs_common.DZcsa;
import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_compress.cs_compress;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_entry.cs_entry;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_lusol.cs_lusol;
import static edu.emory.mathcs.csparsej.tdcomplex.DZcs_util.cs_spalloc;
import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdouble.Dklu_version.KLU_INVALID;
import static edu.ufl.cise.klu.tdcomplex.Zklu_diagnostics.klu_z_condest;
import static edu.ufl.cise.klu.tdcomplex.Zklu_diagnostics.klu_z_flops;
import static edu.ufl.cise.klu.tdcomplex.Zklu_diagnostics.klu_z_rcond;
import static edu.ufl.cise.klu.tdcomplex.Zklu_diagnostics.klu_z_rgrowth;
import static edu.ufl.cise.klu.tdcomplex.Zklu_extract.klu_z_extract;
import static edu.ufl.cise.klu.tdcomplex.Zklu_factor.klu_z_factor;
import static edu.ufl.cise.klu.tdcomplex.Zklu_refactor.klu_z_refactor;
import static edu.ufl.cise.klu.tdcomplex.Zklu_solve.klu_z_solve;
import static edu.ufl.cise.klu.tdcomplex.Zklu_tsolve.klu_z_tsolve;

public class Zklu_simple extends TestCase {

	private static final double DELTA = 1e-09;

	public void test_klu_z_real_matrix_matches_real_solution() {
		int n = 5;
		int[] Ap = {0, 2, 5, 9, 10, 12};
		int[] Ai = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4};
		double[] Ax = interleaved(new double[] {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.});
		double[] b = {8., 0., 45., 0., -3., 0., 3., 0., 19., 0.};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);

		assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
		for (int i = 0; i < n; i++) {
			assertEquals(i + 1.0, b[2 * i], DELTA);
			assertEquals(0.0, b[2 * i + 1], DELTA);
		}
	}

	public void test_klu_z_complex_diagonal() {
		int n = 2;
		int[] Ap = {0, 1, 2};
		int[] Ai = {0, 1};
		double[] Ax = {1., 1., 2., -1.};
		double[] b = {2., 4., 3., 1.};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		Common.scale = 0;
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);

		assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
		assertEquals(3.0, b[0], DELTA);
		assertEquals(1.0, b[1], DELTA);
		assertEquals(1.0, b[2], DELTA);
		assertEquals(1.0, b[3], DELTA);
	}

	public void test_klu_z_multiple_rhs() {
		int n = 2;
		int[] Ap = {0, 2, 4};
		int[] Ai = {0, 1, 0, 1};
		double[] Ax = {
			2., 0.,
			1., -1.,
			0., 1.,
			3., 0.
		};
		double[] b = new double[2 * n * 2];
		set(b, n, 0, 0, 4., 1.);
		set(b, n, 0, 1, 5., -1.);
		set(b, n, 1, 0, 1., 2.);
		set(b, n, 1, 1, -3., 4.);

		double[] expected = b.clone();
		solveDense2(Ax, expected, 0);
		solveDense2(Ax, expected, 1);

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);

		assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 2, b, Common));
		for (int i = 0; i < b.length; i++) {
			assertEquals(expected[i], b[i], DELTA);
		}
	}

	public void test_klu_z_matches_csparsej_complex_lu() {
		int n = 2;
		int[] Ap = {0, 2, 4};
		int[] Ai = {0, 1, 0, 1};
		double[] Ax = {
			2., 0.,
			1., -1.,
			0., 1.,
			3., 0.
		};
		double[] jkB = {4., 1., 5., -1.};
		double[] csB = jkB.clone();

		DZcs triplet = cs_spalloc(n, n, Ap[n], true, true);
		for (int col = 0; col < n; col++) {
			for (int p = Ap[col]; p < Ap[col + 1]; p++) {
				assertTrue(cs_entry(triplet, Ai[p], col, Ax[2 * p], Ax[2 * p + 1]));
			}
		}
		DZcs csc = cs_compress(triplet);
		DZcsa rhs = new DZcsa(csB);
		assertTrue(cs_lusol(2, csc, rhs, 1.0));

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);
		assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, jkB, Common));

		for (int i = 0; i < jkB.length; i++) {
			assertEquals(csB[i], jkB[i], DELTA);
		}
		assertResidualBelow(n, Ap, Ai, Ax, new double[] {4., 1., 5., -1.}, jkB, 1e-12);
	}

	public void test_klu_z_refactor_same_pattern() {
		int n = 2;
		int[] Ap = {0, 1, 2};
		int[] Ai = {0, 1};
		double[] Ax = {1., 1., 2., -1.};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);

		double[] Ax2 = {0., 2., 4., 0.};
		assertEquals(1, klu_z_refactor(Ap, Ai, Ax2, Symbolic, Numeric, Common));

		double[] b = {-2., 6., 8., -4.};
		assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
		assertEquals(3.0, b[0], DELTA);
		assertEquals(1.0, b[1], DELTA);
		assertEquals(2.0, b[2], DELTA);
		assertEquals(-1.0, b[3], DELTA);
	}

	public void test_klu_z_direct_flat_refactor_solve_and_extract_use_new_values() {
		String oldMinN = System.getProperty("jklu.complex.factor.columnKernelMinN");
		String oldPrecompute = System.getProperty(
				"jklu.complex.factor.precomputeUdiagInv");
		try {
			System.setProperty("jklu.complex.factor.columnKernelMinN", "0");
			System.setProperty("jklu.complex.factor.precomputeUdiagInv", "true");
			int n = 2;
			int[] Ap = {0, 2, 4};
			int[] Ai = {0, 1, 0, 1};
			double[] Ax = {
				2., 0.,
				1., 0.,
				1., 0.,
				3., 0.
			};

			KLU_common Common = new KLU_common();
			klu_defaults(Common);
			Common.scale = 0;
			Common.btf = 0;
			Common.ordering = 2;
			KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
			KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
			assertNotNull(Numeric);
			assertNotNull(Numeric.DirectLp);
			assertNotNull(Numeric.UdiagInv);

			double[] Ax2 = {
				4., 0.,
				2., 0.,
				1., 0.,
				5., 0.
			};
			assertEquals(1, klu_z_refactor(Ap, Ai, Ax2, Symbolic, Numeric, Common));
			assertEquals(0.25, Numeric.UdiagInv[0], DELTA);
			assertEquals(0.0, Numeric.UdiagInv[1], DELTA);
			assertEquals(1.0 / 4.5, Numeric.UdiagInv[2], DELTA);
			assertEquals(0.0, Numeric.UdiagInv[3], DELTA);

			double[] b = {11., 0., 19., 0.};
			assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
			assertEquals(2.0, b[0], DELTA);
			assertEquals(0.0, b[1], DELTA);
			assertEquals(3.0, b[2], DELTA);
			assertEquals(0.0, b[3], DELTA);

			int[] Up = new int[n + 1];
			int[] Ui = new int[Numeric.unz];
			double[] Ux = new double[2 * Numeric.unz];
			assertEquals(1, klu_z_extract(Numeric, Symbolic, null, null, null,
					Up, Ui, Ux, null, null, null, null, null, null, null, Common));
			assertEquals(4.0, findExtractedValue(Up, Ui, Ux, 0, 0, true), DELTA);
			assertEquals(1.0, findExtractedValue(Up, Ui, Ux, 1, 0, true), DELTA);
			assertEquals(4.5, findExtractedValue(Up, Ui, Ux, 1, 1, true), DELTA);
		} finally {
			restoreProperty("jklu.complex.factor.columnKernelMinN", oldMinN);
			restoreProperty("jklu.complex.factor.precomputeUdiagInv", oldPrecompute);
		}
	}

	public void test_klu_z_forced_direct_threshold_solves_off_diagonal_structure() {
		String oldMinN = System.getProperty("jklu.complex.factor.columnKernelMinN");
		try {
			System.setProperty("jklu.complex.factor.columnKernelMinN", "0");
			int n = 2;
			int[] Ap = {0, 1, 3};
			int[] Ai = {1, 0, 1};
			double[] Ax = {
				1., 0.,
				1., 0.,
				1., 0.
			};
			double[] b = {2., 0., 3., 0.};

			KLU_common Common = new KLU_common();
			klu_defaults(Common);
			Common.scale = 0;
			KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
			KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
			assertNotNull(Numeric);

			assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
			assertEquals(1.0, b[0], DELTA);
			assertEquals(0.0, b[1], DELTA);
			assertEquals(2.0, b[2], DELTA);
			assertEquals(0.0, b[3], DELTA);
		} finally {
			restoreProperty("jklu.complex.factor.columnKernelMinN", oldMinN);
		}
	}

	public void test_klu_z_tsolve_complex_system() {
		int n = 2;
		int[] Ap = {0, 2, 4};
		int[] Ai = {0, 1, 0, 1};
		double[] Ax = {
			2., 0.,
			1., -1.,
			0., 1.,
			3., 0.
		};
		double[] b = {3., -1., 5., -2.};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);

		assertEquals(1, klu_z_tsolve(Symbolic, Numeric, n, 1, b, Common));
		assertEquals(1.0, b[0], DELTA);
		assertEquals(1.0, b[1], DELTA);
		assertEquals(2.0, b[2], DELTA);
		assertEquals(-1.0, b[3], DELTA);
	}

	public void test_klu_z_row_scaling_direct_complex() {
		int n = 2;
		int[] Ap = {0, 2, 4};
		int[] Ai = {0, 1, 0, 1};
		double[] Ax = {
			100., 0.,
			1., -1.,
			0., 1.,
			3., 0.
		};
		double[] b = {200., -99., 4., -3.};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		Common.scale = 2;
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);
		assertNotNull(Numeric.Rs);

		assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
		assertEquals(2.0, b[0], DELTA);
		assertEquals(-1.0, b[1], DELTA);
		assertEquals(1.0, b[2], DELTA);
		assertEquals(0.0, b[3], DELTA);
	}

	public void test_klu_z_extract_and_diagnostics_use_original_complex_size() {
		int n = 2;
		int[] Ap = {0, 2, 4};
		int[] Ai = {0, 1, 0, 1};
		double[] Ax = {
			2., 0.,
			1., -1.,
			0., 1.,
			3., 0.
		};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);
		assertEquals(n, Numeric.n);
		assertEquals(1, Numeric.nblocks);

		int packedTriangular = n * (n + 1) / 2;
		int[] Lp = new int[n + 1];
		int[] Li = new int[packedTriangular];
		double[] Lx = new double[2 * packedTriangular];
		int[] Up = new int[n + 1];
		int[] Ui = new int[packedTriangular];
		double[] Ux = new double[2 * packedTriangular];
		int[] Fp = new int[n + 1];
		int[] P = new int[n];
		int[] Q = new int[n];
		int[] R = new int[Symbolic.nblocks + 1];

		assertEquals(1, klu_z_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux,
				Fp, new int[0], new double[0], P, Q, null, R, Common));
		assertEquals(packedTriangular, Lp[n]);
		assertEquals(packedTriangular, Up[n]);
		assertEquals(0, Fp[n]);

		assertEquals(1, klu_z_flops(Symbolic, Numeric, Common));
		assertTrue(Common.flops > 0.0);
		assertEquals(1, klu_z_rcond(Symbolic, Numeric, Common));
		assertTrue(Common.rcond >= 0.0);
		assertEquals(1, klu_z_rgrowth(Ap, Ai, Ax, Symbolic, Numeric, Common));
		assertTrue(Common.rgrowth >= 0.0);
		assertEquals(1, klu_z_condest(Ap, Ax, Symbolic, Numeric, Common));
		assertTrue(Common.condest >= 0.0);
	}

	public void test_klu_z_btf_off_diagonal_entries_are_solved_and_extracted() {
		int n = 3;
		int[] Ap = {0, 1, 2, 4};
		int[] Ai = {0, 1, 0, 2};
		double[] Ax = {
			2., 0.,
			3., 0.,
			5., 1.,
			4., 0.
		};
		double[] b = {
			11., 1.,
			6., 0.,
			4., 0.
		};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		Common.scale = 0;
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);
		assertTrue(Numeric.nblocks > 1);
		assertEquals(1, Numeric.nzoff);

		assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
		assertEquals(3.0, b[0], DELTA);
		assertEquals(0.0, b[1], DELTA);
		assertEquals(2.0, b[2], DELTA);
		assertEquals(0.0, b[3], DELTA);
		assertEquals(1.0, b[4], DELTA);
		assertEquals(0.0, b[5], DELTA);

		int[] Fp = new int[n + 1];
		int[] Fi = new int[Numeric.nzoff];
		double[] Fx = new double[2 * Numeric.nzoff];
		assertEquals(1, klu_z_extract(Numeric, Symbolic, null, null, null, null, null, null,
				Fp, Fi, Fx, null, null, null, null, Common));
		assertEquals(1, Fp[n]);
		assertEquals(0, Fi[0]);
		assertEquals(5.0, Fx[0], DELTA);
		assertEquals(1.0, Fx[1], DELTA);
	}

	public void test_klu_z_direct_btf_off_diagonal_entries_refactor_and_extract() {
		String oldMinN = System.getProperty("jklu.complex.factor.columnKernelMinN");
		try {
			System.setProperty("jklu.complex.factor.columnKernelMinN", "0");
			int n = 3;
			int[] Ap = {0, 1, 2, 4};
			int[] Ai = {0, 1, 0, 2};
			double[] Ax = {
				2., 0.,
				3., 0.,
				5., 1.,
				4., 0.
			};

			KLU_common Common = new KLU_common();
			klu_defaults(Common);
			Common.scale = 0;
			KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
			KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
			assertNotNull(Numeric);
			assertNotNull(Numeric.DirectLp);
			assertTrue(Numeric.nblocks > 1);
			assertEquals(1, Numeric.nzoff);

			double[] Ax2 = {
				2., 0.,
				3., 0.,
				7., -2.,
				4., 0.
			};
			assertEquals(1, klu_z_refactor(Ap, Ai, Ax2, Symbolic, Numeric, Common));

			double[] b = {
				13., -2.,
				6., 0.,
				4., 0.
			};
			assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
			assertEquals(3.0, b[0], DELTA);
			assertEquals(0.0, b[1], DELTA);
			assertEquals(2.0, b[2], DELTA);
			assertEquals(0.0, b[3], DELTA);
			assertEquals(1.0, b[4], DELTA);
			assertEquals(0.0, b[5], DELTA);

			int[] Fp = new int[n + 1];
			int[] Fi = new int[Numeric.nzoff];
			double[] Fx = new double[2 * Numeric.nzoff];
			assertEquals(1, klu_z_extract(Numeric, Symbolic, null, null, null,
					null, null, null, Fp, Fi, Fx, null, null, null, null, Common));
			assertEquals(1, Fp[n]);
			assertEquals(0, Fi[0]);
			assertEquals(7.0, Fx[0], DELTA);
			assertEquals(-2.0, Fx[1], DELTA);
		} finally {
			restoreProperty("jklu.complex.factor.columnKernelMinN", oldMinN);
		}
	}

	public void test_klu_z_fragmented_btf_direct_kernel_solves_offdiagonal() {
		String oldMinN = System.getProperty("jklu.complex.factor.columnKernelMinN");
		String oldDisable = System.getProperty(
				"jklu.complex.factor.disableFragmentedBtfColumnKernel");
		try {
			System.setProperty("jklu.complex.factor.columnKernelMinN", "0");
			System.clearProperty("jklu.complex.factor.disableFragmentedBtfColumnKernel");
			int n = 3;
			int[] Ap = {0, 1, 2, 4};
			int[] Ai = {0, 1, 0, 2};
			double[] Ax = {
				2., 0.,
				3., 0.,
				5., 1.,
				4., 0.
			};
			double[] b = {
				11., 1.,
				6., 0.,
				4., 0.
			};

			KLU_common Common = new KLU_common();
			klu_defaults(Common);
			Common.scale = 0;
			KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
			KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
			assertNotNull(Numeric);
			assertTrue(Numeric.nblocks > 1);
			assertNotNull(Numeric.DirectLp);
			assertNotNull(Numeric.DirectUp);

			assertEquals(1, klu_z_solve(Symbolic, Numeric, n, 1, b, Common));
			assertEquals(3.0, b[0], DELTA);
			assertEquals(0.0, b[1], DELTA);
			assertEquals(2.0, b[2], DELTA);
			assertEquals(0.0, b[3], DELTA);
			assertEquals(1.0, b[4], DELTA);
			assertEquals(0.0, b[5], DELTA);
		} finally {
			restoreProperty("jklu.complex.factor.columnKernelMinN", oldMinN);
			restoreProperty("jklu.complex.factor.disableFragmentedBtfColumnKernel",
					oldDisable);
		}
	}

	public void test_klu_z_numeric_uses_sparse_complex_storage() {
		int n = 5;
		int[] Ap = {0, 1, 2, 3, 4, 5};
		int[] Ai = {0, 1, 2, 3, 4};
		double[] Ax = {
			1., 1.,
			2., -1.,
			3., 0.,
			4., 2.,
			5., -2.
		};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);
		assertNotNull(Numeric.DirectLp);
		assertNotNull(Numeric.DirectUp);
		assertNotNull(Numeric.Udiag);

		int stored = Numeric.DirectLp[n] + Numeric.DirectUp[n] + Numeric.Udiag.length / 2;
		assertEquals(n, stored);
		assertTrue("sparse LU storage should not allocate n*n entries", stored < n * n);
	}

	public void test_klu_z_extract_uses_actual_sparse_lu_sizes() {
		int n = 5;
		int[] Ap = {0, 1, 2, 3, 4, 5};
		int[] Ai = {0, 1, 2, 3, 4};
		double[] Ax = {
			1., 1.,
			2., -1.,
			3., 0.,
			4., 2.,
			5., -2.
		};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);
		assertEquals(n, Numeric.lnz);
		assertEquals(n, Numeric.unz);

		int[] Lp = new int[n + 1];
		int[] Li = new int[Numeric.lnz];
		double[] Lx = new double[2 * Numeric.lnz];
		int[] Up = new int[n + 1];
		int[] Ui = new int[Numeric.unz];
		double[] Ux = new double[2 * Numeric.unz];
		assertEquals(1, klu_z_extract(Numeric, Symbolic, Lp, Li, Lx, Up, Ui, Ux,
				null, null, null, null, null, null, null, Common));
		assertEquals(Numeric.lnz, Lp[n]);
		assertEquals(Numeric.unz, Up[n]);
		for (int i = 0; i < n; i++) {
			assertEquals(i, Li[Lp[i]]);
			assertEquals(1.0, Lx[2 * Lp[i]], DELTA);
			assertEquals(0.0, Lx[2 * Lp[i] + 1], DELTA);
			assertEquals(i, Ui[Up[i]]);
		}
	}

	public void test_klu_z_extract_rejects_undersized_sparse_lu_arrays() {
		int n = 2;
		int[] Ap = {0, 1, 2};
		int[] Ai = {0, 1};
		double[] Ax = {1., 0., 2., 0.};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);

		assertEquals(0, klu_z_extract(Numeric, Symbolic, new int[n + 1],
				new int[Numeric.lnz - 1], new double[2 * Numeric.lnz],
				null, null, null, null, null, null, null, null, null, null, Common));
		assertEquals(KLU_INVALID, Common.status);
	}

	public void test_klu_z_factor_rejects_malformed_csc() {
		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(2, new int[] {0, 1, 2}, new int[] {0, 1}, Common);
		assertNotNull(Symbolic);

		assertNull(klu_z_factor(new int[] {0, 2, 1}, new int[] {0, 1},
				new double[] {1., 0., 2., 0.}, Symbolic, Common));
		assertEquals(KLU_INVALID, Common.status);

		klu_defaults(Common);
		assertNull(klu_z_factor(new int[] {0, 1, 2}, new int[] {0, 3},
				new double[] {1., 0., 2., 0.}, Symbolic, Common));
		assertEquals(KLU_INVALID, Common.status);
	}

	public void test_klu_z_solve_rejects_undersized_rhs_buffer() {
		int n = 2;
		int[] Ap = {0, 1, 2};
		int[] Ai = {0, 1};
		double[] Ax = {1., 0., 2., 0.};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);

		assertEquals(0, klu_z_solve(Symbolic, Numeric, n, 1, new double[2 * n - 1],
				Common));
		assertEquals(KLU_INVALID, Common.status);

		assertEquals(0, klu_z_solve(Symbolic, Numeric, n, 1, new double[2 * n],
				1, Common));
		assertEquals(KLU_INVALID, Common.status);
	}

	public void test_klu_z_rgrowth_rejects_malformed_matrix() {
		int n = 2;
		int[] Ap = {0, 1, 2};
		int[] Ai = {0, 1};
		double[] Ax = {1., 0., 2., 0.};

		KLU_common Common = new KLU_common();
		klu_defaults(Common);
		KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		KLU_z_numeric Numeric = klu_z_factor(Ap, Ai, Ax, Symbolic, Common);
		assertNotNull(Numeric);

		assertEquals(0, klu_z_rgrowth(new int[] {0, 2, 1}, Ai, Ax,
				Symbolic, Numeric, Common));
		assertEquals(KLU_INVALID, Common.status);
	}

	private static double[] interleaved(double[] realValues) {
		double[] z = new double[2 * realValues.length];
		for (int i = 0; i < realValues.length; i++) {
			z[2 * i] = realValues[i];
			z[2 * i + 1] = 0.0;
		}
		return z;
	}

	private static void set(double[] b, int d, int rhs, int row, double re, double im) {
		int p = 2 * (rhs * d + row);
		b[p] = re;
		b[p + 1] = im;
	}

	private static void solveDense2(double[] Ax, double[] b, int rhs) {
		int d = 2;
		int b0 = 2 * (rhs * d);
		int b1 = b0 + 2;

		double a00r = Ax[0], a00i = Ax[1];
		double a10r = Ax[2], a10i = Ax[3];
		double a01r = Ax[4], a01i = Ax[5];
		double a11r = Ax[6], a11i = Ax[7];
		double detR = a00r * a11r - a00i * a11i - a01r * a10r + a01i * a10i;
		double detI = a00r * a11i + a00i * a11r - a01r * a10i - a01i * a10r;

		double x0r = b[b0] * a11r - b[b0 + 1] * a11i - a01r * b[b1] + a01i * b[b1 + 1];
		double x0i = b[b0] * a11i + b[b0 + 1] * a11r - a01r * b[b1 + 1] - a01i * b[b1];
		double x1r = a00r * b[b1] - a00i * b[b1 + 1] - b[b0] * a10r + b[b0 + 1] * a10i;
		double x1i = a00r * b[b1 + 1] + a00i * b[b1] - b[b0] * a10i - b[b0 + 1] * a10r;

		double den = detR * detR + detI * detI;
		b[b0] = (x0r * detR + x0i * detI) / den;
		b[b0 + 1] = (x0i * detR - x0r * detI) / den;
		b[b1] = (x1r * detR + x1i * detI) / den;
		b[b1 + 1] = (x1i * detR - x1r * detI) / den;
	}

	private static void assertResidualBelow(int n, int[] Ap, int[] Ai, double[] Ax,
			double[] b, double[] x, double tol) {
		double max = 0.0;
		for (int row = 0; row < n; row++) {
			double rr = b[2 * row];
			double ri = b[2 * row + 1];
			for (int col = 0; col < n; col++) {
				for (int p = Ap[col]; p < Ap[col + 1]; p++) {
					if (Ai[p] == row) {
						double ar = Ax[2 * p];
						double ai = Ax[2 * p + 1];
						double xr = x[2 * col];
						double xi = x[2 * col + 1];
						rr -= ar * xr - ai * xi;
						ri -= ai * xr + ar * xi;
					}
				}
			}
			max = Math.max(max, Math.sqrt(rr * rr + ri * ri));
		}
		assertTrue("complex residual " + max, max < tol);
	}

	private static double findExtractedValue(int[] Ap, int[] Ai, double[] Ax,
			int col, int row, boolean real) {
		for (int p = Ap[col]; p < Ap[col + 1]; p++) {
			if (Ai[p] == row) {
				return Ax[2 * p + (real ? 0 : 1)];
			}
		}
		return Double.NaN;
	}

	private static void restoreProperty(String key, String value) {
		if (value == null) {
			System.clearProperty(key);
		} else {
			System.setProperty(key, value);
		}
	}

}
