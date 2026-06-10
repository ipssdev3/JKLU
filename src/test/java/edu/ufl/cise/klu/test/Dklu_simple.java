package edu.ufl.cise.klu.test;

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicReference;

import junit.framework.TestCase;
import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_numeric;
import edu.ufl.cise.klu.common.KLU_symbolic;

import static edu.ufl.cise.klu.tdouble.Dklu_defaults.klu_defaults;
import static edu.ufl.cise.klu.tdouble.Dklu_analyze.klu_analyze;
import static edu.ufl.cise.klu.tdouble.Dklu_factor.klu_factor;
import static edu.ufl.cise.klu.tdouble.Dklu_solve.klu_solve;

public class Dklu_simple extends TestCase {

	private static final double DELTA = 1e-09 ;

	private static int n = 5 ;
	private static int [ ] Ap = {0, 2, 5, 9, 10, 12} ;
	private static int [ ] Ai = {0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
	private static double [ ] Ax = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
	private static double [ ] b = {8., 45., -3., 3., 19.} ;

	/**
	 * a simple KLU demo; solution is x = (1,2,3,4,5)
	 */
	public static void test_klu_simple() {
		int i;
		KLU_symbolic Symbolic;
		KLU_numeric Numeric;
		KLU_common Common = new KLU_common();

		//Dklu_version.NPRINT = false ;
		//Dklu_internal.NDEBUG = false ;

		klu_defaults (Common);
		Symbolic = klu_analyze (n, Ap, Ai, Common);
		Numeric = klu_factor (Ap, Ai, Ax, Symbolic, Common);
		double[] rhs = b.clone();
		klu_solve (Symbolic, Numeric, 5, 1, rhs, 0, Common);

		for (i = 0 ; i < n ; i++) {
			System.out.printf("x [%d] = %g\n", i, rhs [i]) ;
			assertEquals(i + 1.0, rhs [i], DELTA) ;
		}
	}

	public static void test_klu_solve_concurrent_calls_use_independent_workspace() throws Exception {
		final KLU_common Common = new KLU_common();
		klu_defaults(Common);
		final KLU_symbolic Symbolic = klu_analyze(n, Ap, Ai, Common);
		final KLU_numeric Numeric = klu_factor(Ap, Ai, Ax, Symbolic, Common);
		final AtomicReference<Throwable> failure = new AtomicReference<Throwable>();
		final CountDownLatch start = new CountDownLatch(1);
		Thread[] threads = new Thread[8];

		for (int t = 0; t < threads.length; t++) {
			threads[t] = new Thread(new Runnable() {
				@Override
				public void run() {
					try {
						start.await();
						for (int repeat = 0; repeat < 200; repeat++) {
							double[] rhs = b.clone();
							assertEquals(1, klu_solve(Symbolic, Numeric, n, 1, rhs, 0, Common));
							for (int i = 0; i < n; i++) {
								assertEquals(i + 1.0, rhs[i], DELTA);
							}
						}
					} catch (Throwable e) {
						failure.compareAndSet(null, e);
					}
				}
			});
			threads[t].start();
		}

		start.countDown();
		for (int t = 0; t < threads.length; t++) {
			threads[t].join();
		}
		if (failure.get() != null) {
			throw new AssertionError(failure.get());
		}
	}

}
