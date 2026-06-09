package edu.ufl.cise.klu.test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import edu.ufl.cise.klu.bench.ZkluMatrixMarketBenchmark;
import edu.ufl.cise.klu.solver.KluComplexSparseSet;
import edu.ufl.cise.klu.solver.KluComplexSparseSet.CscMatrix;
import edu.ufl.cise.klu.solver.KluComplexSparseSet.TripletMatrix;
import junit.framework.TestCase;

public class Zklu_sparse_set extends TestCase
{
	private static final double DELTA = 1e-09;

	public void test_sparse_set_solve_complex_system()
	{
		KluComplexSparseSet set = new KluComplexSparseSet(2);
		assertTrue(set.setMatrixElement(1, 1, 2.0, 0.0));
		assertTrue(set.setMatrixElement(2, 1, 1.0, -1.0));
		assertTrue(set.setMatrixElement(1, 2, 0.0, 1.0));
		assertTrue(set.setMatrixElement(2, 2, 3.0, 0.0));

		assertEquals(KluComplexSparseSet.OK, set.factor());
		double[] solution = new double[4];
		assertEquals(KluComplexSparseSet.OK,
			set.solve(new double[] {4.0, 1.0, 5.0, -1.0}, solution));

		assertComplexEquals(2.1923076923076925, 0.03846153846153838, solution, 0);
		assertComplexEquals(0.9230769230769232, 0.38461538461538464, solution, 1);
	}

	public void test_add_matrix_element_is_symmetric()
	{
		KluComplexSparseSet set = new KluComplexSparseSet(2);
		assertTrue(set.addMatrixElement(1, 2, 1.0, -1.0));

		assertComplexEquals(1.0, -1.0, set.getMatrixElement(1, 2), 0);
		assertComplexEquals(1.0, -1.0, set.getMatrixElement(2, 1), 0);
		assertEquals(2, set.getNNZ());
	}

	public void test_primitive_matrix_skips_ground_node()
	{
		KluComplexSparseSet set = new KluComplexSparseSet(2);
		int[] nodes = {0, 1, 2};
		double[] primitive = {
			99.0, 0.0, 99.0, 0.0, 99.0, 0.0,
			99.0, 0.0, 2.0, 1.0, 3.0, -1.0,
			99.0, 0.0, 4.0, 2.0, 5.0, -2.0
		};

		assertTrue(set.addPrimitiveMatrix(nodes, primitive));

		assertComplexEquals(0.0, 0.0, set.getMatrixElement(0, 1), 0);
		assertComplexEquals(2.0, 1.0, set.getMatrixElement(1, 1), 0);
		assertComplexEquals(3.0, -1.0, set.getMatrixElement(1, 2), 0);
		assertComplexEquals(4.0, 2.0, set.getMatrixElement(2, 1), 0);
		assertComplexEquals(5.0, -2.0, set.getMatrixElement(2, 2), 0);
	}

	public void test_increment_zeroise_and_refactor_same_pattern()
	{
		KluComplexSparseSet set = new KluComplexSparseSet(2);
		assertTrue(set.setMatrixElement(1, 1, 2.0, 0.0));
		assertTrue(set.setMatrixElement(1, 2, 1.0, 0.0));
		assertTrue(set.setMatrixElement(2, 1, 1.0, 0.0));
		assertTrue(set.setMatrixElement(2, 2, 3.0, 0.0));
		assertEquals(KluComplexSparseSet.OK, set.factor());

		assertTrue(set.incrementMatrixElement(2, 2, 1.0, 0.0));
		assertTrue(set.zeroiseMatrixElement(1, 2));
		assertEquals(KluComplexSparseSet.OK, set.refactor());

		double[] solution = new double[4];
		assertEquals(KluComplexSparseSet.OK,
			set.solve(new double[] {2.0, 0.0, 9.0, 0.0}, solution));
		assertComplexEquals(1.0, 0.0, solution, 0);
		assertComplexEquals(2.0, 0.0, solution, 1);
	}

	public void test_extraction_export_and_islands()
			throws Exception
	{
		KluComplexSparseSet set = new KluComplexSparseSet(3);
		assertTrue(set.setMatrixElement(1, 1, 2.0, 0.0));
		assertTrue(set.setMatrixElement(2, 2, 3.0, 0.0));
		assertTrue(set.setMatrixElement(3, 3, 4.0, 0.0));
		assertTrue(set.addMatrixElement(1, 2, 1.0, 0.5));

		CscMatrix csc = set.getCompressedMatrix();
		TripletMatrix triplet = set.getTripletMatrix();
		assertEquals(3, csc.n);
		assertEquals(5, csc.colPointers[3]);
		assertEquals(5, triplet.rows.length);
		assertEquals(10, triplet.values.length);

		int[] islands = set.findIslands();
		assertEquals(islands[0], islands[1]);
		assertTrue(islands[0] != islands[2]);

		File matrixFile = new File("target/Zklu_sparse_set_matrix.mtx");
		File vectorFile = new File("target/Zklu_sparse_set_vector.mtx");
			assertTrue(set.saveAsMatrixMarket(matrixFile.getPath(),
				new double[] {1.0, 0.0, 2.0, -1.0, 3.0, 1.0}, vectorFile.getPath()));
			assertFirstLine("%%MatrixMarket matrix coordinate complex general", matrixFile);
			assertFirstLine("%%MatrixMarket matrix array complex general", vectorFile);

			ZkluMatrixMarketBenchmark.Matrix matrix =
					ZkluMatrixMarketBenchmark.readMatrix(matrixFile.getPath());
			double[] rhs = ZkluMatrixMarketBenchmark.readVector(vectorFile.getPath(), matrix.n);
			ZkluMatrixMarketBenchmark.Result result =
					ZkluMatrixMarketBenchmark.run(matrix, rhs, 0, 1);
			assertTrue("Matrix Market replay residual " + result.maxResidual,
					result.maxResidual < 1.0e-9);
		}

	private static void assertComplexEquals(double re, double im, double[] values, int index)
	{
		assertEquals(re, values[2 * index], DELTA);
		assertEquals(im, values[2 * index + 1], DELTA);
	}

	private static void assertFirstLine(String expected, File file)
			throws Exception
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		try
		{
			assertEquals(expected, reader.readLine());
		}
		finally
		{
			reader.close();
		}
	}
}
