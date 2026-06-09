package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;

import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_add_size_t;
import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_malloc_dbl;
import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_malloc_int;
import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_mult_size_t;
import static edu.ufl.cise.klu.tdouble.Dklu_memory.klu_realloc_dbl;

public class Zklu_memory extends Zklu_internal
{
	public static int klu_z_add_size_t(int a, int b, int[] ok)
	{
		return klu_add_size_t(a, b, ok);
	}

	public static int klu_z_mult_size_t(int a, int k, int[] ok)
	{
		return klu_mult_size_t(a, k, ok);
	}

	public static int[] klu_z_malloc_int(int n, KLU_common Common)
	{
		return klu_malloc_int(n, Common);
	}

	public static double[] klu_z_malloc_dbl(int n, KLU_common Common)
	{
		return klu_malloc_dbl(2 * n, Common);
	}

	public static double[] klu_z_realloc_dbl(int nnew, int nold,
			double[] p, KLU_common Common)
	{
		return klu_realloc_dbl(2 * nnew, 2 * nold, p, Common);
	}
}
