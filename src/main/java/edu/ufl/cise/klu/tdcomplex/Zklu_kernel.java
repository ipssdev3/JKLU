package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_z_numeric;

/**
 * Complex numeric kernel marker.
 * 
 * The current direct Java complex implementation stores and factors interleaved
 * complex values directly in sparse row LU state on {@link KLU_z_numeric}.
 * Factorization maintains active column-to-row incidence so pivot selection and
 * elimination traverse only rows that currently contain the pivot column. A
 * future optimization can add the real KLU kernel's symbolic DFS and pruning
 * details without changing the public klu_z_* API.
 */
public class Zklu_kernel extends Zklu_internal
{
	private Zklu_kernel()
	{
	}
}
