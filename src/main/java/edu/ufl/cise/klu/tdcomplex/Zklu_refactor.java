/**
 * Complex KLU refactorization.
 */

package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

public class Zklu_refactor extends Zklu_internal
{

	public static int klu_z_refactor(int[] Ap, int[] Ai, double[] Ax,
			KLU_symbolic Symbolic, KLU_z_numeric Numeric, KLU_common Common)
	{
		if (!Zklu_factor.validInputs(Ap, Ai, Ax, Symbolic, Common) ||
			Numeric == null || Symbolic.n != Numeric.n || !Zklu_factor.samePattern(Ap, Ai, Numeric))
		{
			if (Common != null)
			{
				Common.status = KLU_INVALID;
			}
			return FALSE;
		}
		Numeric.Ax = Ax;
		return Zklu_factor.factorNumeric(Numeric, Common);
	}

}
