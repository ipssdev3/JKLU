package edu.ufl.cise.klu.tdcomplex;

import edu.ufl.cise.klu.common.KLU_common;
import edu.ufl.cise.klu.common.KLU_symbolic;
import edu.ufl.cise.klu.common.KLU_z_numeric;

public class Zklu_sort extends Zklu_internal
{
	public static int klu_z_sort(KLU_symbolic Symbolic, KLU_z_numeric Numeric,
			KLU_common Common)
	{
		if (Common == null)
		{
			return FALSE;
		}
		if (Symbolic == null || Numeric == null || Numeric.LUrowCols == null)
		{
			Common.status = KLU_INVALID;
			return FALSE;
		}
		Common.status = KLU_OK;
		return TRUE;
	}
}
