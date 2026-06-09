/**
 * Complex helpers for JKLU.
 */

package edu.ufl.cise.klu.tdcomplex;

public class Zklu_complex
{

	public static int offset(int i)
	{
		return 2 * i;
	}

	public static double real(double[] x, int i)
	{
		return x[offset(i)];
	}

	public static double imag(double[] x, int i)
	{
		return x[offset(i) + 1];
	}

	public static void set(double[] x, int i, double re, double im)
	{
		int p = offset(i);
		x[p] = re;
		x[p + 1] = im;
	}

	public static void copy(double[] dst, int di, double[] src, int si)
	{
		int d = offset(di);
		int s = offset(si);
		dst[d] = src[s];
		dst[d + 1] = src[s + 1];
	}

	public static void clear(double[] x, int i)
	{
		set(x, i, 0.0, 0.0);
	}

	public static void one(double[] x, int i)
	{
		set(x, i, 1.0, 0.0);
	}

	public static double abs(double[] x, int i)
	{
		return abs(real(x, i), imag(x, i));
	}

	public static double abs(double re, double im)
	{
		double ar = Math.abs(re);
		double ai = Math.abs(im);
		if (ar == 0.0 && ai == 0.0)
		{
			return 0.0;
		}
		if (ar >= ai)
		{
			double r = ai / ar;
			return ar * Math.sqrt(1.0 + r * r);
		}
		double r = ar / ai;
		return ai * Math.sqrt(1.0 + r * r);
	}

	public static boolean isZero(double[] x, int i)
	{
		return real(x, i) == 0.0 && imag(x, i) == 0.0;
	}

	public static void divide(double[] dst, int di, double[] a, int ai, double[] b, int bi)
	{
		divide(dst, di, real(a, ai), imag(a, ai), real(b, bi), imag(b, bi));
	}

	public static void divide(double[] dst, int di, double ar, double ai, double br, double bi)
	{
		int d = offset(di);
		double scalar;
		if (Math.abs(br) >= Math.abs(bi))
		{
			double ratio = bi / br;
			scalar = 1.0 / (br + bi * ratio);
			dst[d] = scalar * (ar + ai * ratio);
			dst[d + 1] = scalar * (ai - ar * ratio);
		}
		else
		{
			double ratio = br / bi;
			scalar = 1.0 / (br * ratio + bi);
			dst[d] = scalar * (ar * ratio + ai);
			dst[d + 1] = scalar * (ai * ratio - ar);
		}
	}

	public static void multiplySubtract(double[] x, int xi, double ar, double ai, double br, double bi)
	{
		int p = offset(xi);
		x[p] -= ar * br - ai * bi;
		x[p + 1] -= ai * br + ar * bi;
	}

	public static void scaleByReal(double[] x, int xi, double scale)
	{
		int p = offset(xi);
		x[p] /= scale;
		x[p + 1] /= scale;
	}

	public static boolean equals(double[] x, int xi, double re, double im, double tol)
	{
		return abs(real(x, xi) - re, imag(x, xi) - im) <= tol;
	}

}
