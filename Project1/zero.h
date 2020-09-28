#pragma once
#ifndef _ZEROIN
#define _ZEROIN
#define DBL_EPSILON 2.22045e-16
#include <cassert>
template <class functor>
inline double zeroin
(
	const double ax,
	const double bx,
	functor& f,
	const double tol
)
{
	if (tol > 0)
	{
		assert("Tolerance must be positive");
	}
	if (bx > ax)
	{
		assert("Left end point of the interval should be strictly less than the right one");
	}
	double b = bx;
	double fb = f(b);
	double a = ax;
	double fa = f(a);
	double c = a;
	double fc = fa;
	for (;;)
	{
		const double prev_step = b - a;
		if (abs(fc) < abs(fb))
		{
			a = b;
			b = c;
			c = a;
			fa = fb;
			fb = fc;
			fc = fa;
		}
		const double tol_act = 2 * DBL_EPSILON * abs(b) + tol / 2;
		double new_step = (c - b) / 2;
		if (abs(new_step) <= tol_act || fb == 0)
		{
			return b;
		}
		if (abs(prev_step) >= tol_act && abs(fa) > abs(fb))
		{
			double p;
			double q;
			const double cb = c - b;
			if (a == c)
			{
				const double t1 = fb / fa;
				p = cb * t1;
				q = 1.0 - t1;
			}
			else
			{
				const double t1 = fb / fc;
				const double t2 = fb / fa;
				q = fa / fc;
				p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
				q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
			}
			if (p > 0)
			{
				q = -q;
			}
			else
			{
				p = -p;
			}
			if (2 * p < (1.5 * cb * q - abs(tol_act * q)) && 2 * p < abs(prev_step * q))
			{
				new_step = p / q;
			}
		}
		if (abs(new_step) < tol_act)
		{
			new_step = new_step > 0 ? tol_act : -tol_act;
		}
		a = b;
		fa = fb;
		b += new_step;
		fb = f(b);
		if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0))
		{
			c = a;
			fc = fa;
		}
	}
}
#endif // !_ZEROIN