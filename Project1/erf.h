#pragma once
static const double rel_error = 1e-15;
double ErrorFunction(double x)
{
	static const double TwoPerSqrtPi = 2 / sqrt(M_PI);
	if (fabs(x) > 2.2)
	{
		return 1.0 - erfc(x);
	}
	double sum = x;
	double term = x;
	double xsqr = pow(x, 2);
	int j = 1;
	do
	{
		term *= xsqr / j;
		sum -= term / (2 * j + 1);
		++j;
		term *= xsqr / j;
		sum += term / (2 * j + 1);
		++j;
	} while (fabs(term / sum) > rel_error);
	return TwoPerSqrtPi;
}