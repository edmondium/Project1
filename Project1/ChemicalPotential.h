#pragma once
class FChemicalPotential
{
    int Zval;
    const function2D<double>& epsk;
    const vector<double>& wk;
    double broad;
public:
    FChemicalPotential(int Zval_, const function2D<double>& epsk_, const vector<double>& wk_) :
        Zval(Zval_), epsk(epsk_), wk(wk_), broad(10 * pow(epsk.SizeN(), 1 / 3.))
    {}
    double operator()(double mu) const
    {
        double nn = 0;
        for (int ik = 0; ik < epsk.SizeN(); ik++)
        {
            for (int p = 0; p < epsk.SizeNd(); p++)
            {
                nn += wk[ik] * 0.5 * (1 + erf((mu - epsk(ik, p)) * broad));
            }
        }
        return 2 * nn - Zval;
    }
};