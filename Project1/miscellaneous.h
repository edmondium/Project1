#pragma once
// This is the parametrization for the effective potential from experiment (for Cu)
double VeffP(double R)
{
    return 29 * exp(-2.3151241717834 * pow(R, 0.81266614122432) + 2.1984250222603e-2 * pow(R, 4.2246376280056))
        - 0.15595606773483 * R - 3.1350051440417e-3 * R * R + 5.1895222293006e-2 * pow(R, 3) - 2.8027608685637e-2 * pow(R, 4);
}
// Given the input density rho, calculates the Hartree potential
// The boundary conditions used are U(0)=0 and U(S)=Zq. The boundary condition at S is only a constant shift
// of potential which is later readjusted by choosing MT zero. So, this boundary condition is not very relevant.
void SolvePoisson(int Zq, const vector<double>& Rmesh, const function1D<double>& rho, vector<double>& Uhartree)
{
    static vector<double> RHS(Rmesh.size());
    for (int i = 0; i < Rmesh.size(); i++)
    {
        RHS[i] = -4 * M_PI * Rmesh[i] * rho[i];
    }
    Uhartree[0] = 0; 
    Uhartree[1] = (Rmesh[1] - Rmesh[0]);
    //NumerovInhom(RHS, RHS.size(), Rmesh[1] - Rmesh[0], Uhartree);
    NumerovInhom(RHS, RHS.size(), Uhartree[1], Uhartree);
    int ilast = Uhartree.size() - 1;
    double U_last = Uhartree[ilast];
    double alpha = (Zq - U_last) / Rmesh[ilast];
    for (int i = 0; i < Rmesh.size(); i++)
    {
        Uhartree[i] += alpha * Rmesh[i];
    }
}
double IntegrateCharge(const vector<double>& R, const function1D<double>& rho)
{
    static function1D<double> temp(R.size());
    for (int i = 0; i < R.size(); i++)
    {
        temp[i] = rho[i] * sqr(R[i]) * 4 * M_PI;
    }
    return integrate4<double>(temp, R[1] - R[0], temp.size());
}
