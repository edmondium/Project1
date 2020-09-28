#pragma once
double extrapolate(double f1, double f2, double x0, double x1, double x2)
{
    return (f1 * (x2 - x0) - f2 * (x1 - x0)) / (x2 - x1);
}
class PartialWave
{
    int Z;
    vector<double> Rmesh;
    vector<double> rhs_MT;
    vector<double> ur;
    vector<double> urp;
    vector<double> temp;
    vector<double> inhom;
    vector<double> dlogPsi;
    vector<double> dlogPsip;
    vector<double> Psi;
    vector<double> Psip;
    vector<double> PsipPsip;
    vector<double> Veff0;
    vector<vector<double> > Psi_l;
    vector<vector<double> > Psip_l;
public:
    PartialWave(int N, double RMuffinTin, int Z_, int lMax) : 
        Z(Z_), Rmesh(N), rhs_MT(N), ur(N), urp(N), temp(N), inhom(N), dlogPsi(lMax + 1), dlogPsip(lMax + 1), Psi(lMax + 1), Psip(lMax + 1), PsipPsip(lMax + 1), Psi_l(lMax + 1), Psip_l(lMax + 1), Veff0(N)
    {
        double dh = RMuffinTin / (N - 1.);
        for (int i = 0; i < N; i++)
        {
            Rmesh[i] = i * dh;
        }
    }
    int Rsize() const
    { 
        return Rmesh.size();
    }
    double R(int i) const
    { 
        return Rmesh[i];
    }
    const vector<double> R() const
    { 
        return Rmesh;
    }
    double psi(int l) const
    { 
        return Psi[l];
    }
    double psip(int l) const
    { 
        return Psip[l];
    }
    double psi(int l, int r) const
    { 
        return Psi_l[l][r];
    }
    double psip(int l, int r) const
    { 
        return Psip_l[l][r];
    }
    double dlog(int l) const
    { 
        return dlogPsi[l];
    }
    double dlogp(int l) const
    { 
        return dlogPsip[l];
    }
    double PP(int l) const
    { 
        return PsipPsip[l];
    }
    double startSol(int Z, int l, double r) const
    {
        return pow(r, l + 1) * (1 - Z * r / (l + 1));
    }
    void SetVeff0(const vector<double>& Veff)
    {
        for (int i = 0; i < Rmesh.size(); i++)
        {
            Veff0[i] = Veff[i];
        }
    }
    double SolveSchrodingerEquation(const vector<double>& Enu)
    {
        for (int l = 0; l < Psi.size(); l++)
        {
            rhs_MT[0] = 0;
            for (int i = 1; i < Rmesh.size(); i++)
            {
                double Veff = Veff0[i] + 0.5 * l * (l + 1) / sqr(Rmesh[i]);
                rhs_MT[i] = 2 * (Veff - Enu[l]);
            }
            double dh = Rmesh[1] - Rmesh[0];
            ur[0] = 0;
            ur[1] = startSol(Z, l, dh);
            Numerov(rhs_MT, ur.size(), dh, ur);
            for (int i = 0; i < Rmesh.size(); i++)
            {
                temp[i] = sqr(ur[i]);
            }
            double norm = 1. / sqrt(integrate4<double>(temp, dh, temp.size()));
            for (int i = 0; i < ur.size(); i++)
            {
                ur[i] *= norm;
            }
            Psi_l[l] = ur;
            for (int ir = 1; ir < Rmesh.size(); ir++)
            {
                Psi_l[l][ir] /= Rmesh[ir];
            }
            Psi_l[l][0] = extrapolate(Psi_l[l][1], Psi_l[l][2], Rmesh[0], Rmesh[1], Rmesh[2]);
            for (int i = 0; i < Rmesh.size(); i++)
            {
                inhom[i] = -2. * ur[i];
            }
            urp[0] = 0;
            urp[1] = startSol(Z, l, dh);
            NumerovGen(rhs_MT, inhom, urp.size(), dh, urp);
            for (int i = 0; i < ur.size(); i++)
            {
                temp[i] = ur[i] * urp[i];
            }
            double alpha = integrate4<double>(temp, dh, temp.size());
            for (int i = 0; i < ur.size(); i++)
            {
                urp[i] -= alpha * ur[i];
            }
            Psip_l[l] = urp;
            for (int ir = 1; ir < Rmesh.size(); ir++)
            {
                Psip_l[l][ir] /= Rmesh[ir];
            }
            Psip_l[l][0] = extrapolate(Psip_l[l][1], Psip_l[l][2], Rmesh[0], Rmesh[1], Rmesh[2]);
            for (int i = 0; i < urp.size(); i++)
            {
                temp[i] = pow(urp[i], 2);
            }
            PsipPsip[l] = integrate4<double>(temp, dh, temp.size());
            int N0 = ur.size() - 1;
            double RMuffinTin = Rmesh[N0];
            double v1 = rhs_MT[N0] * ur[N0];
            double v0 = rhs_MT[N0 - 1] * ur[N0 - 1];
            double w1 = rhs_MT[N0] * urp[N0] + inhom[N0];
            double w0 = rhs_MT[N0 - 1] * urp[N0 - 1] + inhom[N0 - 1];
            double dudr = (ur[N0] - ur[N0 - 1]) / dh + 0.125 * dh * (3 * v1 + v0);
            double dupdr = (urp[N0] - urp[N0 - 1]) / dh + 0.125 * dh * (3 * w1 + w0);
            dlogPsi[l] = RMuffinTin * dudr / ur[N0] - 1;
            dlogPsip[l] = RMuffinTin * dupdr / urp[N0] - 1;
            Psi[l] = ur[N0] / RMuffinTin;
            Psip[l] = urp[N0] / RMuffinTin;
        }
    }
};
