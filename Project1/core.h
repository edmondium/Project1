#ifndef _CORE_
#define _CORE_
class RadialWave
{
public:
    int N;
    vector<double> R;
    vector<double> Solution;
    vector<double> rhs;
    vector<double> Veff;
    vector<double> Veff0;
public:
    RadialWave(const vector<double>& Rmesh) : N(Rmesh.size()), Solution(N), R(N), rhs(N), Veff(N), Veff0(N)
    {
        for (int i = 0; i < N; i++)
        {
            R[i] = Rmesh[N - 1 - i];
        }
        Solution[0] = R[0] * exp(-R[0]);
        Solution[1] = R[1] * exp(-R[1]);
    }
    double operator()(double E)
    {
        double h = R[1] - R[0];
        for (int i = 0; i < N - 1; i++)
        {
            rhs[i] = 2 * (Veff[i] - E);
        }
        rhs[R.size() - 1] = 0;
        Numerov(rhs, Solution.size() - 1, h, Solution);
        int last = Solution.size() - 1;
        Solution[last] = Solution[last - 1] * (2 + pow(h, 2) * rhs[last - 1]) - Solution[last - 2];
        return Solution[last];
    }
    void Density(vector<double>& rho)
    {
        rho.resize(Solution.size());
        int N = Solution.size();
        for (int i = 0; i < Solution.size(); i++)
        {
            rho[i] = pow(Solution[N - 1 - i], 2);
        }
        double norm = 1. / integrate4<double>(rho, R[0] - R[1], rho.size());
        for (int i = 1; i < rho.size(); i++)
        {
            rho[i] *= norm / (4 * M_PI * sqr(R[N - i - 1])); //short version
        }
        rho[0] = 2 * rho[1] - rho[2];
    }
    void SetVeff0(const vector<double>& Veff)
    {
        for (int i = 0; i < R.size(); i++)
        {
            Veff0[i] = Veff[N - 1 - i];
        }
    }
    void AddCentrifugal(int l)
    {
        for (int i = 0; i < R.size() - 1; i++)
        {
            Veff[i] = Veff0[i] + 0.5 * l * (l + 1) / sqr(R[i]);
        }
    }
    double Vcore(int i)
    { 
        return Veff0[N - 1 - i];
    }
};
void FindCoreStates(const vector<int>& core, int Z, double dEz, RadialWave& wave, int& Nc, double& Ec, function1D<double>& coreRho)
{
    static vector<double> drho(coreRho.size());
    for (int ir = 0; ir < coreRho.size(); ir++)
    {
        coreRho[ir] = 0;
    }
    Nc = 0;
    Ec = 0;
    for (int l = 0; l < core.size(); l++)
    {
        wave.AddCentrifugal(l);
        double x = -0.5 * pow(Z, 2) / sqr(l + 1) - 3.;
        double v0 = wave(x);
        double v1 = v0;
        int n = 0;
        int j = 0;
        while (n < core[l] && x < 10.)
        {
            x += dEz;
            v1 = wave(x);
            if (v0 * v1 < 0)
            {
                double Energy = zeroin(x - dEz, x, wave, 1e-10);
                int dN = 2 * (2 * l + 1);
                wave.Density(drho);
                for (int ir = 0; ir < coreRho.size(); ir++)
                {
                    coreRho[ir] += drho[ir] * 2 * (2 * l + 1);
                }
                Ec += 2 * (2 * l + 1) * Energy;
                Nc += dN;
                clog << "Found core state for n = " << n + l << ", l = " << l << " at " << Energy << endl;
                n++;
                v0 = v1;
            }
        }
    }
}
#endif //_CORE_