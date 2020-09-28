#pragma once
class LAPW
{
    int Ksize;
    int ksize;
    int Rsize;
    int lMax;
    double Vol;
    double RMuffinTin;
    const vector<dvector3>& Km;
    function2D<double> omegal;
    function2D<double> C1;
    function2D<double> C2_1;
    function2D<double> C2_2;
    function2D<double> Olap_I;
    function2D<double> Olap;
    function2D<double> Ham;
    function2D<double> temp0;
    function2D<double> temp1;
    vector<function2D<double> > C2l;
    vector<vector<vector<double> > > weigh0;
    vector<vector<vector<double> > > weigh1;
    vector<vector<vector<double> > > weigh2;
    vector<vector<double> > weighI;
public:
    LAPW(int Ksize_, int ksize_, int Rsize_, int lMax_, double Vol_, double RMuffinTin_, const vector<dvector3>& Km_) :
        Ksize(Ksize_), ksize(ksize_), Rsize(Rsize_), lMax(lMax_), Vol(Vol_), RMuffinTin(RMuffinTin_), Km(Km_),
        omegal(Ksize, lMax + 1), C1(Ksize, lMax + 1), C2l(lMax + 1), C2_1(Ksize, Ksize), C2_2(Ksize, Ksize),
        weigh0(ksize), weigh1(ksize), weigh2(ksize), weighI(ksize),
        Olap_I(Ksize, Ksize), Olap(Ksize, Ksize), Ham(Ksize, Ksize), temp0(Ksize, Ksize), temp1(Ksize, Ksize)
    {
        for (int l = 0; l <= lMax; l++)
        {
            C2l[l].resize(Ksize, Ksize);
        }
        for (int ik = 0; ik < ksize; ik++)
        {
            weighI[ik].resize(Ksize);
            weigh0[ik].resize(lMax + 1);
            weigh1[ik].resize(lMax + 1);
            weigh2[ik].resize(lMax + 1);
            for (int il = 0; il <= lMax; il++)
            {
                weigh0[ik][il].resize(Ksize);
                weigh1[ik][il].resize(Ksize);
                weigh2[ik][il].resize(Ksize);
            }
        }
    }
    void ComputeInterstitialOverlap()
    {
        for (int i = 0; i < Ksize; i++)
        {
            Olap_I(i, i) = 1 - 4 * M_PI * sqr(RMuffinTin) * RMuffinTin / (3. * Vol);
            for (int j = i + 1; j < Ksize; j++)
            {
                double KKl = (Km[i] - Km[j]).length();
                Olap_I(i, j) = -4 * M_PI * sqr(RMuffinTin) * bessel_j(1, KKl * RMuffinTin) / (KKl * Vol);
                Olap_I(j, i) = Olap_I(i, j);
            }
        }
    }
    void ComputeEigensystem(const dvector3& k, const PartialWave& wave, const vector<double>& Enu, double VKSi, function<double>& Energy)
    {
        for (int iK = 0; iK < Ksize; iK++)
        {
            for (int il = 0; il <= lMax; il++)
            {
                double Dl, jl, jlDl;
                dlog_bessel_j(il, (k + Km[iK]).length() * RMuffinTin, Dl, jl, jlDl);
                omegal(iK, il) = -wave.psi(il) / wave.psip(il) * (Dl - wave.dlog(il)) / (Dl - wave.dlogp(il));
                C1(iK, il) = sqrt(4 * M_PI * (2 * il + 1) / Vol) * (jlDl - jl * wave.dlogp(il)) / (wave.psi(il) * (-wave.dlogp(il) + wave.dlog(il)));
            }
        }
        // This part of the code needs most of the time. Should be very optimized
        for (int iK = 0; iK < Ksize; iK++)
        {
            for (int jK = 0; jK < Ksize; jK++)
            {
                dvector3 qi(k + Km[iK]);
                dvector3 qj(k + Km[jK]);
                double qi_len = qi.length();
                double qj_len = qj.length();
                double argv = (qi_len * qj_len == 0) ? 1. : qi * qj / (qi_len * qj_len);
                double olapMT = 0;
                double hamMT = 0;
                for (int il = 0; il <= lMax; il++)
                {
                    double tC2l = C1(iK, il) * C1(jK, il) * Legendre(il, argv);
                    double toop = (1. + omegal(iK, il) * omegal(jK, il) * wave.PP(il));
                    olapMT += tC2l * toop;
                    hamMT += tC2l * (0.5 * (omegal(iK, il) + omegal(jK, il)) + toop * Enu[il]);
                    C2l[il](iK, jK) = tC2l;
                }
                Olap(iK, jK) = olapMT + Olap_I(iK, jK);
                Ham(iK, jK) = (0.25 * (qi * qi + qj * qj) + VKSi) * Olap_I(iK, jK) + hamMT;
            }
        }
        Eigensystem(Ksize, Energy, Olap, Ham); // When many K-points are used, this lapack call takes most of the time
    }
    void ComputeWeightsForDensity(int ik, const dvector3& k)
    {
        for (int il = 0; il <= lMax; il++)
        {
            for (int iK = 0; iK < Ksize; iK++)
            {
                for (int jK = 0; jK < Ksize; jK++)
                {
                    C2_1(iK, jK) = C2l[il](iK, jK) * (omegal(iK, il) + omegal(jK, il));
                    C2_2(iK, jK) = C2l[il](iK, jK) * omegal(iK, il) * omegal(jK, il);
                }
            }
            temp0.product("N", "N", Ham, C2l[il]);
            temp1.product("N", "T", temp0, Ham);
            for (int p = 0; p < Ksize; p++)
            {
                weigh0[ik][il][p] = temp1(p, p);
            }
            temp0.product("N", "N", Ham, C2_1);
            temp1.product("N", "T", temp0, Ham);
            for (int p = 0; p < Ksize; p++)
            {
                weigh1[ik][il][p] = temp1(p, p);
            }
            temp0.product("N", "N", Ham, C2_2);
            temp1.product("N", "T", temp0, Ham);
            for (int p = 0; p < Ksize; p++)
            {
                weigh2[ik][il][p] = temp1(p, p);
            }
        }
        temp0.product("N", "N", Ham, Olap_I);
        temp1.product("N", "T", temp0, Ham);
        for (int p = 0; p < Ksize; p++)
        {
            weighI[ik][p] = temp1(p, p);
        }
    }
    void ComputeMTDensity(function1D<double>& nMTRho, const function2D<double>& Energy, double mu, const vector<double>& wk, const PartialWave& wave)
    {
        nMTRho = 0;
        for (int il = 0; il <= lMax; il++)
        {
            double w0 = 0;
            double w1 = 0;
            double w2 = 0;
            for (int ik = 0; ik < ksize; ik++)
            {
                double sum0 = 0;
                double sum1 = 0;
                double sum2 = 0;
                for (int p = 0; p < Ksize; p++)
                {
                    if (Energy(ik, p) <= mu)
                    {
                        sum0 += weigh0[ik][il][p];
                        sum1 += weigh1[ik][il][p];
                        sum2 += weigh2[ik][il][p];
                    }
                }
                w0 += sum0 * wk[ik];
                w1 += sum1 * wk[ik];
                w2 += sum2 * wk[ik];
            }
            for (int ir = 0; ir < Rsize; ir++)
            {
                nMTRho[ir] += (w0 * sqr(wave.psi(il, ir)) + w1 * wave.psi(il, ir) * wave.psip(il, ir) + w2 * sqr(wave.psip(il, ir))) / (4 * M_PI);
            }
        }
        for (int ir = 0; ir < Rsize; ir++)
        {
            nMTRho[ir] *= 2;
        }
    }
    double ComputeInterstitialCharge(const function2D<double>& Energy, double mu, const vector<double>& wk)
    {
        double sIntRho = 0;
        for (int ik = 0; ik < ksize; ik++)
        {
            for (int p = 0; p < Ksize; p++)
            {
                if (Energy(ik, p) <= mu)
                {
                    sIntRho += weighI[ik][p] * wk[ik];
                }
            }
        }
        sIntRho *= 2;
        return sIntRho;
    }
    void PrintBandStructure(int ik, const function<double>& Energy, ostream& out)
    {
        out << setw(10) << ik / (ksize - 1.) << " ";
        for (int iK = 0; iK < Ksize; iK++)
        {
            out << setw(12) << Energy[iK] << " ";
        }
        out << endl;
    }
};