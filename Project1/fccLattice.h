#pragma once
class FccLattice
{
    double LatConst;
    double Volume;
    dvector3 a0;
    dvector3 a1;
    dvector3 a2;
    dvector3 b0;
    dvector3 b1;
    dvector3 b2;
    dvector3 GammaPoint;
    dvector3 LPoint;
    dvector3 KPoint;
    dvector3 XPoint;
    dvector3 WPoint;
    vector<dvector3> Kmesh;
    vector<dvector3> kmesh;
    vector<double> wkp;
public:
    double Vol() const
    { 
        return Volume;
    }
    int Ksize() const
    { 
        return Kmesh.size();
    }
    int ksize() const
    { 
        return kmesh.size();
    }
    double wk(int ik) const
    { 
        return wkp[ik];
    }
    const vector<double>& wk() const
    { 
        return wkp;
    }
    const vector<dvector3>& Km()
    { 
        return Kmesh;
    }
    const dvector3& K(int i) const
    { 
        return Kmesh[i];
    }
    const dvector3& k(int i) const
    { 
        return kmesh[i];
    }
    FccLattice(double LatConst_) : LatConst(LatConst_)
    {
        a0 = dvector3(0.5 * LatConst, 0.5 * LatConst, 0);
        a1 = dvector3(0.5 * LatConst, 0, 0.5 * LatConst);
        a2 = dvector3(0, 0.5 * LatConst, 0.5 * LatConst);
        Volume = fabs(VProduct(a0, a1, a2));
        clog << "Volume is " << Volume << endl;
        b0 = (2 * M_PI / Volume) * cross(a1, a0);
        b1 = (2 * M_PI / Volume) * cross(a0, a2);
        b2 = (2 * M_PI / Volume) * cross(a2, a1);
        double brs = 2 * M_PI / LatConst; //wavenumber
        GammaPoint = dvector3(0, 0, 0);
        LPoint = dvector3(0.5 * brs, 0.5 * brs, 0.5 * brs);
        KPoint = dvector3(0.75 * brs, 0.75 * brs, 0);
        XPoint = dvector3(1 * brs, 0, 0);
        WPoint = dvector3(1 * brs, 0.5 * brs, 0); //symmetry points
    }
    void GenerateReciprocalVectors(int q, double CutOffK)
    {
        list<dvector3> Kmesh0;
        for (int n = -q; n < q; n++)
        {
            for (int l = -q; l < q; l++)
            {
                for (int m = -q; m < q; m++)
                {
                    Kmesh0.push_back(n * b0 + l * b1 + m * b2);
                }
            }
        }
        Kmesh0.sort(cmp);
        int Ksize = 0;
        for (list<dvector3>::const_iterator l = Kmesh0.begin(); l != Kmesh0.end(); l++, Ksize++)
        {
            if (l->length() > CutOffK)
            {
                break;
            }
        }
        Kmesh.resize(Ksize);
        int j = 0;
        for (list<dvector3>::const_iterator l = Kmesh0.begin(); l != Kmesh0.end() && j < Ksize; l++, j++)
        {
            Kmesh[j] = *l;
        }
        clog << "K-mesh size = " << Kmesh.size() << endl;
    }
    void ChoosePointsInFBZ(int nkp, int type = 0)
    {
        if (type == 0)
        {
            list<dvector3> kp;
            GenerateAllKPoints(nkp, b0, b1, b2, kp);
            clog << "Number of all k-points = " << kp.size() << endl;
            ChooseIrreducible(kp, kmesh, wkp);
            clog << "Number of irreducible k-points = " << kmesh.size() << endl;
        }
        else
        {
            int tm = 4 * static_cast<int>(pow(nkp, 3) / 4.);
            int nkp = tm;
            clog << "Number of k-points = " << nkp << endl;
            kmesh.resize(nkp);
            int N0 = kmesh.size() / 4;
            for (int i = 0; i < N0; i++)
            {
                kmesh[i] = GammaPoint + (XPoint - GammaPoint) * i / (N0 - 1.);
            }
            for (int i = 0; i < N0; i++)
            {
                kmesh[N0 + i] = XPoint + (LPoint - XPoint) * i / (N0 - 1.);
            }
            for (int i = 0; i < N0; i++)
            {
                kmesh[N0 * 2 + i] = LPoint + (GammaPoint - LPoint) * i / (N0 - 1.);
            }
            for (int i = 0; i < N0; i++)
            {
                kmesh[N0 * 3 + i] = GammaPoint + (KPoint - GammaPoint) * i / (N0 - 1.);
            }
        }
    }
    double RMuffinTin() const
    { 
        return 0.5 * a0.length();
    }
};