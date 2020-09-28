//#include <cstdio>
/*int main()
{
    printf("hello from %s!\n", "Project1");
    return 0;
}*/
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <list>
#include <algorithm>
#include "numerov.h"
#include "integrate.h"
#include "kmesh.h"
#include "bessel.h"
#include "zero.h"
#include "erf.h"
#include "function.h"
#include "blas.h"
#include "exchcorr.h"
#include "core.h"
#include "util.h"
#include "fccLattice.h"
#include "PartialWave.h"
#include "ChemicalPotential.h"
#include "LAPW.h"
#include "miscellaneous.h"
using namespace std;
int main(int argc, char* argv[])
{
    int Z = 29;                     // Number of electrons in the atom
    double LatConst = 6.8219117;    // Lattic constant
    int lMax = 5;                   // Maximum l considered in calculation
    int N = 1001;                   // Number of points in radial mesh
    int nkp = 8;                    // Number of k-points in first Brillouin zone: (nkp x nkp x nkp)
    double CutOffK = 3.5;           // Largest length of reciprocal vectors K (only shorter vectors are taken into account)
    double dEz = 0.1;               // Step in serching for core states
    double mu_min = 0.0;            // Interval where chemical potential is looking for
    double mu_max = 1.5;
    double Mix = 0.2;               // Linear mixing parameter for charge
    double Vmix = 0.2;              // Linear mixing parameter for potential
    double precision = 1e-5;        // accuracy of energy
    int nitt = 2000;                 // maximum number of iterations
    bool read = false;              // weather to read input potential
    /* Core states */
    int lMaxCore = 2;
    vector<int> core(lMaxCore + 1);
    core[0] = 3; // 1s-3s is in core
    core[1] = 2; // 1p, 2p is in core
    core[2] = 0; // no d in core
    int i = 0;
    while (++i < argc)
    {
        string str(argv[i]);
        if (str == "-Z" && i < argc - 1)
        {
            Z = atoi(argv[++i]);
        }
        if (str == "-dE" && i < argc - 1)
        {
            dEz = atof(argv[++i]);
        }
        if (str == "-lmax" && i < argc - 1)
        {
            lMax = atoi(argv[++i]);
        }
        if (str == "-N" && i < argc - 1)
        {
            N = atoi(argv[++i]);
        }
        if (str == "-nkp" && i < argc - 1)
        {
            nkp = atoi(argv[++i]);
        }
        if (str == "-CutOffK" && i < argc - 1)
        {
            CutOffK = atof(argv[++i]);
        }
        if (str == "-mumin" && i < argc - 1)
        {
            mu_min = atof(argv[++i]);
        }
        if (str == "-mumax" && i < argc - 1)
        {
            mu_max = atof(argv[++i]);
        }
        if (str == "-Mix" && i < argc - 1)
        {
            Mix = atof(argv[++i]);
        }
        if (str == "-Vmix" && i < argc - 1)
        {
            Vmix = atof(argv[++i]);
        }
        if (str == "-precision" && i < argc - 1)
        {
            precision = atof(argv[++i]);
        }
        if (str == "-nitt" && i < argc - 1)
        {
            nitt = atoi(argv[++i]);
        }
        if (str == "-read")
        {
            read = true;
        }
        if (str == "-h" || str == "--help")
        {
            clog << "**************** LAPW program for fcc lattice ********" << endl;
            clog << "**      Developed more by Edmond Febrinicko Armay   **" << endl;
            clog << "**      Adopted from Kristjan Haule, 21.11.2005     **" << endl;
            clog << "******************************************************" << endl;
            clog << endl;
            clog << "lapw [-dE double] [] []" << endl;
            clog << "Options:   -Z          Number of electrons (" << Z << ")" << endl;
            clog << "           -dE         Step in searching for states (" << dEz << ")" << endl;
            clog << "           -lmax       Maximum l used in calculation (" << lMax << ")" << endl;
            clog << "           -N          Number of points in radial mesh (" << N << ")" << endl;
            clog << "           -nkp        Number of k-points from irreducible Brillouin zone used (nkp x nkp x nkp) (" << nkp << ")" << endl;
            clog << "           -CutOffK    Largest length of reciprocal vectors K (only shorter vectors are taken into account) (" << CutOffK << ")" << endl;
            clog << "           -mumin      Start looking for chemical potential (" << mu_min << ")" << endl;
            clog << "           -mumax      Stop looking for chemical potential (" << mu_max << ")" << endl;
            clog << "           -Mix        Linear mixing parameter for charge (" << Mix << ")" << endl;
            clog << "           -Vmix       Linear mixing parameter for potential (" << Vmix << ")" << endl;
            clog << "           -precision  Total energy accuracy required (" << precision << ")" << endl;
            clog << "           -nitt       Maximum number of iterations (" << nitt << ")" << endl;
            clog << "           -read       Weather to read from file Potential_input.txt (" << read << ")" << endl;
            clog << "******************************************************" << endl;
            return 0;
        }
    }
    clog.precision(10);
    int Zcor = 0;
    for (int l = 0; l < core.size(); l++)
    {
        Zcor += 2 * (2 * l + 1) * core[l];
    }
    int Zval = Z - Zcor;
    clog << "Z core = " << Zcor << " and Z valence = " << Zval << endl;
    // Generates and stores momentum points
    FccLattice fcc(LatConst);                           // Information about lattice
    double RMuffinTin = fcc.RMuffinTin();               // Muffin-Tin radius chosen such that spheres touch
    double VMT = 4 * M_PI * pow(RMuffinTin, 3) / 3.;    // Volume of MT
    double Vinter = fcc.Vol() - VMT;                    // Volume of the interstitial region
    clog << "Muffin-Tin radius = " << RMuffinTin << endl;
    clog << "Volume of the MT sphere    = " << VMT << endl;
    clog << "Volume of the interstitial = " << Vinter << endl;
    // For solving Schrödinger equation
    PartialWave wave(N, RMuffinTin, Z, lMax);   // This is for partial waves in the valence band
    RadialWave wavec(wave.R());                 // This is for core states
    fcc.GenerateReciprocalVectors(4, CutOffK);  // Reciprocal bravais lattice is built, K points taken into account only for |K|<CutOff
    fcc.ChoosePointsInFBZ(nkp, 0);              // Choose the path in the first Brillouin zone or the k-points in the irreducible Brillouin zone
    //ExchangeCorrelation XC(3);                  // (look http://physics.nist.gov/PhysRefData/DFTdata/Tables/ptable.html)
    ExchangeCorrelation XC(5);
    vector<double> Enu(lMax + 1);               // Linearization energies.
    Enu[0] = 0.11682;                           // Most of high energy partial waves should be centered around mu
    Enu[1] = 0.18794;                           // In general, it is a good idea to change them through iteration
    Enu[2] = 0.211145;
    for (int il = 3; il < Enu.size(); il++)
    {
        Enu[il] = 0.3;
    }
    double VKSi = 0;                            // Potential in the interstitials
    double mu = 0;                              // Chemical potential
    vector<double> Uhartree(wave.Rsize()), Vxc(wave.Rsize()), Veff(wave.Rsize());               // Hartree and exchange-correlation and KS potential
    function1D<double> TotRho(wave.Rsize()), nTotRho(wave.Rsize()), drho(wave.Rsize());         // total density, input and output
    function2D<double> Energy(fcc.ksize(), fcc.Ksize());                                        // E(k,p)- bands
    LAPW lapw(fcc.Ksize(), fcc.ksize(), wave.Rsize(), lMax, fcc.Vol(), RMuffinTin, fcc.Km());   // basic class for LAPW calculation
    function1D<double>  MTRho(wave.Rsize()), coreRho(wave.Rsize());                             // partial densities, core and valence MT
    // Starting guess for the Hartree and exchange correlation potential
    for (int i = 1; i < wave.Rsize(); i++)
    {
        Veff[i] = -VeffP(wave.R(i)) / wave.R(i);
    }
    Veff[0] = extrapolate(Veff[1], Veff[2], wave.R(0), wave.R(1), wave.R(2));
    // Read potential from file
    if (read)
    {
        ifstream inp("Potential_input.txt");
        double r;
        double v;
        double v0;
        int ii = 0;
        while (inp >> r && inp >> v && inp >> v0 && ii < wave.Rsize())
        {
            Veff[ii++] = -v / r;
        }
        Veff[0] = extrapolate(Veff[1], Veff[2], wave.R(0), wave.R(1), wave.R(2));
    }
    double zeroMT = Veff[wave.Rsize() - 1]; // adjusting MT zero
    for (int i = 0; i < wave.Rsize(); i++)
    {
        Veff[i] -= zeroMT;
    }
    double pEc = 0; // previous core energy
    /////////////////// Starting calculation ////////////////////////////////
    lapw.ComputeInterstitialOverlap(); // Overlap in the interstitials can be calculated outside
    //////////////////// Main loop ///////////////////////////
    for (int itt = 0; itt < nitt; itt++)
    {
        clog << "****** Iteration number " << itt << " ************" << endl;
        /////////////////// Potential part /////////////////////////////
        if (itt > 0)
        {
            SolvePoisson(Z, wave.R(), TotRho, Uhartree); // Poisson gives Hartree potential
            for (int i = 0; i < wave.Rsize(); i++)
            {
                Vxc[i] = XC.Vxc(rs(TotRho[i])); // Adding exchange-correlation part
            }
            // This is total KS effective potential
            for (int i = 1; i < wave.Rsize(); i++)
            {
                Veff[i] = (1 - Vmix) * Veff[i] + Vmix * ((-Z + Uhartree[i]) / wave.R(i) + Vxc[i]);
            }
            Veff[0] = extrapolate(Veff[1], Veff[2], wave.R(0), wave.R(1), wave.R(2));
            // Zero energy is chosen at each iteration such that potential vanishes on MT sphere
            double VMTzero = Veff[wave.Rsize() - 1]; // New MT zero
            for (int i = 0; i < wave.Rsize(); i++)
            {
                Veff[i] -= VMTzero; // Shift potential to new MT zero
            }
        }
        ///////////// Schrödinger equation for MT region //////////////////
        wave.SetVeff0(Veff);
        wave.SolveSchrodingerEquation(Enu);// Energy Enu depends on l
        //////////////// Core part ////////////////////////////////////////
        wavec.SetVeff0(Veff);
        int Nc;
        double Ec;
        FindCoreStates(core, Z, dEz, wavec, Nc, Ec, coreRho);
        clog << "core Z = " << Nc << ", and energy = " << Ec << endl;
        ///////////// Main LAPW loop over k points /////////////////////////
        for (int ik = 0; ik < fcc.ksize(); ik++)
        {
            dvector3 k = fcc.k(ik); // current k-point
            lapw.ComputeEigensystem(k, wave, Enu, VKSi, Energy[ik]);
            lapw.ComputeWeightsForDensity(ik, k);
        }
        /////////////////// New chemical potential ////////////////////////////
        FChemicalPotential chemp(Zval, Energy, fcc.wk());
        mu = zeroin(mu_min, mu_max, chemp, 1e-13);
        clog << "Chemical potential found at " << COLOR(BLUE, mu) << endl;
        //////////////////// New valence density /////////////////////////////
        lapw.ComputeMTDensity(MTRho, Energy, mu, fcc.wk(), wave);
        double sIntRho = lapw.ComputeInterstitialCharge(Energy, mu, fcc.wk());
        double sMTRho = IntegrateCharge(wave.R(), MTRho);
        cout << "MTRho = " << sMTRho + sIntRho << endl;
        double scoreRho = IntegrateCharge(wave.R(), coreRho);
        ////////////////// New total charge /////////////////////////////////
        for (int i = 0; i < wave.Rsize(); i++)
        {
            nTotRho[i] = MTRho[i] + coreRho[i];
        }
        clog << "Weight in the MT sphere = " << sMTRho << " and in the interstitials = " << sIntRho << " and in core = " << scoreRho << endl;
        double renorm = Z / (sMTRho + sIntRho + scoreRho);
        clog << "Total charge found = " << scoreRho + sMTRho + sIntRho << " should be " << Z << " -> renormalizing charge by " << renorm << endl;
        ///////////////// Renormalization of charge ////////////////////////
        for (int i = 0; i < wave.Rsize(); i++)
        {
            nTotRho[i] *= renorm;
        }
        //////////////////// Charge difference //////////////////////////////
        for (int i = 0; i < wave.Rsize(); i++)
        {
            drho[i] = fabs(TotRho[i] - nTotRho[i]);
        }
        double ChargeDifference = integrate4<double>(drho, wave.R(1) - wave.R(0), drho.size());
        if (itt == 0)
        {
            Mix = 1;
        }
        //////////// Linear mixing. Could be improved with Broyden, Vanderbild or Johannson mixing /////////////
        for (int i = 0; i < wave.Rsize(); i++)
        {
            TotRho[i] = TotRho[i] * (1 - Mix) + nTotRho[i] * Mix;
        }
        ////////////// Convergence criteria ////////////////////////////////////////////////////
        clog << "Core energy difference = " << COLOR(RED, fabs(Ec - pEc)) << ", charge difference = " << COLOR(GREEN, ChargeDifference) << endl;
        if (fabs(Ec - pEc) < precision)
        {
            break;
        }
        pEc = Ec;
    }
    ////// Printing of band structure /////////////////////////
    ofstream bands("bands.txt");
    fcc.ChoosePointsInFBZ(nkp, 1);// Points for plotting band structure
    function1D<double> epsk(fcc.Ksize());
    if (nitt == 0)
    {
        wave.SetVeff0(Veff);
        wave.SolveSchrodingerEquation(Enu);
    }
    for (int ik = 0; ik < fcc.ksize(); ik++)
    {
        dvector3 k = fcc.k(ik); // current k-point
        lapw.ComputeEigensystem(k, wave, Enu, VKSi, epsk);
        lapw.PrintBandStructure(ik, epsk, bands);
    }
    //////////////////// Some printing //////////////////////////////
    ofstream out("charge1.txt");
    for (int ir = 0; ir < wave.Rsize(); ir++)
    {
        out << wave.R(ir) << "  " << TotRho[ir] * 4 * M_PI * sqr(wave.R(ir)) << "  " << MTRho[ir] * 4 * M_PI * sqr(wave.R(ir)) << "  " << coreRho[ir] * 4 * M_PI * sqr(wave.R(ir)) << endl;
    }
    ofstream outc("Charge2.txt");
    for (int ir = 0; ir < wave.Rsize(); ir++)
    {
        outc << wave.R(ir) << "  " << TotRho[ir] << "  " << MTRho[ir] << "  " << coreRho[ir] << endl;
    }
    ofstream outu("potential1.txt");
    for (int i = 0; i < wave.Rsize(); i++)
    {
        outu << wave.R(i) << " " << Veff[i] << " " << Uhartree[i] << " " << Vxc[i] << endl;
    }
    ofstream outp("Potential2.txt");
    for (int i = 0; i < wave.Rsize(); i++)
    {
        outp << wave.R(i) << " " << -Veff[i] * wave.R(i) << "  " << VeffP(wave.R(i)) << endl;
    }
    return 0;
}