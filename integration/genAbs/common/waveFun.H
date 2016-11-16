/*---------------------------------------------------------------------------*\
   IH-Cantabria 2015 (http://www.ihcantabria.com/en/)
   IHFOAM 2015 (http://ihfoam.ihcantabria.com/) 

   Author(s):  Javier Lopez Lara (jav.lopez@unican.es)
               Gabriel Barajas   (barajasg@unican.es)
\*---------------------------------------------------------------------------*/

#ifndef waveFun_H
#define waveFun_H

namespace otherFun
{
    double interpolation (double x1, double x2, double y1, double y2, double xInt);
}

namespace StokesIFun
{
    double waveLength (double h, double T);
    double eta (double H, double Kx, double x, double Ky, double y, double omega, double t, double phase);
    double U (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z);
    double W (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z);
}

namespace StokesIIFun
{
    double eta (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase);
    double U (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z);
    double W (double H, double h, double Kx, double x, double Ky, double y, double omega, double t, double phase, double z);
    double timeLag (double H, double h, double Kx, double x, double Ky, double y, double T, double phase);
}

namespace Elliptic
{
    void ellipticIntegralsKE (double m, double* K, double* E);
    void JacobiSnCnDn (double u, double m, double* Sn, double* Cn, double* Dn);
}

namespace cnoidalFun
{
    double eta (double H, double m, double kx, double ky, double T, double x, double y, double t);
    double etaCnoidal1D (double H, double m, double t, double T, double K, double E);
    double etaMeanSq (double H, double m, double T);
    double d1EtaDx (double H, double m, double uCnoidal, double L, double K, double E);
    double d2EtaDx (double H, double m, double uCnoidal, double L, double K, double E);
    double d3EtaDx (double H, double m, double uCnoidal, double L, double K, double E);
    int calculations (double H, double d, double T, double* mOut, double* LOut);
    double U (double H, double h, double m, double kx, double ky, double T, double x, double y, double t, double z);
    double W (double H, double h, double m, double kx, double ky, double T, double x, double y, double t, double z);
}

namespace stokesVFun
{
    double A11 (double h, double k);
    double A13 (double h, double k);
    double A15 (double h, double k);
    double A22 (double h, double k);
    double A24 (double h, double k);
    double A33 (double h, double k);
    double A35 (double h, double k);
    double A44 (double h, double k);
    double A55 (double h, double k);
    double B22 (double h, double k);
    double B24 (double h, double k);
    double B33 (double h, double k);
    double B33k (double h, double k);
    double B35 (double h, double k);
    double B35k (double h, double k);
    double B44 (double h, double k);
    double B55 (double h, double k);
    double B55k (double h, double k);
    double C1 (double h, double k);
    double C1k (double h, double k);
    double C2 (double h, double k);
    double C2k (double h, double k);
    double C3 (double h, double k);
    double C4 (double h, double k);
    int StokesVNR (double H, double d, double T, double* kOut, double* LambdaOut, double* f1Out, double* f2Out );
    int StokesExtendedSolver (double H, double d, double T, double* kOut, double* LambdaOut, double* LambdaErrOut );
    double eta (double h, double kx, double ky, double lambda, double T, double x, double y, double t, double phase);
    double U (double d, double kx, double ky, double lambda, double T, double x, double y, double t, double phase, double z);
    double W (double d, double kx, double ky, double lambda, double T, double x, double y, double t, double phase, double z);
}

namespace BoussinesqFun
{
    double eta (double H, double h, double x, double y, double theta, double t, double X0);
    double Deta1 (double H, double h, double x, double y, double theta, double t, double X0);
    double Deta2 (double H, double h, double x, double y, double theta, double t, double X0);
    double Deta3 (double H, double h, double x, double y, double theta, double t, double X0);
    double U (double H, double h, double x, double y, double theta, double t, double X0, double z);
    double W (double H, double h, double x, double y, double theta, double t, double X0, double z);
}

#endif
