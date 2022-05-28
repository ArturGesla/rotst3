#pragma once
#include "Eigen/Dense"
namespace rs1DStab
{
    void rs1DStab(double Re = 1.0, int n = 401, double Ro =1, double zmax = 30, double sigmaR=1, double sigmaI=0,
    double alphaR=1, double alphaI=0, double betaR=0, double betaI=0, int nev=10);
}
