#pragma once
#include "Eigen/Dense"
namespace bwp
{
    void bwpStab(double Re = 1.7, int n = 801,  double zmax = 30, double sigmaR=-0.218, double sigmaI=0,
    double alphaR=0.34, double alphaI=0.0776, double betaR=-0.1174, double betaI=0, int nev=10);
}
