#pragma once
#include "Eigen/Dense"
namespace bwp
{
    // void bwpStab(double Re = 21.7, double Ro=1.0, int n = 801,  double zmax = 30, double sigmaR=-0.218, double sigmaI=0,
    // double alphaR=0.34, double alphaI=0.0776, double betaR=-0.1174, double betaI=0, int nev=10);
    
    // void bwpStab(double Re = 47.4, double Ro=0.8, int n = 1601,  double zmax = 40, double sigmaR=-0.218, double sigmaI=0,
    // double alphaR=0.406, double alphaI=0.141, double betaR=0.0495, double betaI=0, int nev=40);
    
    
    void bwpStab(double Re = 47.4, double Ro=0.4, int n = 1601,  double zmax = 40, double sigmaR=-0.218, double sigmaI=0,
    double alphaR=0.403, double alphaI=0.191, double betaR=0.157, double betaI=0, int nev=40);
    
    // void bwpStab(double Re = 27.4, int n = 801,  double zmax = 30, double sigmaR=-0.218, double sigmaI=0,
    // double alphaR=0.487, double alphaI=0.0, double betaR=0.115, double betaI=0, int nev=10);
}
