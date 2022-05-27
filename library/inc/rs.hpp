#pragma once
#include "Eigen/Dense"
using namespace Eigen;

namespace rs1D
{
    VectorXd rs1D(int n = 101, double Ro = 0, double zmax = 30, bool verbose = false);
}
