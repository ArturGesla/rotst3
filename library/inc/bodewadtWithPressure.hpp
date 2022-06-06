#pragma once

#include "Eigen/Dense"
using namespace Eigen;

namespace bwp
{
    VectorXd bwpBaseFlow(int n = 101, double zmax = 30, bool verbose = false);
}
