#include "rs.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>

#include "Eigen/Dense"
#include "Eigen/Sparse"

using namespace Eigen;
namespace rs1D
{

    void calculateJacobianAndRhs(SparseMatrix<double> &jac, VectorXd &rhs, const VectorXd &u,
                                 int n, int neq, double Ro, double zmax);
    void newtonIteration(const SparseMatrix<double> &jac, const VectorXd &rhs, VectorXd &u, bool verbose = false);
    void saveResult(VectorXd &u);
    void checkJacobian();

    VectorXd rs1D(int n, double Ro, double zmax, bool verbose)
    {
        {
            std::cout << "Rotst3 | Base flow 1D | based on Lingowood 1996" << std::endl;
            std::cout << "n: " << n << " Ro: " << Ro << " zmax: " << zmax << std::endl;
        }
        int neq = 3;
        SparseMatrix<double> jac(neq * n, neq * n);
        VectorXd rhs = VectorXd::Zero(n * neq);
        VectorXd u = VectorXd::Zero(n * neq);

        for (int i = 0; i < u.size() / 3; i++)
        {
            u(i * 3 + 1) = 1;
            u(i * 3 + 2) = 1;
        }

        // check whether jacobian is correct
        checkJacobian();

        int nnewt = 15;

        for (int i = 0; i < nnewt; i++)
        {
            if (verbose == true)
                std::cout << "== it: " << i << std::endl;
            calculateJacobianAndRhs(jac, rhs, u, n, neq, Ro, zmax);
            newtonIteration(jac, rhs, u, verbose);
        }

        if (rhs.norm() > 1e-10)
            throw std::invalid_argument("Base flow not converged.");

        saveResult(u);
        return u;
    }
    void newtonIteration(const SparseMatrix<double> &jac, const VectorXd &rhs, VectorXd &u, bool verbose)
    {
        VectorXd du;
        VectorXd mrhs = -rhs;
        SparseLU<SparseMatrix<double>> solver;
        solver.compute(jac);
        if (solver.info() != Success)
        {
            // decomposition failed
            std::cout << "LU failed" << std::endl;
            return;
        }
        du = solver.solve(mrhs);
        if (solver.info() != Success)
        {
            // solving failed
            std::cout << "solve failed" << std::endl;
            return;
        }
        u = u + du;
        if (verbose == true)
        {
            std::cout << "norm of prev rhs: " << rhs.norm() << std::endl;
            std::cout << "eqs solved with acc: " << (jac * du + rhs).norm() << std::endl;
        }
        // {
        //     SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver2;
        //     solver2.compute(jac);
        //     std::cout << "rank of jac: " << solver2.rank() << std::endl;
        // }
    }
    void calculateJacobianAndRhs(SparseMatrix<double> &jac, VectorXd &rhs, const VectorXd &u,
                                 int n, int neq, double Ro, double zmax)
    {
        // base flow
        // int n = 11;
        // int neq = 3;
        // VectorXd u = VectorXd::Zero(3 * n);
        // double zmax = 20;
        double hz = zmax / double(n - 1);

        // rhs = VectorXd::Zero(neq * n);

        std::vector<Triplet<double>> tripletList;

        // parms
        // double Ro = -1;
        double Co = 2.0 - Ro - Ro * Ro;

        // BC
        double fBottom = 0;
        double gBottom = 0;
        double hBottom = 0;

        double fTop = 0;
        double gTop = 1;

        // Boundaries | Bottom
        {
            int i = 0;
            int iif = i * neq;
            int iig = iif + 1;
            int iih = iig + 1;

            rhs(iif) = u(iif) - fBottom;
            rhs(iig) = u(iig) - gBottom;
            rhs(iih) = u(iih) - hBottom;

            tripletList.push_back(Triplet<double>(iif, iif, 1));
            tripletList.push_back(Triplet<double>(iig, iig, 1));
            tripletList.push_back(Triplet<double>(iih, iih, 1));
        }

        // main body
        {
            for (int i = 1; i < n - 1; i++)
            {
                int iif = i * neq;
                int iig = iif + 1;
                int iih = iig + 1;
                int iifzp = iif + neq;
                int iifzm = iif - neq;
                int iigzp = iig + neq;
                int iigzm = iig - neq;
                int iihzm = iih - neq;

                rhs(iif) = Ro * (u(iif) * u(iif) + u(iih) * (u(iifzp) - u(iifzm)) / 2.0 / hz - (u(iig) * u(iig) - 1)) - Co * (u(iig) - 1) - (u(iifzp) - 2 * u(iif) + u(iifzm)) / hz / hz;
                rhs(iig) = Ro * (2 * u(iif) * u(iig) + u(iih) * (u(iigzp) - u(iigzm)) / 2.0 / hz) + Co * u(iif) - (u(iigzp) - 2 * u(iig) + u(iigzm)) / hz / hz;
                rhs(iih) = 2 * (u(iifzm) + u(iif)) / 2.0 + (u(iih) - u(iihzm)) / hz;

                // first
                {
                    double value = Ro * (2 * u(iif)) - (-2) / hz / hz;
                    tripletList.push_back(Triplet<double>(iif, iif, value));
                }
                {
                    double value = Ro * (u(iih) * (1) / 2.0 / hz) - (1) / hz / hz;
                    tripletList.push_back(Triplet<double>(iif, iifzp, value));
                }
                {
                    double value = Ro * (u(iih) * (-1) / 2.0 / hz) - (1) / hz / hz;
                    tripletList.push_back(Triplet<double>(iif, iifzm, value));
                }
                {
                    double value = Ro * (-(2 * u(iig))) - Co * (1);
                    tripletList.push_back(Triplet<double>(iif, iig, value));
                }
                {
                    double value = Ro * (1 * (u(iifzp) - u(iifzm)) / 2.0 / hz);
                    tripletList.push_back(Triplet<double>(iif, iih, value));
                }

                // second
                {
                    double value = Ro * (2 * u(iif) * 1) - (-2) / hz / hz;
                    tripletList.push_back(Triplet<double>(iig, iig, value));
                }
                {
                    double value = Ro * (u(iih) * (1) / 2.0 / hz) - (1) / hz / hz;
                    tripletList.push_back(Triplet<double>(iig, iigzp, value));
                }
                {
                    double value = Ro * (u(iih) * (-1) / 2.0 / hz) - (1) / hz / hz;
                    tripletList.push_back(Triplet<double>(iig, iigzm, value));
                }
                {
                    double value = Ro * (2 * 1 * u(iig)) + Co * 1;
                    tripletList.push_back(Triplet<double>(iig, iif, value));
                }
                {
                    double value = Ro * (1 * (u(iigzp) - u(iigzm)) / 2.0 / hz);
                    tripletList.push_back(Triplet<double>(iig, iih, value));
                }

                // third
                tripletList.push_back(Triplet<double>(iih, iih, 1 / hz));
                tripletList.push_back(Triplet<double>(iih, iihzm, -1 / hz));
                tripletList.push_back(Triplet<double>(iih, iif, 1));
                tripletList.push_back(Triplet<double>(iih, iifzm, 1));
            }
        }

        // Boundaries | Top
        {
            int i = n - 1;
            int iif = i * neq;
            int iig = iif + 1;
            int iih = iig + 1;
            int iifzm = iif - neq;
            int iihzm = iih - neq;

            rhs(iif) = u(iif) - fTop;
            rhs(iig) = u(iig) - gTop;
            rhs(iih) = 2 * (u(iifzm) + u(iif)) / 2.0 + (u(iih) - u(iihzm)) / hz;

            tripletList.push_back(Triplet<double>(iif, iif, 1));
            tripletList.push_back(Triplet<double>(iig, iig, 1));

            // third
            tripletList.push_back(Triplet<double>(iih, iih, 1 / hz));
            tripletList.push_back(Triplet<double>(iih, iihzm, -1 / hz));
            tripletList.push_back(Triplet<double>(iih, iif, 1));
            tripletList.push_back(Triplet<double>(iih, iifzm, 1));
        }

        jac.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void checkJacobian()
    {
        double eps = 1e-6;

        int n = 4;
        int neq = 3;
        SparseMatrix<double> jac(neq * n, neq * n);
        VectorXd rhs = VectorXd::Zero(n * neq);
        VectorXd u = VectorXd::Ones(n * neq);

        double Ro = -1;
        double zmax = 21;

        calculateJacobianAndRhs(jac, rhs, u, n, neq, Ro, zmax);

        int N = n * neq;

        // std::cout<<jac;

        for (int j = 0; j < N; j++)
        {
            VectorXd u1 = u;
            calculateJacobianAndRhs(jac, rhs, u1, n, neq, Ro, zmax);
            VectorXd rhs1 = rhs;

            VectorXd u2 = u;
            u2(j) += eps;
            calculateJacobianAndRhs(jac, rhs, u2, n, neq, Ro, zmax);
            VectorXd rhs2 = rhs;

            VectorXd correctColumn = (rhs2 - rhs1) / eps;

            for (int i = 0; i < N; i++)
            {
                if (std::abs(jac.coeff(i, j) - correctColumn(i)) > eps * 3)
                {
                    std::cout << rhs1 << std::endl;
                    std::cout << std::endl;
                    std::cout << rhs2 << std::endl;
                    std::cout << std::endl;
                    std::cout << jac << std::endl;

                    std::cout << "i: " << i << " j: " << j << " jac: " << jac.coeff(i, j) << " corr: " << correctColumn(i) << std::endl;
                    throw std::invalid_argument("lol");
                    return;
                }
            }
        }
        std::cout << "Jacobian for base flow correct." << std::endl;
    }

    void saveResult(VectorXd &u)
    {

        std::cout << "Saving to: " << std::filesystem::current_path() << "data.txt" << std::endl;

        std::ofstream myfile("data.txt");

        for (int i = 0; i < u.size() / 3; i++)
        {
            int ii = i * 3;
            myfile << u(ii) << "\t" << u(ii + 1) << "\t" << u(ii + 2) << "\n";
        }
        myfile.close();
    }

} // namespace