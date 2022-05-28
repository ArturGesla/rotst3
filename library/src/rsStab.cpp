#include "rsStab.hpp"
#include "rs.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>
#include <complex>

#include "Eigen/Dense"
#include "Eigen/Sparse"

// for ev search
#include "arcomp.h"
#include "lcmatrxeEigen.h"
//#include "lcmatrxf.h"
#include "arlnsmat.h"
#include "arlgcomp.h"
#include "lcompsol.h"

using namespace Eigen;

namespace rs1DStab
{

    void calculateJacobianAndRhsAndB(const VectorXd &phi, SparseMatrix<std::complex<double>> &jac,
                                     SparseMatrix<std::complex<double>> &B, VectorXcd &rhs,
                                     const VectorXd &u, int n, int neq, double Ro, double zmax,
                                     double Re, std::complex<double> alpha, std::complex<double> beta);
    // void newtonIteration(const SparseMatrix<double> &jac, const VectorXd &rhs, VectorXd &u);
    // void saveResult(VectorXd &u);
    void checkJacobian();
    void findev(SparseMatrix<std::complex<double>> jac, SparseMatrix<std::complex<double>> B,
                int nev, double sigmaR, double sigmaI);

    // void rs1DStab(double Re, int n, double Ro, double zmax)
    void rs1DStab(double Re, int n, double Ro, double zmax, double sigmaR, double sigmaI,
                  double alphaR, double alphaI, double betaR, double betaI, int nev)
    {
        {
            std::cout << "Rotst3 |Stability analysis 1D | based on Lingowood 1996" << std::endl;
            std::cout << "n: " << n << " Ro: " << Ro << " zmax: " << zmax << " Re: " << Re << std::endl;
        }
        int neq = 6;
        SparseMatrix<std::complex<double>> jac(neq * n, neq * n);
        SparseMatrix<std::complex<double>> B(neq * n, neq * n);
        VectorXcd rhs = VectorXd::Zero(n * neq);
        VectorXd phi = VectorXd::Zero(n * neq);

        VectorXd u(n * 3);
        u = rs1D::rs1D(n, Ro, zmax, false);
        // std::complex<double> alpha(1.0, 0.0);
        // std::complex<double> beta(0.0, 0.0);
        std::complex<double> alpha(alphaR, alphaI);
        std::complex<double> beta(betaR, betaI);

        // std::complex<double> alpha(0.34, 0.0776);
        //         std::complex<double> beta(-0.1174, 0);

        calculateJacobianAndRhsAndB(phi, jac, B, rhs, u, n, neq, Ro, zmax, Re, alpha, beta);

        // check whether jacobian is correct
        checkJacobian();

        // J*phi=omega*phi*B
        findev(jac, B, nev, sigmaR, sigmaI);
    }
    void findev(SparseMatrix<std::complex<double>> jac, SparseMatrix<std::complex<double>> Bmat,
                int nev, double sigmaR, double sigmaI)
    {

        // Defining variables;

        int n;                          // Dimension of the problem.
        int nnza, nnzb;                 // Number of nonzero elements in A and B.
        int *irowa, *irowb;             // pointers to arrays that store the row
                                        // indices of the nonzeros in A and B.
        int *pcola, *pcolb;             // pointers to arrays of pointers to the
                                        // beginning of each column of A and B in
                                        // valA and ValB.
        arcomplex<double> rho;          // parameter used in CompMatrixE.
        arcomplex<double> *valA, *valB; // pointers to arrays that store the
                                        // nonzero elements of A and B.

        // Creating complex matrices A and B.

        n = jac.rows();
        rho = arcomplex<double>(10.0, 0.0);
        CompMatrixE(n, rho, nnza, valA, irowa, pcola, jac);
        ARluNonSymMatrix<arcomplex<double>, double> A(n, nnza, valA, irowa, pcola);

        CompMatrixE(n, rho, nnzb, valB, irowb, pcolb, Bmat);
        ARluNonSymMatrix<arcomplex<double>, double> B(n, nnzb, valB, irowb, pcolb);

        {
            SparseQR<SparseMatrix<std::complex<double>>, COLAMDOrdering<int>> solver2;
            solver2.compute(jac);
            std::cout << "rank of jac: " << solver2.rank() << " n: " << jac.rows() << std::endl;
        }
        // std::cout<<Bmat<<std::endl;

        // for (int i = 0; i < Bmat.nonZeros(); i++)
        // {
        //     std::cout<<*(valB+i)<<std::endl;
        // }

        // Defining what we need: the four eigenvectors nearest to sigma.

        ARluCompGenEig<double> dprob(nev, A, B, arcomplex<double>(sigmaR, sigmaI));

        // Finding eigenvalues and eigenvectors.

        dprob.FindEigenvectors();

        // Printing solution.

        Solution(A, B, dprob);

        

        // save to file
        int nconv = dprob.ConvergedEigenvalues();

        {
            // Printing eigenvalues.
            std::ofstream myfile("evs.dat", std::ios_base::app);
            // std::cout << "Eigenvalues:" << std::endl;
            for (int i = 0; i < nconv; i++)
            {
                // std::cout << "  lambda[" << (i+1) << "]: " << dprob.Eigenvalue(i) << std::endl;
                myfile << dprob.Eigenvalue(i).real() << "\t" << dprob.Eigenvalue(i).imag() << "\n";
            }
            // std::cout << std::endl;
            myfile.close();
        }

        {
            // Printing vectors.
            std::ofstream myfile("evcs.dat");
            // std::cout << "Eigenvalues:" << std::endl;
            for (int j = 0; j < nconv; j++)
            {
                // std::cout << "  lambda[" << (i+1) << "]: " << dprob.Eigenvalue(i) << std::endl;
                // myfile << dprob.RawEigenvector(i).real() << "\t" << dprob.RawEigenvector(i).imag() << "\n";
                VectorXcd eigenvector(n);
                for (size_t i = 0; i < n; i++)
                {
                    eigenvector(i) = dprob.Eigenvector(j, i);
                }
                myfile << eigenvector.real().transpose() << std::endl;
                std::cout<<"Mag of evc "<<j<<" :"<<eigenvector.real().norm()<<"\t "<<eigenvector.imag().norm()<<std::endl;
            }
            // std::cout << std::endl;
            myfile.close();
        }
    }
    // void newtonIteration(const SparseMatrix<double> &jac, const VectorXd &rhs, VectorXd &u)
    // {
    //     VectorXd du;
    //     VectorXd mrhs = -rhs;
    //     SparseLU<SparseMatrix<double>> solver;
    //     solver.compute(jac);
    //     if (solver.info() != Success)
    //     {
    //         // decomposition failed
    //         std::cout << "LU failed" << std::endl;
    //         return;
    //     }
    //     du = solver.solve(mrhs);
    //     if (solver.info() != Success)
    //     {
    //         // solving failed
    //         std::cout << "solve failed" << std::endl;
    //         return;
    //     }
    //     u = u + du;
    //     {
    //         std::cout << "norm of prev rhs: " << rhs.norm() << std::endl;
    //         std::cout << "eqs solved with acc: " << (jac * du + rhs).norm() << std::endl;
    //     }
    //     // {
    //     //     SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> solver2;
    //     //     solver2.compute(jac);
    //     //     std::cout << "rank of jac: " << solver2.rank() << std::endl;
    //     // }
    // }
    void calculateJacobianAndRhsAndB(const VectorXd &phi, SparseMatrix<std::complex<double>> &jac,
                                     SparseMatrix<std::complex<double>> &B, VectorXcd &rhs,
                                     const VectorXd &u, int n, int neq, double Ro, double zmax,
                                     double Re, std::complex<double> alpha, std::complex<double> beta)
    {
        // base flow
        // int n = 11;
        // int neq = 3;
        // VectorXd u = VectorXd::Zero(3 * n);
        // double zmax = 20;

        if (u.size() * 2 != phi.size())
        {
            throw std::invalid_argument("U and phi size don't match.");
        }
        double hz = zmax / double(n - 1);

        // rhs = VectorXd::Zero(neq * n);

        std::vector<Triplet<std::complex<double>>> tripletList;
        std::vector<Triplet<std::complex<double>>> tripletListB;

        // parms
        // double Ro = -1;
        double Co = 2.0 - Ro - Ro * Ro;

        std::complex<double> i(0, 1);
        std::complex<double> gammaSq = alpha * alpha + beta * beta;
        std::complex<double> alphaBar = alpha - i * Re / Ro;
        std::complex<double> gammaBarSq = alpha * alphaBar + beta * beta;
        std::complex<double> omega(0, 0);

        // BC
        double phiBC = 0;

        // Boundaries | Bottom
        {
            int i = 0;
            int phi1 = i * neq;
            int phi2 = i * neq + 1;
            int phi3 = i * neq + 2;
            int phi4 = i * neq + 3;
            int phi5 = i * neq + 4;
            int phi6 = i * neq + 5;

            rhs(phi1) = phi(phi1) - phiBC;
            rhs(phi2) = phi(phi2) - phiBC;
            rhs(phi3) = phi(phi3) - phiBC;
            rhs(phi4) = phi(phi4) - phiBC;
            rhs(phi5) = phi(phi5) - phiBC;
            rhs(phi6) = phi(phi6) - phiBC;

            tripletList.push_back(Triplet<std::complex<double>>(phi1, phi1, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi2, phi2, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi3, phi3, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi4, phi4, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi5, phi5, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi6, phi6, 1));
        }

        // main body
        {
            for (int iip = 1; iip < n - 1; iip++)
            {

                int iif = iip * 3; // neq for u is 3
                int iig = iif + 1;
                int iih = iig + 1;

                int iifzp = iif + 3;
                int iifzm = iif - 3;
                int iigzp = iig + 3;
                int iigzm = iig - 3;
                int iihzp = iih + 3;
                int iihzm = iih - 3;

                int phi1 = iip * neq;
                int phi2 = iip * neq + 1;
                int phi3 = iip * neq + 2;
                int phi4 = iip * neq + 3;
                int phi5 = iip * neq + 4;
                int phi6 = iip * neq + 5;

                double dUdz = (u(iifzp) - u(iifzm)) / 2.0 / hz;
                double dVdz = (u(iigzp) - u(iigzm)) / 2.0 / hz;
                double dWdz = (u(iihzp) - u(iihzm)) / 2.0 / hz;

                rhs(phi1) = -(phi(phi1 + neq) - phi(phi1 - neq)) / 2.0 / hz +
                            phi(phi2);
                rhs(phi2) = -(phi(phi2 + neq) - phi(phi2 - neq)) / 2.0 / hz +
                            (gammaSq + i * Re * (alpha * u(iif) + beta * u(iig) - omega) + Ro * u(iif)) * phi(phi1) +
                            Ro * u(iih) * phi(phi2) +
                            (alphaBar * dUdz + beta * dVdz) * Re * phi(phi3) +
                            i * gammaBarSq * Re * phi(phi4) +
                            -(2 * Ro * u(iig) + Co) * phi(phi5);
                rhs(phi3) = -(phi(phi3 + neq) - phi(phi3 - neq)) / 2.0 / hz +
                            -i * phi(phi1);
                rhs(phi4) = -(phi(phi4 + neq) - phi(phi4 - neq)) / 2.0 / hz +
                            -i * Ro * u(iih) * phi(phi1) / Re +
                            -i * phi(phi2) / Re +
                            -(gammaSq + i * Re * (alpha * u(iif) + beta * u(iig) - omega) + Ro * dWdz) * phi(phi3) / Re;
                rhs(phi5) = -(phi(phi5 + neq) - phi(phi5 - neq)) / 2.0 / hz +
                            phi(phi6);
                rhs(phi6) = -(phi(phi6 + neq) - phi(phi6 - neq)) / 2.0 / hz +
                            (2 * Ro * u(iig) + Co) * phi(phi1) +
                            (alphaBar * dVdz - beta * dUdz) * Re * phi(phi3) +
                            beta * Ro * phi(phi4) +
                            (gammaSq + i * Re * (alpha * u(iif) + beta * u(iig) - omega) + Ro * u(iif)) * phi(phi5) +
                            Ro * u(iih) * phi(phi6);

                // first
                {
                    tripletList.push_back(Triplet<std::complex<double>>(phi1, phi1 + neq, -1 / 2.0 / hz));
                    tripletList.push_back(Triplet<std::complex<double>>(phi1, phi1 - neq, 1 / 2.0 / hz));
                    std::complex<double> value = 1;
                    tripletList.push_back(Triplet<std::complex<double>>(phi1, phi2, value));
                }
                // second
                {
                    tripletList.push_back(Triplet<std::complex<double>>(phi2, phi2 + neq, -1 / 2.0 / hz));
                    tripletList.push_back(Triplet<std::complex<double>>(phi2, phi2 - neq, 1 / 2.0 / hz));

                    std::complex<double> value;

                    value = (gammaSq + i * Re * (alpha * u(iif) + beta * u(iig) - omega) + Ro * u(iif));
                    tripletList.push_back(Triplet<std::complex<double>>(phi2, phi1, value));
                    value = Ro * u(iih);
                    tripletList.push_back(Triplet<std::complex<double>>(phi2, phi2, value));
                    value = (alphaBar * dUdz + beta * dVdz) * Re;
                    tripletList.push_back(Triplet<std::complex<double>>(phi2, phi3, value));
                    value = i * gammaBarSq * Re;
                    tripletList.push_back(Triplet<std::complex<double>>(phi2, phi4, value));
                    value = -(2 * Ro * u(iig) + Co);
                    tripletList.push_back(Triplet<std::complex<double>>(phi2, phi5, value));

                    value = i * Re;
                    tripletListB.push_back(Triplet<std::complex<double>>(phi2, phi1, value));
                }
                // third
                {
                    tripletList.push_back(Triplet<std::complex<double>>(phi3, phi3 + neq, -1 / 2.0 / hz));
                    tripletList.push_back(Triplet<std::complex<double>>(phi3, phi3 - neq, 1 / 2.0 / hz));

                    std::complex<double> value;

                    value = -i;
                    tripletList.push_back(Triplet<std::complex<double>>(phi3, phi1, value));
                }
                // fourth
                {
                    tripletList.push_back(Triplet<std::complex<double>>(phi4, phi4 + neq, -1 / 2.0 / hz));
                    tripletList.push_back(Triplet<std::complex<double>>(phi4, phi4 - neq, 1 / 2.0 / hz));

                    std::complex<double> value;

                    value = -i * Ro * u(iih) * 1.0 / Re;
                    tripletList.push_back(Triplet<std::complex<double>>(phi4, phi1, value));
                    value = -i * 1.0 / Re;
                    tripletList.push_back(Triplet<std::complex<double>>(phi4, phi2, value));
                    value = -(gammaSq + i * Re * (alpha * u(iif) + beta * u(iig) - omega) + Ro * dWdz) * 1.0 / Re;
                    tripletList.push_back(Triplet<std::complex<double>>(phi4, phi3, value));

                    value = -i;
                    tripletListB.push_back(Triplet<std::complex<double>>(phi4, phi3, value));
                }
                // fifth
                {
                    tripletList.push_back(Triplet<std::complex<double>>(phi5, phi5 + neq, -1 / 2.0 / hz));
                    tripletList.push_back(Triplet<std::complex<double>>(phi5, phi5 - neq, 1 / 2.0 / hz));

                    std::complex<double> value;

                    value = 1.0;
                    tripletList.push_back(Triplet<std::complex<double>>(phi5, phi6, value));
                }
                // sixth
                {
                    tripletList.push_back(Triplet<std::complex<double>>(phi6, phi6 + neq, -1 / 2.0 / hz));
                    tripletList.push_back(Triplet<std::complex<double>>(phi6, phi6 - neq, 1 / 2.0 / hz));

                    std::complex<double> value;

                    value = (2 * Ro * u(iig) + Co);
                    tripletList.push_back(Triplet<std::complex<double>>(phi6, phi1, value));
                    value = (alphaBar * dVdz - beta * dUdz) * Re;
                    tripletList.push_back(Triplet<std::complex<double>>(phi6, phi3, value));
                    value = beta * Ro;
                    tripletList.push_back(Triplet<std::complex<double>>(phi6, phi4, value));
                    value = (gammaSq + i * Re * (alpha * u(iif) + beta * u(iig) - omega) + Ro * u(iif));
                    tripletList.push_back(Triplet<std::complex<double>>(phi6, phi5, value));
                    value = Ro * u(iih);
                    tripletList.push_back(Triplet<std::complex<double>>(phi6, phi6, value));

                    value = i * Re;
                    tripletListB.push_back(Triplet<std::complex<double>>(phi6, phi5, value));
                }
            }
        }

        // Boundaries | Top
        {

            int i = n - 1;
            int phi1 = i * neq;
            int phi2 = i * neq + 1;
            int phi3 = i * neq + 2;
            int phi4 = i * neq + 3;
            int phi5 = i * neq + 4;
            int phi6 = i * neq + 5;

            rhs(phi1) = phi(phi1) - phiBC;
            rhs(phi2) = phi(phi2) - phiBC;
            rhs(phi3) = phi(phi3) - phiBC;
            rhs(phi4) = phi(phi4) - phiBC;
            rhs(phi5) = phi(phi5) - phiBC;
            rhs(phi6) = phi(phi6) - phiBC;

            tripletList.push_back(Triplet<std::complex<double>>(phi1, phi1, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi2, phi2, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi3, phi3, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi4, phi4, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi5, phi5, 1));
            tripletList.push_back(Triplet<std::complex<double>>(phi6, phi6, 1));
        }

        jac.setFromTriplets(tripletList.begin(), tripletList.end());
        B.setFromTriplets(tripletListB.begin(), tripletListB.end());
    }

    void checkJacobian()
    {
        double eps = 1e-6;

        int n = 4;
        int neq = 6;
        SparseMatrix<std::complex<double>> jac(neq * n, neq * n);
        SparseMatrix<std::complex<double>> B(neq * n, neq * n);

        VectorXcd rhs = VectorXd::Zero(n * neq);
        VectorXd phi = VectorXd::Ones(n * neq);

        double Ro = 1;
        double zmax = 21;

        //       VectorXd u(n * 3);
        // u = rs1D::rs1D(n, Ro, zmax);
        VectorXd u = VectorXd::Ones(n * 3);

        std::complex<double> alpha(1, 1);
        std::complex<double> beta(1, 1);

        double Re = 10;

        calculateJacobianAndRhsAndB(phi, jac, B, rhs, u, n, 6, Ro, zmax, Re, alpha, beta);

        int N = n * neq;

        // std::cout<<jac;
        double maxDiff = 0;

        for (int j = 0; j < N; j++)
        {
            VectorXd phi1 = phi;
            calculateJacobianAndRhsAndB(phi1, jac, B, rhs, u, n, 6, Ro, zmax, Re, alpha, beta);
            VectorXcd rhs1 = rhs;

            VectorXd phi2 = phi;
            phi2(j) += eps;
            calculateJacobianAndRhsAndB(phi2, jac, B, rhs, u, n, 6, Ro, zmax, Re, alpha, beta);
            VectorXcd rhs2 = rhs;

            VectorXcd correctColumn = (rhs2 - rhs1) / eps;
            for (int i = 0; i < N; i++)
            {
                double diff = std::abs(jac.coeff(i, j) - correctColumn(i));
                maxDiff = std::max(diff, maxDiff);
                if (diff > eps * 3)
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
        std::cout << "Jacobian for stability analysis correct. Max diff: " << maxDiff << std::endl;
    }

    // void saveResult(VectorXd &u)
    // {

    //     std::cout << "Saving to: " << std::filesystem::current_path() << "data.txt" << std::endl;

    //     std::ofstream myfile("data.txt");

    //     for (int i = 0; i < u.size() / 3; i++)
    //     {
    //         int ii = i * 3;
    //         myfile << u(ii) << "\t" << u(ii + 1) << "\t" << u(ii + 2) << "\n";
    //     }
    //     myfile.close();
    // }

} // namespace