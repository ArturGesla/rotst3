#include "bodewadtWithPressure.hpp"

#include <iostream>
#include <fstream>
//#include <filesystem>
#include <complex>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "unsupported/Eigen/SparseExtra"

// for ev search
#include "arcomp.h"
#include "lcmatrxeEigen.h"
//#include "lcmatrxf.h"
#include "arlnsmat.h"
#include "arlgcomp.h"
#include "lcompsol.h"

using namespace Eigen;

namespace bwp
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
    void bwpStab(double Re, double Ro, int n, double zmax, double sigmaR, double sigmaI,
                 double alphaR, double alphaI, double betaR, double betaI, int nev)
    {
        //double Ro = 1; // bodewadt

        {
            std::cout << "Rotst3 |Stability analysis 1D | based on Lingowood 1996" << std::endl;
            std::cout << "n: " << n << " Ro: " << Ro << " zmax: " << zmax << " Re: " << Re << std::endl;
        }
        int neq = 4; // conti, r,th,z mom
        SparseMatrix<std::complex<double>> jac(neq * n, neq * n);
        SparseMatrix<std::complex<double>> B(neq * n, neq * n);
        VectorXcd rhs = VectorXd::Zero(n * neq);
        VectorXd u = VectorXd::Ones(n * neq);

        VectorXd U(n * neq);
        U = bwp::bwpBaseFlow(n, zmax, false, Ro);
        // std::complex<double> alpha(1.0, 0.0);
        // std::complex<double> beta(0.0, 0.0);
        std::complex<double> alpha(alphaR, alphaI);
        std::complex<double> beta(betaR, betaI);

        // std::complex<double> alpha(0.34, 0.0776);
        // std::complex<double> beta(-0.1174, 0);

        // check whether jacobian is correct
        checkJacobian();

        calculateJacobianAndRhsAndB(u, jac, B, rhs, U, n, neq, Ro, zmax, Re, alpha, beta);

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

        //what is life
        jac=jac/std::complex<double>(0.0,1.0);
        Bmat=Bmat/std::complex<double>(0.0,1.0); 

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
        saveMarket(jac, "jac.dat");
        saveMarket(Bmat, "Bmat.dat");
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
                std::cout << "Mag of evc " << j << " :" << eigenvector.real().norm() << "\t " << eigenvector.imag().norm() << std::endl;
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
    void calculateJacobianAndRhsAndB(const VectorXd &u, SparseMatrix<std::complex<double>> &jac,
                                     SparseMatrix<std::complex<double>> &B, VectorXcd &rhs,
                                     const VectorXd &U, int n, int neq, double Ro, double zmax,
                                     double Re, std::complex<double> alpha, std::complex<double> beta)
    {
        // base flow
        // int n = 11;
        // int neq = 3;
        // VectorXd u = VectorXd::Zero(3 * n);
        // double zmax = 20;

        if (u.size() != U.size())
        {
            throw std::invalid_argument("U and u size don't match.");
        }
        double hz = zmax / double(n - 1);

        // rhs = VectorXd::Zero(neq * n);

        std::vector<Triplet<std::complex<double>>> tripletList;
        std::vector<Triplet<std::complex<double>>> tripletListB;

        // parms
        // double Ro = -1;
        double Co = 2.0 - Ro - Ro * Ro;

        std::complex<double> i(0, 1);
        double r = Re;

        // BC
        double uBc = 0;

        // Boundaries | Bottom
        {
            int i = 0;
            int iip = i * neq;     // p,conti
            int iif = i * neq + 1; // f,rmom
            int iig = i * neq + 2; // g,thmom
            int iih = i * neq + 3; // h,zmom

            rhs(iip) = u(iip) - uBc;
            rhs(iif) = u(iif) - uBc;
            rhs(iig) = u(iig) - uBc;
            rhs(iih) = u(iih) - uBc;

            tripletList.push_back(Triplet<std::complex<double>>(iip, iip, 1));
            tripletList.push_back(Triplet<std::complex<double>>(iif, iif, 1));
            tripletList.push_back(Triplet<std::complex<double>>(iig, iig, 1));
            tripletList.push_back(Triplet<std::complex<double>>(iih, iih, 1));
        }
        // main body
        {
            // ii - point
            for (int ii = 1; ii < n - 1; ii++)
            {

                int iiF = ii * neq; // from external data
                int iiG = iiF + 1;
                int iiH = iiG + 1;

                int iiFzp = iiF + neq;
                int iiFzm = iiF - neq;
                int iiGzp = iiG + neq;
                int iiGzm = iiG - neq;
                int iiHzp = iiH + neq;
                int iiHzm = iiH - neq;

                int iip = ii * neq;
                int iipzp = iip + neq;
                int iipzm = iip - neq;

                int iif = ii * neq + 1;
                int iifzp = iif + neq;
                int iifzm = iif - neq;

                int iig = ii * neq + 2;
                int iigzp = iig + neq;
                int iigzm = iig - neq;

                int iih = ii * neq + 3;
                int iihzp = iih + neq;
                int iihzm = iih - neq;

                rhs(iip) = (u(iif) + u(iifzm)) / 2.0 * i * alpha + (u(iif) + u(iifzm)) / 2.0 / r +
                           +(u(iig) + u(iigzm)) / 2.0 * i * beta + (u(iih) - u(iihzm)) / hz;

                // rhs(iip) = (u(iif)) / 1.0 * i * alpha + (u(iif)) / 1.0 / r +
                //            +(u(iig)) / 1.0 * i * beta + (u(iihzp) - u(iihzm)) / 2.0 / hz;

                rhs(iif) = r * U(iiF) * u(iif) * i * alpha + u(iif) * U(iiF) + U(iiG) * u(iif) * i * r * beta +
                           +U(iiH) * (u(iifzp) - u(iifzm)) / 2.0 / hz + r * u(iih) * (U(iiFzp) - U(iiFzm)) / 2.0 / hz +
                           -2.0 * U(iiG) * u(iig) +
                           -(-i * alpha * (u(iip) + u(iipzp)) / 2.0 - u(iif) / r / r + i * u(iif) * alpha / r +
                             -alpha * alpha * u(iif) - beta * beta * u(iif) +
                             +(u(iifzp) - 2 * u(iif) + u(iifzm)) / hz / hz - 2.0 / r * u(iig) * i * beta);

                rhs(iig) = r * U(iiF) * u(iig) * i * alpha + u(iif) * U(iiG) + U(iiG) * u(iig) * i * r * beta +
                           +U(iiH) * (u(iigzp) - u(iigzm)) / 2.0 / hz + r * u(iih) * (U(iiGzp) - U(iiGzm)) / 2.0 / hz +
                           +U(iiF) * u(iig) + u(iif) * U(iiG) +
                           -(-i * beta * (u(iip) + u(iipzp)) / 2.0 - u(iig) / r / r + i * u(iig) * alpha / r +
                             -alpha * alpha * u(iig) - beta * beta * u(iig) +
                             +(u(iigzp) - 2 * u(iig) + u(iigzm)) / hz / hz + 2.0 / r * u(iif) * i * beta);

                rhs(iih) = r * U(iiF) * u(iih) * i * alpha + U(iiG) * u(iih) * i * r * beta +
                           +U(iiH) * (u(iihzp) - u(iihzm)) / 2.0 / hz + u(iih) * (U(iiGzp) - U(iiGzm)) / 2.0 / hz +
                           -(-(u(iipzp) - u(iip)) / hz + i * u(iih) * alpha / r +
                             -alpha * alpha * u(iih) - beta * beta * u(iih) +
                             +(u(iihzp) - 2 * u(iih) + u(iihzm)) / hz / hz);

                // first
                {
                    std::complex<double> value;

                    value = (1) / 2.0 * i * alpha + (1) / 2.0 / r;
                    tripletList.push_back(Triplet<std::complex<double>>(iip, iif, value));

                    value = (1) / 2.0 * i * alpha + (1) / 2.0 / r;
                    tripletList.push_back(Triplet<std::complex<double>>(iip, iifzm, value));

                    value = (1) / 2.0 * i * beta;
                    tripletList.push_back(Triplet<std::complex<double>>(iip, iig, value));

                    value = (1) / 2.0 * i * beta;
                    tripletList.push_back(Triplet<std::complex<double>>(iip, iigzm, value));

                    value = (1) / hz;
                    tripletList.push_back(Triplet<std::complex<double>>(iip, iih, value));

                    value = (-1) / hz;
                    tripletList.push_back(Triplet<std::complex<double>>(iip, iihzm, value));
                }
                // {
                //     std::complex<double> value;

                //     value = (1) / 1.0 * i * alpha + (1) / 1.0 / r;
                //     tripletList.push_back(Triplet<std::complex<double>>(iip, iif, value));

                //     value = (1) / 1.0 * i * beta;
                //     tripletList.push_back(Triplet<std::complex<double>>(iip, iig, value));

                //     value = (1) / 2.0 / hz;
                //     tripletList.push_back(Triplet<std::complex<double>>(iip, iihzp, value));

                //     value = (-1) / 2.0 / hz;
                //     tripletList.push_back(Triplet<std::complex<double>>(iip, iihzm, value));
                // }
                // second
                {
                    std::complex<double> value;

                    value = -(-i * alpha * 1.0 / 2.0);
                    tripletList.push_back(Triplet<std::complex<double>>(iif, iip, value));

                    value = -(-i * alpha * 1.0 / 2.0);
                    tripletList.push_back(Triplet<std::complex<double>>(iif, iipzp, value));

                    value = r * U(iiF) * 1.0 * i * alpha + 1.0 * U(iiF) + U(iiG) * 1.0 * i * r * beta +
                            -(-1 / r / r + i * 1.0 * alpha / r +
                              -alpha * alpha * 1.0 - beta * beta * 1.0 +
                              +(-2) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iif, iif, value));

                    value = -2.0 * U(iiG) * 1.0 +
                            -(-2.0 / r * 1.0 * i * beta);
                    tripletList.push_back(Triplet<std::complex<double>>(iif, iig, value));

                    value = r * 1.0 * (U(iiFzp) - U(iiFzm)) / 2.0 / hz;
                    tripletList.push_back(Triplet<std::complex<double>>(iif, iih, value));

                    value = U(iiH) * (1) / 2.0 / hz - ((1) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iif, iifzp, value));

                    value = U(iiH) * (-1) / 2.0 / hz - ((1) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iif, iifzm, value));

                    tripletListB.push_back(Triplet<std::complex<double>>(iif, iif, i*r));
                }

                // third
                {
                    std::complex<double> value;

                    value = -(-i * beta * 1.0 / 2.0);
                    tripletList.push_back(Triplet<std::complex<double>>(iig, iip, value));

                    value = -(-i * beta * 1.0 / 2.0);
                    tripletList.push_back(Triplet<std::complex<double>>(iig, iipzp, value));

                    value = 1.0 * U(iiG) +
                            1.0 * U(iiG) +
                            -(2.0 / r * 1.0 * i * beta);
                    tripletList.push_back(Triplet<std::complex<double>>(iig, iif, value));

                    value = r * U(iiF) * 1.0 * i * alpha + U(iiG) * 1.0 * i * r * beta +
                            +U(iiF) * 1.0 +
                            -(-1 / r / r + i * 1.0 * alpha / r +
                              -alpha * alpha * 1.0 - beta * beta * 1.0 +
                              +(-2) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iig, iig, value));

                    value = r * 1.0 * (U(iiGzp) - U(iiGzm)) / 2.0 / hz;
                    tripletList.push_back(Triplet<std::complex<double>>(iig, iih, value));

                    value = U(iiH) * (1) / 2.0 / hz - ((1) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iig, iigzp, value));

                    value = U(iiH) * (-1) / 2.0 / hz - ((1) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iig, iigzm, value));

                    tripletListB.push_back(Triplet<std::complex<double>>(iig, iig, i*r));
                }
                // fourth
                {
                    std::complex<double> value;

                    value = -(-(1) / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iih, iipzp, value));

                    value = -(-(-1) / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iih, iip, value));

                    value = r * U(iiF) * 1.0 * i * alpha + U(iiG) * 1.0 * i * r * beta +
                            +1.0 * (U(iiGzp) - U(iiGzm)) / 2.0 / hz +
                            -(i * 1.0 * alpha / r +
                              -alpha * alpha * 1.0 - beta * beta * 1.0 +
                              +(-2) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iih, iih, value));

                    value = U(iiH) * (1) / 2.0 / hz - (+(1) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iih, iihzp, value));

                    value = U(iiH) * (-1) / 2.0 / hz - (+(1) / hz / hz);
                    tripletList.push_back(Triplet<std::complex<double>>(iih, iihzm, value));

                    tripletListB.push_back(Triplet<std::complex<double>>(iih, iih, i*r));
                }
            }
        }

        // Boundaries | Top
        {

            int ii = n - 1;
            int iip = ii * neq;
            int iipzp = iip + neq;
            int iipzm = iip - neq;

            int iif = ii * neq + 1;
            int iifzp = iif + neq;
            int iifzm = iif - neq;

            int iig = ii * neq + 2;
            int iigzp = iig + neq;
            int iigzm = iig - neq;

            int iih = ii * neq + 3;
            int iihzp = iih + neq;
            int iihzm = iih - neq;

            rhs(iip) = (u(iif) + u(iifzm)) / 2.0 * i * alpha + (u(iif) + u(iifzm)) / 2.0 / r +
                      +(u(iig) + u(iigzm)) / 2.0 * i * beta + (u(iih) - u(iihzm)) / hz;
             //rhs(iip)=u(iip)-uBc;
            rhs(iif) = u(iif) - uBc;
            rhs(iig) = u(iig) - uBc;
            rhs(iih) = u(iih) - uBc;

             //tripletList.push_back(Triplet<std::complex<double>>(iip, iip, 1));
            //   first
            {
                std::complex<double> value;

                value = (1) / 2.0 * i * alpha + (1) / 2.0 / r;
                tripletList.push_back(Triplet<std::complex<double>>(iip, iif, value));

                value = (1) / 2.0 * i * alpha + (1) / 2.0 / r;
                tripletList.push_back(Triplet<std::complex<double>>(iip, iifzm, value));

                value = (1) / 2.0 * i * beta;
                tripletList.push_back(Triplet<std::complex<double>>(iip, iig, value));

                value = (1) / 2.0 * i * beta;
                tripletList.push_back(Triplet<std::complex<double>>(iip, iigzm, value));

                value = (1) / hz;
                tripletList.push_back(Triplet<std::complex<double>>(iip, iih, value));

                value = (-1) / hz;
                tripletList.push_back(Triplet<std::complex<double>>(iip, iihzm, value));
            } //shouldn't this give not fullrank?
            tripletList.push_back(Triplet<std::complex<double>>(iif, iif, 1));
            tripletList.push_back(Triplet<std::complex<double>>(iig, iig, 1));
            tripletList.push_back(Triplet<std::complex<double>>(iih, iih, 1));
        }
        // boundaries addon
        // {
        //     int i = n - 1;
        //     int iip = i * neq;     // p,conti
        //     int iif = i * neq + 1; // f,rmom
        //     int iig = i * neq + 2; // g,thmom
        //     int iih = i * neq + 3; // h,zmom

        //     rhs(0) = u(iip) - uBc;

        //     tripletList.push_back(Triplet<std::complex<double>>(0, iip, 1));
        // }

        jac.setFromTriplets(tripletList.begin(), tripletList.end());
        B.setFromTriplets(tripletListB.begin(), tripletListB.end());
    }

    void checkJacobian()
    {
        double eps = 1e-8;

        int n = 4;
        int neq = 4;
        SparseMatrix<std::complex<double>> jac(neq * n, neq * n);
        SparseMatrix<std::complex<double>> B(neq * n, neq * n);

        VectorXcd rhs = VectorXd::Zero(n * neq);
        VectorXd phi = VectorXd::Ones(n * neq);

        double Ro = 1;
        double zmax = 21;

        //       VectorXd u(n * 3);
        // u = rs1D::rs1D(n, Ro, zmax);
        VectorXd u = VectorXd::Random(n * neq);

        std::complex<double> alpha(1, 1);
        std::complex<double> beta(1, 1);

        double Re = 10;
        calculateJacobianAndRhsAndB(phi, jac, B, rhs, u, n, neq, Ro, zmax, Re, alpha, beta);

        int N = n * neq;

        // std::cout<<jac;
        double maxDiff = 0;

        // std::cout<<B<<std::endl;

        for (int j = 0; j < N; j++)
        {
            VectorXd phi1 = phi;
            calculateJacobianAndRhsAndB(phi1, jac, B, rhs, u, n, neq, Ro, zmax, Re, alpha, beta);
            VectorXcd rhs1 = rhs;

            VectorXd phi2 = phi;
            phi2(j) += eps;
            calculateJacobianAndRhsAndB(phi2, jac, B, rhs, u, n, neq, Ro, zmax, Re, alpha, beta);
            VectorXcd rhs2 = rhs;

            VectorXcd correctColumn = (rhs2 - rhs1) / eps;
            for (int i = 0; i < N; i++)
            {
                double diff = std::abs(jac.coeff(i, j) - correctColumn(i));
                maxDiff = std::max(diff, maxDiff);
                if (diff > eps * 10)
                {
                    std::cout << rhs1 << std::endl;
                    std::cout << std::endl;
                    std::cout << rhs2 << std::endl;
                    std::cout << std::endl;
                    std::cout << jac << std::endl;

                    std::cout << "i: " << i << " j: " << j << " jac: " << jac.coeff(i, j) << " corr: " << correctColumn(i) << std::endl;
                    std::cout << "diff: " << diff << std::endl;
                    std::cout.precision(17);

                    std::cout << std::fixed << rhs1(i) << std::endl;
                    std::cout << std::fixed << rhs2(i) << std::endl;

                    std::cout << std::fixed << phi1(i) << std::endl;
                    std::cout << std::fixed << phi2(i) << std::endl;

                    std::cout << std::scientific << rhs1(i) - rhs2(i) << std::endl;
                    std::cout << std::fixed << (rhs1(i).real() - rhs2(i).real()) / eps << std::endl;
                    std::cout << std::fixed << 1e-12 / 1e-12 << std::endl;

                    throw std::invalid_argument("Jacobian totally not correct.");
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