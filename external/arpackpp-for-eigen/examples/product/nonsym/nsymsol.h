/*
   ARPACK++ v1.2 2/18/2000
   c++ interface to ARPACK code.

   MODULE NSymSol.h
   Template functions that exemplify how to print information 
   about nonsymmetric standard  eigenvalue problems.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NSYMSOL_H
#define NSYMSOL_H

#include <cmath>
#include "blas1c.h"
#include "lapackc.h"
#include "matprod.h"
#include "arsnsym.h"

template<class ARMATRIX, class ARFLOAT>
void Solution(ARMATRIX &A, ARNonSymStdEig<ARFLOAT, ARMATRIX> &Prob)
/*
  Prints eigenvalues and eigenvectors of nonsymmetric eigen-problems
  on standard "std::cout" stream.
*/

{

  int     i, n, nconv, mode;
  ARFLOAT *Ax;
  ARFLOAT *ResNorm;

  n     = Prob.GetN();
  nconv = Prob.ConvergedEigenvalues();
  mode  = Prob.GetMode();

  std::cout << std::endl << std::endl << "Testing ARPACK++ class ARNonSymStdEig \n";
  std::cout << "Real nonsymmetric eigenvalue problem: A*x - lambda*x" << std::endl;
  switch (mode) {
  case 1:
    std::cout << "Regular mode" << std::endl << std::endl;
    break;
  case 3:
    std::cout << "Shift and invert mode" << std::endl << std::endl;
  }

  std::cout << "Dimension of the system            : " << n             << std::endl;
  std::cout << "Number of 'requested' eigenvalues  : " << Prob.GetNev() << std::endl;
  std::cout << "Number of 'converged' eigenvalues  : " << nconv         << std::endl;
  std::cout << "Number of Arnoldi vectors generated: " << Prob.GetNcv() << std::endl;
  std::cout << "Number of iterations taken         : " << Prob.GetIter() << std::endl;
  std::cout << std::endl;

  if (Prob.EigenvaluesFound()) {

    // Printing eigenvalues.

    std::cout << "Eigenvalues:" << std::endl;
    for (i=0; i<nconv; i++) {
      std::cout << "  lambda[" << (i+1) << "]: " << Prob.EigenvalueReal(i);
      if (Prob.EigenvalueImag(i)>=0.0) {
        std::cout << " + " << Prob.EigenvalueImag(i) << " I" << std::endl;
      }
      else {
        std::cout << " - " << fabs(Prob.EigenvalueImag(i)) << " I" << std::endl;
      }
    }
    std::cout << std::endl;
  }

  if (Prob.EigenvectorsFound()) {

    // Printing the residual norm || A*x - lambda*x ||
    // for the nconv accurately computed eigenvectors.

    Ax      = new ARFLOAT[n];
    ResNorm = new ARFLOAT[nconv+1];

    for (i=0; i<nconv; i++) {

      if (Prob.EigenvalueImag(i)==0.0) { // Eigenvalue is real.

        A.MultMv(Prob.RawEigenvector(i), Ax);
        axpy(n, -Prob.EigenvalueReal(i), Prob.RawEigenvector(i), 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1)/fabs(Prob.EigenvalueReal(i));

      }
      else {                 // Eigenvalue is complex.

        A.MultMv(Prob.RawEigenvector(i), Ax);
        axpy(n, -Prob.EigenvalueReal(i), Prob.RawEigenvector(i), 1, Ax, 1);
        axpy(n, Prob.EigenvalueImag(i), Prob.RawEigenvector(i+1), 1, Ax, 1);
        ResNorm[i] = nrm2(n, Ax, 1);
        A.MultMv(Prob.RawEigenvector(i+1), Ax);
        axpy(n, -Prob.EigenvalueImag(i), Prob.RawEigenvector(i), 1, Ax, 1);
        axpy(n, -Prob.EigenvalueReal(i), Prob.RawEigenvector(i+1), 1, Ax, 1);
        ResNorm[i] = lapy2(ResNorm[i], nrm2(n, Ax, 1))/
                     lapy2(Prob.EigenvalueReal(i), Prob.EigenvalueImag(i));
        ResNorm[i+1] = ResNorm[i];
        i++;

      }
    }

    for (i=0; i<nconv; i++) {
      std::cout << "||A*x(" << (i+1) << ") - lambda(" << (i+1);
      std::cout << ")*x(" << (i+1) << ")||: " << ResNorm[i] << "\n";
    }
    std::cout << "\n";

    delete[] Ax;
    delete[] ResNorm;

  }

} // Solution


#endif // NSYMSOL_H
