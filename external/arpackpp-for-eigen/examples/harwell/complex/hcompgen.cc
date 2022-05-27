/*
   ARPACK++ v1.2 2/20/2000
   c++ interface to ARPACK code.

   MODULE MCompGen.cc.
   Example program that illustrates how to solve a complex
   generalized eigenvalue problem using the ARluCompGenEig class.

   1) Problem description:

      In this example we try to solve A*x = B*x*lambda in
      regular or shift and invert mode, where A is read
      from a file in Harwell-Boeing format.

   2) Included header files:

      File             Contents
      -----------      -------------------------------------------
      arerror.h        The ArpackError class definition.
      arlnsmat.h       The ARluNonSymMatrix class definition.
      arlgcomp.h       The ARluCompGenEig class definition.
      lcompsol.h       The Solution function.

   3) ARPACK Authors:

      Richard Lehoucq
      Kristyn Maschhoff
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include <iostream>
#include <string>
#include "arerror.h"
#include "arlnsmat.h"
#include "arlgcomp.h"
#include "lcompsol.h"


void PrintHelp()
/*
  Prints information about hcompgen usage.
*/
{

  std::cout << "ARPACK++ version 1.2 feb 2000" << std::endl;
  std::cout << "hcompgen: a generalized complex eigenvalue problems solver" << std::endl;
  std::cout << "usage:    hcompgen [parameters] file1 file2" << std::endl;
  std::cout << "parameters:" << std::endl;
  std::cout << "      -n <number of desired eigenvalues>" << std::endl;
  std::cout << "      -c <number of Arnoldi vectors per iteration>" << std::endl;
  std::cout << "      -l <maximum number of iterations>" << std::endl;
  std::cout << "      -s <real part of the shift>" << std::endl;
  std::cout << "      -i <imaginary part of the shift>" << std::endl;
  std::cout << "      -t <stopping criterion>" << std::endl;
  std::cout << "      -u <LU pivot threshold>" << std::endl;
  std::cout << "      -o <column ordering for factorization>" << std::endl;
  std::cout << "      -w <desired portion of the spectrum. " << std::endl;
  std::cout << "          acceptable values: LM, SM, LR, SR, LI, SI>" << std::endl;
  std::cout << std::endl;

} // PrintHelp.


bool ReadParameters(int n, char* v[], int &nev, int &ncv, int &maxit,
                    int &order, bool &shift, double &sigmar,
                    double &sigmai, double &tol, double &thresh,
                    std::string &which, std::string &fileA, std::string &fileB)
/*
  Reads parameters from the command line.
*/
{

  int  i;
  bool ok;

  // Defining default values.

  nev     = 5;
  ncv     = 0;
  maxit   = 0;
  order   = 1;
  shift   = false;
  sigmar  = 0.0;
  sigmai  = 0.0;
  tol     = 0.0;
  thresh  = 0.1;
  which   = "LM";
  fileA   = " ";
  fileB   = " ";
  ok      = true;

  // Returning if the number of parameters is even.

  if ((n==1)||(!(n%2))) {
    ok = false;
  }
  else {

    // Reading parameters.

    i = 1;
    while ((i<(n-2)) && (ok)) {

      if ((v[i][0] != '-') || (strlen(v[i]) != 2)) {
        ok = false;
        i += 2;
      }
      else {
        switch (v[i++][1]) {
        case 'n':
          nev = atoi(v[i++]);
          break;
        case 's':
          sigmar = atof(v[i++]);
          shift  = true;
          break;
        case 'i':
          sigmai = atof(v[i++]);
          shift  = true;
          break;
        case 'w':
          which = v[i++];
          break;
        case 'c':
          ncv = atoi(v[i++]);
          break;
        case 't':
          tol = atof(v[i++]);
          break;
        case 'l':
          maxit = atoi(v[i++]);
          break;
        case 'o':
          order = atoi(v[i++]);
          break;
        case 'u':
          thresh = atof(v[i++]);
          break;
        default :
          std::cout << "unrecognized parameter: -" << v[i-1][1] << std::endl;
          ok = false; 
        }
      }
    }
  }

  if (ok) {
    fileA = v[i++];
    fileB = v[i];
  }
  else {
    PrintHelp();
  }

  return ok;

} // ReadParameters.


int main(int argc, char* argv[])
{

  // Defining variables.

  int    nev;
  int    ncv;
  int    maxit;
  int    order;
  bool   shift;
  double sigmar;
  double sigmai;
  double tol;
  double thresh;
  std::string  which;
  std::string  fileA;
  std::string  fileB;

  // Reading parameters.

  if (!ReadParameters(argc, argv, nev, ncv, maxit, order, shift, sigmar,
                      sigmai, tol, thresh, which, fileA, fileB)) {
    return 0;
  }

  // Reading and storing matrices A and B.

  ARluNonSymMatrix<arcomplex<double>, double> A(fileA, thresh, order);
  ARluNonSymMatrix<arcomplex<double>, double> B(fileB, thresh, order);

  // Defining the eigenvalue problem.

  ARluCompGenEig<double> dprob(nev, A, B, which, ncv, tol, maxit);

  // Defining the shift.

  if (shift) dprob.SetShiftInvertMode(arcomplex<double>(sigmar,sigmai));

  // Finding eigenvalues and eigenvectors.

  try {
    dprob.FindEigenvectors();
  }
  catch (ArpackError) { return 0; }

  // Printing solution.

  Solution(A, B, dprob);

} // main.

