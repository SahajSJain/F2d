#include <eigen3/Eigen/Dense>
// #include <string> //Wheres Wallace????make
#include "./headers/CLASSES.h"
#include "./headers/COEFFVARS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/FLOWVARS.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/INPUT.h"
#include "headers/INPUT.h"
#include "omp.h"

// #include <iostream>
// #include <iostream>
using namespace INPUTDATA;
using namespace FLOWVARS;

using namespace FIELDVARS;
using namespace Eigen;
using namespace COEFFVARS;
using namespace CLASSES;
using namespace std;
void PointSOR(Eigen::ArrayXXd &p, Eigen::ArrayXXd &f, CLASSES::SolverCoeffs &A,
              int &kSOR, double &error) {
  kSOR = 0;
  error = 1000000.0;
  // Res *= 0;
  int kgarb = 0;
  double Ds = 0;
  double pt = 0;
  double resp = 0;
  double resl = 0;
  if (printflag == 0)
    outputout << "--------------------------------------" << "\n";

  if (printflag == 0)
    outputout << "\n"
              << "SOR - Pressure Solver" << "\n";
  if (printflag == 0)
    outputout << "\n"
              << "Overrelaxation factor = " << overrelax << "\n";

  if (printflag == 0)
    outputout << "\n"
              << "Iteration \t Error" << "\n";
  if (printflag == 0)
    outputout << "--------------------------------------" << "\n";
  while (kSOR < kmax && error > eps_pres) // k is the number of total iterations
  {
    kSOR = kSOR + kiner;
    for (kgarb = 0; kgarb < kiner; kgarb++) {
#pragma omp parallel for private(j, i) collapse(2)
      for (j = j_s; j <= j_e; j++) {
        for (i = i_s; i <= i_e; i++) {
          Ds = A.W(i, j) * p(i - 1, j) + A.E(i, j) * p(i + 1, j) +
               A.S(i, j) * p(i, j - 1) + A.N(i, j) * p(i, j + 1);
          pt = (Ds - f(i, j)) * A.Pinv(i, j);
          p(i, j) = (1 - overrelax) * p(i, j) + overrelax * pt;
        }
      }
      PresBC(p);
    }
    // Res *= 0;
    resp = 0.0;

#pragma omp parallel for private(j, i) collapse(2) reduction(+ : resp)
    for (j = j_s; j <= j_e; j++) {
      for (i = i_s; i <= i_e; i++) {
        // Ds = A.W(i, j) * p(i - 1, j) + A.E(i, j) * p(i + 1, j) +
        //      A.S(i, j) * p(i, j - 1) + A.N(i, j) * p(i, j + 1);
        resl = p(i, j) * A.P(i, j) -
               (A.W(i, j) * p(i - 1, j) + A.E(i, j) * p(i + 1, j) +
                A.S(i, j) * p(i, j - 1) + A.N(i, j) * p(i, j + 1) - f(i, j));
        resp = resp + resl * resl;
      }
    }
    error = (sqrt(resp)) / (Nx * Ny);
    if (printflag == 0)
      outputout << kSOR << "\t" << error << "\n";
  }
  if (printflag == 0)
    outputout << "--------------------------------------" << "\n";
}
