#include "eigen3/Eigen/Dense"
// #include <string> //Wheres Wallace????make
#include "./headers/CLASSES.h"
#include "./headers/COEFFVARS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/IBMVARS.h"
#include "./headers/INPUT.h"

#include "omp.h"

#include <iostream>
// #include<iostream>
using namespace INPUTDATA;
using namespace FIELDVARS;
using namespace Eigen;
using namespace COEFFVARS;
using namespace CLASSES;
using namespace std;
using namespace IBMVARS;

void VelocityPredictor(Eigen::ArrayXXd &us, Eigen::ArrayXXd &vs,
                       Eigen::ArrayXXd &RHSu, Eigen::ArrayXXd &RHSv,
                       CLASSES::SolverCoeffs &B, int &kmom, double &errormom) {
  kmom = 1;
  errormom = 100000.0;
  if (printflag == 0) {
    outputout << "--------------------------------------" << "\n"
              << "\n"
              << "Velocity Projection/Prediction" << "\n"
              << "\n"
              << "Iteration \t Error" << "\n"
              << "--------------------------------------" << "\n";
  }
  int kgarb = 0;
  velkiner = 3;
  double resp = 0;
  double resl = 0;
  for (kmom = 0; kmom <= adkmax; kmom = kmom + velkiner) {
    for (kgarb = 0; kgarb < velkiner; kgarb++) {
#pragma omp parallel for private(j, i) collapse(2)
      for (j = j_s; j <= j_e; j++) {
        for (i = i_s; i <= i_e; i++) {
          us(i, j) = -RHSu(i, j) + B.W(i, j) * us(i - 1, j) +
                     B.E(i, j) * us(i + 1, j) + B.S(i, j) * us(i, j - 1) +
                     B.N(i, j) * us(i, j + 1);
          us(i, j) = (1. - blank.blank(i, j)) * (us(i, j) * B.Pinv(i, j)) + blank.blank(i, j);

          vs(i, j) = -RHSv(i, j) + B.W(i, j) * vs(i - 1, j) +
                     B.E(i, j) * vs(i + 1, j) + B.S(i, j) * vs(i, j - 1) +
                     B.N(i, j) * vs(i, j + 1);
          vs(i, j) = (1. - blank.blank(i, j)) * (vs(i, j) * B.Pinv(i, j)) + blank.blank(i, j);
        }
      }
      VelBC(us, vs); // impose boundary conditions on us and vs
    }
    resp = 0.0;
#pragma omp parallel for private(j, i) collapse(2) reduction(+ : resp)
    for (j = j_s; j <= j_e; j++) {
      for (i = i_s; i <= i_e; i++) {
        resl = -RHSu(i, j) + B.W(i, j) * us(i - 1, j) +
               B.E(i, j) * us(i + 1, j) + B.S(i, j) * us(i, j - 1) +
               B.N(i, j) * us(i, j + 1) - B.P(i, j) * us(i, j);
        resp = resp + (1-blank.blank(i, j))*resl * resl;
        resl = -RHSv(i, j) + B.W(i, j) * vs(i - 1, j) +
               B.E(i, j) * vs(i + 1, j) + B.S(i, j) * vs(i, j - 1) +
               B.N(i, j) * vs(i, j + 1) - B.P(i, j) * vs(i, j);
        resp = resp + (1-blank.blank(i, j))*resl * resl;
      }
    }
    errormom = 0.5 * (sqrt(resp)) / (Nx * Ny);
    if (printflag == 0)
      outputout << kmom << "\t" << errormom << "\n";

    if (errormom < eps_ad && kmom > 3)
      break;
  }
  if (printflag == 0)
    outputout << "--------------------------------------" << "\n";
}