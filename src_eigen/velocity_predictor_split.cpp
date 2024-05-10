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
#define resl1                                                                  \
  (-RHSu(i, j) + B.W(i, j) * us(i - 1, j) + B.E(i, j) * us(i + 1, j) +         \
   B.S(i, j) * us(i, j - 1) + B.N(i, j) * us(i, j + 1) - B.P(i, j) * us(i, j))
#define resl2                                                                  \
  (-RHSv(i, j) + B.W(i, j) * vs(i - 1, j) + B.E(i, j) * vs(i + 1, j) +         \
   B.S(i, j) * vs(i, j - 1) + B.N(i, j) * vs(i, j + 1) - B.P(i, j) * vs(i, j))
#define utilde                                                                 \
  (-RHSu(i, j) + B.W(i, j) * us(i - 1, j) + B.E(i, j) * us(i + 1, j) +         \
   B.S(i, j) * us(i, j - 1) + B.N(i, j) * us(i, j + 1))
#define vtilde                                                                 \
  (-RHSv(i, j) + B.W(i, j) * vs(i - 1, j) + B.E(i, j) * vs(i + 1, j) +         \
   B.S(i, j) * vs(i, j - 1) + B.N(i, j) * vs(i, j + 1))
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
  for (kmom = 0; kmom <= adkmax; kmom = kmom + velkiner) {
    for (kgarb = 0; kgarb < velkiner; kgarb++) {
#pragma omp parallel for private(j, i) shared(us, blank, RHSu, B) collapse(2)
      for (j = j_s; j <= j_e; j++) {
        for (i = i_s; i <= i_e; i++) {
          us(i, j) = (1. - blank.blank(i, j)) * (utilde * B.Pinv(i, j)) +
                     blank.blank(i, j);
        }
      }

#pragma omp parallel for private(j, i) shared(vs, blank, RHSv, B) collapse(2)
        for (j = j_s; j <= j_e; j++) {
          for (i = i_s; i <= i_e; i++) {

            vs(i, j) = (1. - blank.blank(i, j)) * (vtilde * B.Pinv(i, j)) +
                       blank.blank(i, j);
          }
        }
        VelBC(us, vs); // impose boundary conditions on us and vs
      }
      resp = 0.0;
#pragma omp parallel for private(j, i) shared(us, blank, RHSu, B)        \
    collapse(2) reduction(+ : resp)
      for (j = j_s; j <= j_e; j++) {
        for (i = i_s; i <= i_e; i++) {
          resp = resp + (1 - blank.blank(i, j)) * resl1 * resl1;
        }
      }
#pragma omp parallel for private(j, i)  shared(vs, blank, RHSv, B)   collapse(2) reduction(+ : resp)
      for (j = j_s; j <= j_e; j++) {
        for (i = i_s; i <= i_e; i++) {
          resp = resp + (1 - blank.blank(i, j)) * resl2 * resl2;
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