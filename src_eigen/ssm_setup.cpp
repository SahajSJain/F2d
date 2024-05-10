// my_program.cpp
// libraries to include
#include "./headers/CLASSES.h"
#include "./headers/COEFFVARS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/IBMVARS.h"
#include "./headers/INPUT.h"
#include "omp.h"
#include <cmath>
#include "eigen3/Eigen/Dense"
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;
using namespace INPUTDATA;
using namespace CLASSES;
using namespace FIELDVARS;
using namespace COEFFVARS;
using namespace IBMVARS;
using namespace Eigen;

inline bool
checkinorout(double &xpoint,
             double &ypoint) // level based function that tells you wheter you
                             // are in or out. Is one on the inside
{
  bool inorout;
  double io_radius, io_radius_1;
  switch (ssm_bodytype) {
  case 2: // ellipse
    io_radius = hypot(xpoint - ssm_x1, ypoint - ssm_y1);
    io_radius_1 = hypot(xpoint - ssm_x2, ypoint - ssm_y2);
    inorout = ((io_radius + io_radius_1) < (2.0 * ssm_rad));
    break;
  case 3: // rectangle
    io_radius = abs(xpoint - ssm_x1);
    io_radius_1 = abs(ypoint - ssm_y1);
    inorout = ((io_radius < ssm_x2) && (io_radius_1 < ssm_y2));
    break;
  default: // circle
    io_radius = hypot(xpoint - ssm_x1, ypoint - ssm_y1);
    inorout = (io_radius < ssm_rad);
    break;
  }
  return inorout;
}

void ssm_setup() {

  blank.assign(Nx_c, Ny_c);  // initialize blank face
  SSMADC.assign(Nx_c, Ny_c); // Advection Diffusion Coefficients
  SSMDEq.assign(Nx_c, Ny_c); // Poisson Coefficients

  SSMADC = ADC;
  SSMDEq = DEq;
  n_u.assign(Nx_c, Ny_c); //
  n_v.assign(Nx_c, Ny_c); //
  n_p.assign(Nx_c, Ny_c); //
  //   n_u0.assign(Nx_c, Ny_c); // //these are not needed
  //   n_v0.assign(Nx_c, Ny_c); //
  //   n_p0.assign(Nx_c, Ny_c); //
  string ssm_bodytype_string = "";
  switch (ssm_bodytype) {
  case 2: // ellipse
    ssm_bodytype_string = "Ellipse";
    ssm_xc = 0.5 * (ssm_x1 + ssm_x2);
    ssm_yc = 0.5 * (ssm_y1 + ssm_y2);
    break;
  case 3: // rectangle
    ssm_bodytype_string = "Rectangle";
    ssm_xc = (ssm_x1);
    ssm_yc = (ssm_y1);
    break;
  default: // circle
    ssm_bodytype_string = "Circle";
    ssm_xc = (ssm_x1);
    ssm_yc = (ssm_y1);

    break;
  }
  outputout << "\n"
            << "SSM Body type:" << ssm_bodytype << "\n";
  outputout << "\n"
            << "SSM x1:" << ssm_x1 << "\n";
  outputout << "\n"
            << "SSM y1:" << ssm_y1 << "\n";
  outputout << "\n"
            << "SSM x2:" << ssm_x2 << "\n";
  outputout << "\n"
            << "SSM y2:" << ssm_y2 << "\n";
  outputout << "\n"
            << "SSM rad:" << ssm_rad << "\n";

  blank.blank *= 0; // default is no ssm body
  // check in which domain: 1 is outside fluid, 0 is inside solid
  if (ssmflag) { // create ssm body only if flag on
#pragma omp parallel for private(j, i) collapse(2)
    for (j = j_s; j <= j_e; j++) {
      for (i = i_s; i <= i_e; i++) {
        blank.blank(i, j) = checkinorout(xc(i), yc(j));
      }
    }
  }
  //   cout << "\n" << blank.blank << "\n";
  // check if in same domain and change ssmadc and ssmdeq
  int boundarynumber = 0;
#pragma omp parallel for private(j, i) collapse(2) reduction(+ : boundarynumber)
  for (j = j_s; j <= j_e; j++) {
    for (i = i_s; i <= i_e; i++) {
      blank.E(i, j) = (blank.blank(i, j) == blank.blank(i + 1, j));
      blank.W(i, j) = (blank.blank(i, j) == blank.blank(i - 1, j));
      blank.N(i, j) = (blank.blank(i, j) == blank.blank(i, j + 1));
      blank.S(i, j) = (blank.blank(i, j) == blank.blank(i, j - 1));
      blank.ibnode(i, j) = (!(blank.E(i, j) && blank.W(i, j) && blank.S(i, j) &&
                              blank.N(i, j))) &&
                           blank.blank(i, j);
      boundarynumber = boundarynumber + blank.ibnode(i, j);

      // modify SSMADC for dirichlet  bc
      // ap=bp + sum((1-id_nb)*b_nb)
      // a_nb=id_nb*b_nb
      SSMADC.E(i, j) = ADC.E(i, j) * blank.E(i, j);
      SSMADC.W(i, j) = ADC.W(i, j) * blank.W(i, j);
      SSMADC.N(i, j) = ADC.N(i, j) * blank.N(i, j);
      SSMADC.S(i, j) = ADC.S(i, j) * blank.S(i, j);
      SSMADC.P(i, j) = ADC.P(i, j) + ((1 - blank.E(i, j)) * ADC.E(i, j) +
                                      (1 - blank.W(i, j)) * ADC.W(i, j) +
                                      (1 - blank.N(i, j)) * ADC.N(i, j) +
                                      (1 - blank.S(i, j)) * ADC.S(i, j));
      SSMADC.Pinv(i, j) = 1 / SSMADC.P(i, j);

      // modify SSMADC for dirichlet  bc
      // ap=bp + sum((1-id_nb)*b_nb)
      // a_nb=id_nb*b_nb
      SSMADC.E(i, j) = ADC.E(i, j) * blank.E(i, j);
      SSMADC.W(i, j) = ADC.W(i, j) * blank.W(i, j);
      SSMADC.N(i, j) = ADC.N(i, j) * blank.N(i, j);
      SSMADC.S(i, j) = ADC.S(i, j) * blank.S(i, j);
      SSMADC.P(i, j) = ADC.P(i, j) + ((1 - blank.E(i, j)) * ADC.E(i, j) +
                                      (1 - blank.W(i, j)) * ADC.W(i, j) +
                                      (1 - blank.N(i, j)) * ADC.N(i, j) +
                                      (1 - blank.S(i, j)) * ADC.S(i, j));
      SSMADC.Pinv(i, j) = 1 / SSMADC.P(i, j);

      // modify SSMADC for neumann  bc
      // ap=bp - sum((1-id_nb)*b_nb)
      // a_nb=id_nb*b_nb
      SSMDEq.E(i, j) = DEq.E(i, j) * blank.E(i, j);
      SSMDEq.W(i, j) = DEq.W(i, j) * blank.W(i, j);
      SSMDEq.N(i, j) = DEq.N(i, j) * blank.N(i, j);
      SSMDEq.S(i, j) = DEq.S(i, j) * blank.S(i, j);
      SSMDEq.P(i, j) = DEq.P(i, j) - ((1 - blank.E(i, j)) * DEq.E(i, j) +
                                      (1 - blank.W(i, j)) * DEq.W(i, j) +
                                      (1 - blank.N(i, j)) * DEq.N(i, j) +
                                      (1 - blank.S(i, j)) * DEq.S(i, j));
      SSMDEq.Pinv(i, j) = 1 / SSMDEq.P(i, j);
    }
  }

  // cout<<endl<<DEq.N<<endl;
  // cout<<endl<<SSMDEq.N<<endl;

  outputout << "\n"
            << "Number of boundary Nodes: " << boundarynumber << "\n";
  //   cout << "\n" << blank.ibnode << "\n";

  ssmlist.assign(boundarynumber);
  issm = 0;
  // Dont do this in parallel. Might screw up
  for (j = j_s; j <= j_e; j++) {
    for (i = i_s; i <= i_e; i++) {
      if (blank.ibnode(i, j)) {
        ssmlist.i(issm) = i;
        ssmlist.j(issm) = j;
        ssmlist.xc(issm) = xc(i);
        ssmlist.yc(issm) = yc(j);
        issm++;
      }
    }
  }
}
