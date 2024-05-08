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
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;
using namespace INPUTDATA;
using namespace CLASSES;
using namespace FIELDVARS;
using namespace COEFFVARS;
using namespace IBMVARS;
using namespace Eigen;

void ssm_bc_dir(Eigen::ArrayXXd &u,
                CLASSES::NeighbourField &U) // aply dirichlet bc
{
  U.E = u;
  U.N = u;
  U.S = u;
  U.W = u;
  U.P = u;
#pragma omp parallel for private(issm,i,j) shared(U, ssmlist)
  for (issm = 0; issm < ssmlist.Nb; issm++) {
    i = ssmlist.i(issm);
    j = ssmlist.j(issm);
    // dirichlet bc is u_g=-u_nb where g is ghost node value and nb is neighbour
    U.E(i, j) = -u(i - 1, j); // east; (i,j) is east of (i-1,j).
    U.W(i, j) = -u(i + 1, j); // west; (i,j) is west of (i+1,j).
    U.N(i, j) = -u(i, j - 1); // north; (i,j) is north of (i,j-1).
    U.S(i, j) = -u(i, j + 1); // south; (i,j) is east of (i,j+1).
    U.P(i, j) =
        0.25 * (U.E(i, j) + U.W(i, j) + U.N(i, j) +
                U.S(i, j)); // dont need this but will ensure that
                            // diffusion term is enforced zero on boundary cells
    u(i, j) =
        -((1. - blank.E(i, j)) * u(i + 1, j) + // ADHOC Measure
          (1. - blank.W(i, j)) * u(i - 1, j) +
          (1. - blank.N(i, j)) * u(i, j + 1) +
          (1. - blank.S(i, j)) * u(i, j - 1)) /
        (4. - (blank.E(i, j) + blank.W(i, j) + blank.N(i, j) + blank.S(i, j)));
  }
}

void ssm_bc_neum(Eigen::ArrayXXd &u,
                 CLASSES::NeighbourField &U) // aply dirichlet bc
{
  U.E = u;
  U.N = u;
  U.S = u;
  U.W = u;
  U.P = u;
#pragma omp parallel for private(issm,i,j) shared(U, ssmlist)
  for (issm = 0; issm < ssmlist.Nb; issm++) {
    i = ssmlist.i(issm);
    j = ssmlist.j(issm);
    // dirichlet bc is u_g=-u_nb where g is ghost node value and nb is neighbour
    U.E(i, j) = u(i - 1, j); // east; (i,j) is east of (i-1,j).
    U.W(i, j) = u(i + 1, j); // west; (i,j) is west of (i+1,j).
    U.N(i, j) = u(i, j - 1); // north; (i,j) is north of (i,j-1).
    U.S(i, j) = u(i, j + 1); // south; (i,j) is east of (i,j+1).
    U.P(i, j) =
        0.25 * (U.E(i, j) + U.W(i, j) + U.N(i, j) +
                U.S(i, j)); // dont need this but will ensure that
                            // diffusion term is enforced zero on boundary cells
    u(i, j) =
        ((1. - blank.E(i, j)) * u(i + 1, j) + // ADHOC Measure
         (1. - blank.W(i, j)) * u(i - 1, j) +
         (1. - blank.N(i, j)) * u(i, j + 1) +
         (1. - blank.S(i, j)) * u(i, j - 1)) /
        (4. - (blank.E(i, j) + blank.W(i, j) + blank.N(i, j) + blank.S(i, j)));
  }
}

void ssm_flux_bc(Eigen::ArrayXXd &ci,
                 Eigen::ArrayXXd &cj) { // apply zero flux on body
#pragma omp parallel for private(issm,i,j) shared(ci, cj, ssmlist)
  for (issm = 0; issm < ssmlist.Nb; issm++) {
    i = ssmlist.i(issm);
    j = ssmlist.j(issm);
    // dirichlet bc is u_g=-u_nb where g is ghost node value and nb is neighbour
    ci(i, j) = 0;
    ci(i - 1, j) = 0;
    cj(i, j) = 0;
    cj(i, j - 1) = 0;
  }
}

void ssm_adhoc_bc(Eigen::ArrayXXd &u,
                 Eigen::ArrayXXd &v,
                 Eigen::ArrayXXd &p) { // apply zero flux on body
#pragma omp parallel for private(j, i) collapse(2)
  for (j = j_s; j <= j_e; j++) {
    for (i = i_s; i <= i_e; i++) {
        if(((1 - blank.ibnode(i, j)) && (blank.blank(i, j)))) //statement means if inside ib body and not a boundary node
{
    u(i,j)=1;
    v(i,j)=1;
    p(i,j)=1;
}
    }
  }
}

void getforce(Eigen::ArrayXXd &p ) { // get force on body
#pragma omp parallel for private(issm,i,j) shared(p)
  for (issm = 0; issm < ssmlist.Nb; issm++) {
    i = ssmlist.i(issm);
    j = ssmlist.j(issm);
    // dirichlet bc is u_g=-u_nb where g is ghost node value and nb is neighbour
        #define pE p(i+1,j)
        #define pW p(i-1,j)
        #define pN p(i,j+1)
        #define pS p(i,j-1)
        #define idE blank.E(i, j)
        #define idW blank.W(i, j)
        #define idN blank.N(i, j)
        #define idS blank.S(i, j)
        ssmlist.fx(issm)=dely(j)*(-(1.0-idE)*pE+(1.0-idW)*pW);
        ssmlist.fy(issm)=delx(i)*(-(1.0-idN)*pN+(1.0-idS)*pS);
  }
  ssmlist.totalforce(ssm_xc,ssm_yc);
  forceout<<soltime<<"\t"<<ssmlist.forcex<<"\t"<<ssmlist.forcey<<"\t"<<ssmlist.moment<<"\n";
}