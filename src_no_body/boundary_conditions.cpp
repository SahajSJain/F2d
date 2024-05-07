#include <eigen3/Eigen/Dense>
// #include <string> //Wheres Wallace????make
#include "./headers/FUNCTIONS.h"
#include "./headers/INPUT.h"
#include "omp.h"

#include "./headers/FIELDVARS.h"
#include <iostream>
using namespace INPUTDATA;
using namespace FIELDVARS;
using namespace Eigen;
using namespace std;

void VelBC(Eigen::ArrayXXd &u,
           Eigen::ArrayXXd &v) { // applies boundary conditions on velocities
                                 // #pragma omp parallel
  // cout << ubc_n;

  {

#pragma omp parallel for private(j)
    for (j = 0; j < Ny_c; j++) {
      u(0, j) = (2.0 * ubc_w - u(1, j)) * ubctype_w +
                (1.0 - ubctype_w) * u(1, j); // left u
      v(0, j) = (2.0 * vbc_w - v(1, j)) * vbctype_w +
                (1.0 - vbctype_w) * v(1, j); // left v
    }
  }

  {

#pragma omp parallel for private(j)
    for (j = 0; j < Ny_c; j++) {
      u(Nx_c - 1, j) = (2.0 * ubc_e - u(Nx_c - 2, j)) * ubctype_e +
                       (1.0 - ubctype_e) * u(Nx_c - 2, j); // right u
      v(Nx_c - 1, j) = (2.0 * vbc_e - v(Nx_c - 2, j)) * vbctype_e +
                       (1.0 - vbctype_e) * v(Nx_c - 2, j); // right v
    }
  }

  {
#pragma omp parallel for private(i)
    for (i = 0; i < Nx_c; i++) {
      u(i, 0) = (2.0 * ubc_s - u(i, 1)) * ubctype_s +
                (1.0 - ubctype_s) * u(i, 1); // bottom row
      v(i, 0) = (2.0 * vbc_s - v(i, 1)) * vbctype_s +
                (1.0 - vbctype_s) * v(i, 1); // bottom row
    }
  }

  {
#pragma omp parallel for private(i)
    for (i = 0; i < Nx_c; i++) {
      u(i, Ny_c - 1) = (2.0 * ubc_n - u(i, Ny_c - 2)) * ubctype_n +
                       (1.0 - ubctype_n) * u(i, Ny_c - 2); // top row
      v(i, Ny_c - 1) = (2.0 * vbc_n - v(i, Ny_c - 2)) * vbctype_n +
                       (1.0 - vbctype_n) * v(i, Ny_c - 2); // top row
    }
  }
  // cout<<"\n"<<"ubctypen"<<ubctype_n;
  // cout<<"\n"<<(2.0 * ubc_n - 0) * ubctype_n +
  //  (1.0 - ubctype_n) * 0;
}

/// APPLY BOUNDARY CONDITIONS ON FLUXES
void FaceFluxBC(Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj) {
  // West boundary conditions //same as west
  // impose boundary condition
  switch (bctype_w) {
  case 3: // Inflow
#pragma omp parallel for private(j)
    for (j = 1; j < Ny_f; j++) {
      ci(0, j) = ubc_w;
    }
    break;
  case 4: // Outflow
    break;
  default: // wall or symmetric
#pragma omp parallel for private(j)
    for (j = 0; j < Ny_f; j++)
      ci(0, j) = 0; // left u
    break;
  }

  // EAST boundary conditions
  switch (bctype_e) {
  case 3: // Inflow
#pragma omp parallel for private(j)
    for (j = 1; j < Ny_f; j++) {
      ci(Nx_f - 1, j) = ubc_e;
    }
    break;
  case 4: // Outflow
    break;
  default: // wall or symmetric
#pragma omp parallel for private(j)
    for (j = 0; j < Ny_f; j++)
      ci(Nx_f - 1, j) = 0; // left u
    break;
  }

  // South boundary conditions
  switch (bctype_s) {
  case 3: // Inflow
#pragma omp parallel for private(i)
    for (i = 1; i < Nx_f; i++) {
      cj(i, 0) = ubc_s;
    }
    break;
  case 4: // Outflow //do nothing
    break;
  default: // wall or symmetric
#pragma omp parallel for private(i)
    for (i = 0; i < Nx_f; i++)
      cj(i, 0) = 0; // bottom u
    break;
  }

  // North boundary conditions
  switch (bctype_n) {
  case 3: // Inflow
#pragma omp parallel for private(i)
    for (i = 1; i < Nx_f; i++) {
      cj(i, Ny_f - 1) = ubc_n;
    }
    break;
  case 4: // Outflow //do nothing
    break;
  default: // wall or symmetric
#pragma omp parallel for private(i)
    for (i = 0; i < Nx_f; i++)
      cj(i, Ny_f - 1) = 0; // bottom u
    break;
  }
  FaceFluxBC_correction(ci,cj) ;
}
void PresBC(Eigen::ArrayXXd &p) {
  // #pragma omp parallel
  {
#pragma omp parallel for private(j)
    for (j = 0; j < Ny_c; j++) {
      p(0, j) = p(1, j);               // left u
      p(Nx_c - 1, j) = p(Nx_c - 2, j); // right u
    }
#pragma omp parallel for private(i)
    for (i = 0; i < Nx_c; i++) {
      p(i, 0) = p(i, 1);               // bottom row
      p(i, Ny_c - 1) = p(i, Ny_c - 2); // top row
    }
  }
}

void FaceFluxBC_correction(Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj) {
  fluxintegral = 0; // flow outgoing is positive
  outflowarea = 0;
  inflowintegral = 0;
  outflowintegral = 0;
  flowintegralpost = 0;
  inflowintegralpost = 0;
  outflowintegralpost = 0;

  // cout<<"\n"<<ubc_w<<"\n";

  // collect
  switch (bctype_w) {
  case 3: // Inflow
#pragma omp parallel for private(j) shared(ci)                                 \
    reduction(+ : fluxintegral, inflowintegral)
    for (j = 1; j < Ny_f; j++) {
      fluxintegral = fluxintegral - ci(0, j) * dely(j);
      inflowintegral = inflowintegral - ci(0, j) * dely(j);
    }
    break;
  case 4: // Outflow
#pragma omp parallel for private(j) shared(ci)                                 \
    reduction(+ : fluxintegral, outflowintegral, outflowarea)
    for (j = 1; j < Ny_f; j++) {
      fluxintegral = fluxintegral - ci(0, j) * dely(j);
      outflowintegral = outflowintegral - ci(0, j) * dely(j);
      outflowarea = outflowarea + dely(j);
    }
    break;
  default: // wall or symmetric
    break;
  }

  // EAST boundary conditions
  switch (bctype_e) {
  case 3: // Inflow
#pragma omp parallel for private(j) shared(ci)                                 \
    reduction(+ : fluxintegral, inflowintegral)
    for (j = 1; j < Ny_f; j++) {
      fluxintegral = fluxintegral + ci(Nx_f - 1, j) * dely(j);
      inflowintegral = inflowintegral + ci(Nx_f - 1, j) * dely(j);
    }
    break;
  case 4: // Outflow
#pragma omp parallel for private(j) shared(ci)                                 \
    reduction(+ : fluxintegral, outflowintegral, outflowarea)
    for (j = 1; j < Ny_f; j++) {
      fluxintegral = fluxintegral + ci(Nx_f - 1, j) * dely(j);
      outflowarea = outflowarea + dely(j);
      outflowintegral = outflowintegral + ci(Nx_f - 1, j) * dely(j);
    }
    break;
  default: // wall or symmetric
    break;
  }

  // South boundary conditions
  switch (bctype_s) {
  case 3: // Inflow
#pragma omp parallel for private(i) shared(cj)                                 \
    reduction(+ : fluxintegral, inflowintegral)
    for (i = 1; i < Nx_f; i++) {
      fluxintegral = fluxintegral - cj(i, 0) * delx(i);
      inflowintegral = inflowintegral - cj(i, 0) * delx(i);
    }
    break;
  case 4: // Outflow //do nothing
#pragma omp parallel for private(i) shared(cj)                                 \
    reduction(+ : fluxintegral, outflowintegral, outflowarea)
    for (i = 1; i < Nx_f; i++) {
      fluxintegral = fluxintegral - cj(i, 0) * delx(i);
      outflowarea = outflowarea + delx(i);
      outflowintegral = outflowintegral - cj(i, 0) * delx(i);
    }
    break;
  default: // wall or symmetric
    break;
  }

  // North boundary conditions
  switch (bctype_n) {
  case 3: // Inflow
#pragma omp parallel for private(i) shared(cj)                                 \
    reduction(+ : fluxintegral, inflowintegral)
    for (i = 1; i < Nx_f; i++) {
      fluxintegral = fluxintegral + cj(i, Ny_f - 1) * delx(i);
      inflowintegral = inflowintegral + cj(i, Ny_f - 1) * delx(i);
    }
    break;
  case 4: // Outflow //do nothing
#pragma omp parallel for private(i) shared(cj)                                 \
    reduction(+ : fluxintegral, outflowintegral, outflowarea)
    for (i = 1; i < Nx_f; i++) {
      fluxintegral = fluxintegral + cj(i, Ny_f - 1) * delx(i);
      outflowarea = outflowarea + delx(i);
      outflowintegral = outflowintegral + cj(i, Ny_f - 1) * delx(i);
    }
    break;
  default: // wall or symmetric
    break;
  }

  // Redistribute mass over outflow area
  double ofa_inv = 1 / outflowarea;
  // cout<<"\n"<<"OFA INV"<<ofa_inv<<"\n";
  // cout<<"\n"<<"Redistributions"<<fluxintegral * ofa_inv<<"\n";

  // West
  if (bctype_w == 4) {
#pragma omp parallel for private(j) shared(ci, fluxintegral, ofa_inv)
    for (j = 1; j < Ny_f; j++)
      ci(0, j) = ci(0, j) + fluxintegral * ofa_inv ;
  }
  // east
  if (bctype_e == 4) {
#pragma omp parallel for private(j) shared(ci, fluxintegral, ofa_inv)
    for (j = 1; j < Ny_f; j++)
      ci(Nx_f - 1, j) = ci(Nx_f - 1, j) - fluxintegral * ofa_inv ;
  }
  // South
  if (bctype_s == 4) {
#pragma omp parallel for private(i) shared(cj, fluxintegral, ofa_inv)
    for (i = 1; i < Nx_f; i++) 
      cj(i, 0) = cj(i, 0) + fluxintegral * ofa_inv ;
  }
  // North
  if (bctype_n == 4) {
#pragma omp parallel for private(i) shared(cj, fluxintegral, ofa_inv)
    for (i = 1; i < Nx_f; i++)
      cj(i, Ny_f - 1) = cj(i, Ny_f - 1) - fluxintegral * ofa_inv ;
    }

    // Recalculate to confirm
    switch (bctype_w) {
    case 3: // Inflow
#pragma omp parallel for private(j) shared(ci)                                 \
    reduction(+ : flowintegralpost, inflowintegralpost)
      for (j = 1; j < Ny_f; j++) {
        flowintegralpost = flowintegralpost - ci(0, j) * dely(j);
        inflowintegralpost = inflowintegralpost - ci(0, j) * dely(j);
      }
      break;
    case 4: // Outflow
#pragma omp parallel for private(j) shared(ci)                                 \
    reduction(+ : flowintegralpost, outflowintegralpost)
      for (j = 1; j < Ny_f; j++) {
        flowintegralpost = flowintegralpost - ci(0, j) * dely(j);
        outflowintegralpost = outflowintegralpost - ci(0, j) * dely(j);
      }
      break;
    default: // wall or symmetric
      break;
    }

    // EAST boundary conditions
    switch (bctype_e) {
    case 3: // Inflow
#pragma omp parallel for private(j) shared(ci)                                 \
    reduction(+ : flowintegralpost, inflowintegralpost)
      for (j = 1; j < Ny_f; j++) {
        flowintegralpost = flowintegralpost + ci(Nx_f - 1, j) * dely(j);
        inflowintegralpost = inflowintegralpost + ci(Nx_f - 1, j) * dely(j);
      }
      break;
    case 4: // Outflow
#pragma omp parallel for private(j) shared(ci)                                 \
    reduction(+ : flowintegralpost, outflowintegralpost)
      for (j = 1; j < Ny_f; j++) {
        flowintegralpost = flowintegralpost + ci(Nx_f - 1, j) * dely(j);
        outflowintegralpost = outflowintegralpost + ci(Nx_f - 1, j) * dely(j);
      }
      break;
    default: // wall or symmetric
      break;
    }

    // South boundary conditions
    switch (bctype_s) {
    case 3: // Inflow
#pragma omp parallel for private(i) shared(cj)                                 \
    reduction(+ : flowintegralpost, inflowintegralpost)
      for (i = 1; i < Nx_f; i++) {
        flowintegralpost = flowintegralpost - cj(i, 0) * delx(i);
        inflowintegralpost = inflowintegralpost - cj(i, 0) * delx(i);
      }
      break;
    case 4: // Outflow //do nothing
#pragma omp parallel for private(i) shared(cj)                                 \
    reduction(+ : flowintegralpost, outflowintegralpost)
      for (i = 1; i < Nx_f; i++) {
        flowintegralpost = flowintegralpost - cj(i, 0) * delx(i);
        outflowintegralpost = outflowintegralpost - cj(i, 0) * delx(i);
      }
      break;
    default: // wall or symmetric
      break;
    }

    // North boundary conditions
    switch (bctype_n) {
    case 3: // Inflow
#pragma omp parallel for private(i) shared(cj)                                 \
    reduction(+ : flowintegralpost, inflowintegralpost)
      for (i = 1; i < Nx_f; i++) {
        fluxintegral = fluxintegral + cj(i, Ny_f - 1) * delx(i);
        inflowintegralpost = inflowintegralpost + cj(i, Ny_f - 1) * delx(i);
      }
      break;
    case 4: // Outflow //do nothing
#pragma omp parallel for private(i) shared(cj)                                 \
    reduction(+ : flowintegralpost, outflowintegralpost)
      for (i = 1; i < Nx_f; i++) {
        fluxintegral = fluxintegral + cj(i, Ny_f - 1) * delx(i);
        outflowintegralpost = outflowintegralpost + cj(i, Ny_f - 1) * delx(i);
      }
      break;
    default: // wall or symmetric
      break;
    }
  }