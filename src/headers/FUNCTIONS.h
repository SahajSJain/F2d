#pragma once

//#include <eigen3/Eigen/Dense>
#include "./CLASSES.h"
#include <eigen3/Eigen/Dense>
// read_input_data.cpp
void read_input_data();
// write_output_data.cpp
void write_output_data();
// setup_field_variables.cpp
void setup_field_variables();
void setup_coeff_field();
void setup_primitives();
void setup_boundaryconditions();
// boundary_conditions.cpp
void VelBC(Eigen::ArrayXXd &u, Eigen::ArrayXXd &v);
void FaceFluxBC(Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj);
void PresBC(Eigen::ArrayXXd &p);
void FaceFluxBC_correction(Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj);
// calculators.cpp
void HalfCalc(Eigen::ArrayXXd &u, Eigen::ArrayXXd &u_i, Eigen::ArrayXXd &u_j);
void ConvCalc(Eigen::ArrayXXd &C, Eigen::ArrayXXd &u,
              Eigen::ArrayXXd &u_i, Eigen::ArrayXXd &u_j,
              Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj);
void DiffCalc(Eigen::ArrayXXd &D, CLASSES::NeighbourField &U, CLASSES::SolverCoeffs &A);
void MassCalc(Eigen::ArrayXXd &C, Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj);
void GradPFace(Eigen::ArrayXXd &p, Eigen::ArrayXXd &pgx_ci, Eigen::ArrayXXd &pgy_cj);
void GradPCenter(CLASSES::NeighbourField &P, Eigen::ArrayXXd &pgx, Eigen::ArrayXXd &pgy);

// veelocity_predictor.cpp
void VelocityPredictor(Eigen::ArrayXXd &us, Eigen::ArrayXXd &vs,
                       Eigen::ArrayXXd &RHSu, Eigen::ArrayXXd &RHSv,
                       CLASSES::SolverCoeffs &B, int &kmom, double &errormom);
// solvers.cpp
void PointSOR(Eigen::ArrayXXd &p,
              Eigen::ArrayXXd &f,
              CLASSES::SolverCoeffs &A,
              int &kSOR, double &error);
//fractionalstep.cpp
void fractional_step();

//tecplot_printer
void tecplot_printer();

//fractional_step_looper.cpp
void fractional_step_looper();

//debug_export.cp
void debug_export();

//ssm_Setup.cpp
void ssm_setup();

//ssm_bc.cpp
void ssm_bc_dir(Eigen::ArrayXXd &u,
                CLASSES::NeighbourField &U); // aply dirichlet bc

void ssm_bc_neum(Eigen::ArrayXXd &u,
                 CLASSES::NeighbourField &U); // aply dirichlet bc

void ssm_flux_bc(Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj); //flux bc

void ssm_adhoc_bc(Eigen::ArrayXXd &u,
                 Eigen::ArrayXXd &v,
                 Eigen::ArrayXXd &p);

void getforce(Eigen::ArrayXXd &p ); //force calculation
