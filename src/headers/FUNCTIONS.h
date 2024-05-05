#ifndef FUNCTIONS_H // include guard
#define FUNCTIONS_H
#include <eigen3/Eigen/Dense>
#include "./CLASSES.h"

// read_input_data.cpp
void read_input_data();
// write_output_data.cpp
void write_output_data();
// setup_field_variables.cpp
void setup_field_variables();
void setup_coeff_field();
void setup_primitives();
// boundary_conditions.cpp
void VelBC(Eigen::ArrayXXd &u, Eigen::ArrayXXd &v);
void FaceFlux(Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj);
void PresBC(Eigen::ArrayXXd &p);
// calculators.cpp
void HalfCalc(Eigen::ArrayXXd &u, Eigen::ArrayXXd &u_i, Eigen::ArrayXXd &u_j);
void ConvCalc(Eigen::ArrayXXd &C, Eigen::ArrayXXd &u,
              Eigen::ArrayXXd &u_i, Eigen::ArrayXXd &u_j,
              Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj);
void DiffCalc(Eigen::ArrayXXd &D, Eigen::ArrayXXd &u, CLASSES::SolverCoeffs &A);
void MassCalc(Eigen::ArrayXXd &C, Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj);

void GradPFace(Eigen::ArrayXXd &p, Eigen::ArrayXXd &pgx_ci, Eigen::ArrayXXd &pgy_cj);
void GradPCenter(Eigen::ArrayXXd &p, Eigen::ArrayXXd &pgx, Eigen::ArrayXXd &pgy);

// veelocity_predictor.cpp
void VelocityPredictor(Eigen::ArrayXXd &us, Eigen::ArrayXXd &vs,
                       Eigen::ArrayXXd &RHSu, Eigen::ArrayXXd &RHSv,
                       Eigen::MatrixXd &Ru, Eigen::MatrixXd &Rv,
                       CLASSES::SolverCoeffs &B, int &kmom, double &errormom);
// solvers.cpp
void PointSOR(Eigen::ArrayXXd &p,
              Eigen::ArrayXXd &f,
              CLASSES::SolverCoeffs &A,
              Eigen::MatrixXd &Res,
              int &kSOR, double &error);
//fractionalstep.cpp
void fractional_step();
#endif
