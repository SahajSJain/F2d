#pragma once
#include "./CLASSES.h"
#include <eigen3/Eigen/Dense>
namespace COEFFVARS{
    inline CLASSES::SolverCoeffs ADC; // Coeffs for Advections diffusion equation
    inline CLASSES::SolverCoeffs DEq; // Coeffs for
    inline Eigen::VectorXd GCOx(1), GCOy(1); //inverse of dx(i)+dx(i+1); used for gradient calculations
    inline Eigen::VectorXd CoE(1), CoW(1), CoN(1), CoS(1); //Coefficients for Half node calculations
}