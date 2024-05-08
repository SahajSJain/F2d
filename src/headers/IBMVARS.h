#pragma once
#include "./CLASSES.h"
namespace IBMVARS{
    inline CLASSES::NeighbourField n_u, n_v, n_p;
    // inline CLASSES::NeighbourField n_u0, n_v0, n_p0;   //not needed
    inline CLASSES::SolverCoeffs SSMADC; // Coeffs for Advections diffusion equation
    inline CLASSES::SolverCoeffs SSMDEq; // Coeffs for
    inline CLASSES::BlankArray blank; 
    inline CLASSES::BoundaryArray ssmlist;
}
