#include <eigen3/Eigen/Dense>
// #include <string> //Wheres Wallace????make
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/COEFFVARS.h"
#include "./headers/CLASSES.h"

// #include <iostream>
// #include<iostream>
using namespace INPUTDATA;
using namespace FIELDVARS;
using namespace Eigen;
using namespace COEFFVARS;
using namespace CLASSES;
#include <iostream>

// Calculates things on faces
void HalfCalc(Eigen::ArrayXXd &u, Eigen::ArrayXXd &u_i, Eigen::ArrayXXd &u_j)
{

    u_i *= 0;
    u_j *= 0;
    // std::cout<<u_j;
    for (j = jf_s; j <= jf_e; j++)
    {
        for (i = if_s; i <= if_e; i++)
        {
            u_i(i, j) = u(i, j) * CoW(i) + u(i + 1, j) * CoE(i);
            u_j(i, j) = u(i, j) * CoS(j) + u(i, j + 1) * CoN(j);
        }
    }
}

// Calculates convection
void ConvCalc(Eigen::ArrayXXd &C, Eigen::ArrayXXd &u,
              Eigen::ArrayXXd &u_i, Eigen::ArrayXXd &u_j,
              Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj)
{
    HalfCalc(u, u_i, u_j);
    C *= 0;
    for (j = j_s; j <= j_e; j++)
    {
        for (i = i_s; i <= i_e; i++)
        {
            C(i, j) = (ci(i, j) * u_i(i, j) - ci(i - 1, j) * u_i(i - 1, j)) * delx_inv(i) +
                      (cj(i, j) * u_j(i, j) - cj(i, j - 1) * u_j(i, j - 1)) * dely_inv(j);
        }
    }
    // ci=ci;
    // cj=cj;
}

void DiffCalc(Eigen::ArrayXXd &D, Eigen::ArrayXXd &u, CLASSES::SolverCoeffs &A)
{
    D *= 0;
    for (j = j_s; j <= j_e; j++)
    {
        for (i = i_s; i <= i_e; i++)
        {
            D(i, j) = A.W(i, j) * u(i - 1, j) +
                      A.E(i, j) * u(i + 1, j) +
                      A.S(i, j) * u(i, j - 1) +
                      A.N(i, j) * u(i, j + 1) -
                      A.P(i, j) * u(i, j);
        }
    }
    // std::cout<<A.W;
}

void MassCalc(Eigen::ArrayXXd &C, Eigen::ArrayXXd &ci, Eigen::ArrayXXd &cj)
{
    C *= 0;
    for (j = j_s; j <= j_e; j++)
    {
        for (i = i_s; i <= i_e; i++)
        {
            C(i, j) = (ci(i, j) - ci(i - 1, j)) * delx_inv(i) +
                      (cj(i, j) - cj(i, j - 1)) * dely_inv(j);
        }
    }
}

void GradPFace(Eigen::ArrayXXd &p, Eigen::ArrayXXd &pgx_ci, Eigen::ArrayXXd &pgy_cj)
{
    pgx_ci *= 0;
    pgy_cj *= 0;
    for (j = jf_s; j <= jf_e; j++)
    {
        for (i = if_s; i <= if_e; i++)
        {
            pgx_ci(i, j) = (p(i + 1, j) - p(i, j)) * dx_inv(i);
            pgy_cj(i, j) = (p(i, j + 1) - p(i, j)) * dy_inv(j);
        }
    }
}

void GradPCenter(Eigen::ArrayXXd &p, Eigen::ArrayXXd &pgx, Eigen::ArrayXXd &pgy)
{
    pgx *= 0;
    pgy *= 0;

    for (j = j_s; j <= j_e; j++)
    {
        for (i = i_s; i <= i_e; i++)
        {
            pgx(i, j) = (p(i + 1, j) - p(i - 1, j)) * GCOx(i);
            pgy(i, j) = (p(i, j + 1) - p(i, j - 1)) * GCOy(i);
        }
    }
}


// void Equator(CLASSES::FractionalStepArrays &A,CLASSES::FractionalStepArrays &B)
// {//Set A=B
//     A.u=B.u;
//     A.v=B.v;
//     A.us=B.us;
//     A.vs=B.vs;

// }
