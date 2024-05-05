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

void VelocityPredictor(Eigen::ArrayXXd &us, Eigen::ArrayXXd &vs,
                       Eigen::ArrayXXd &RHSu, Eigen::ArrayXXd &RHSv,
                       Eigen::MatrixXd &Ru, Eigen::MatrixXd &Rv,
                       CLASSES::SolverCoeffs &B, int &kmom, double &errormom)
{
    kmom = 1;
    errormom = 100000.0;
    // int kmin = 5;
    // int kmommax = adkmax;
    // double errormin = eps_ad;
    Ru *= 0;
    Rv *= 0;
    for (kmom = 1; kmom <= adkmax; kmom++)
    {
        kmom = kmom + 1;
        for (j = j_s; j <= j_e; j++)
        {
            for (i = i_s; i <= i_e; i++)
            {
                us(i, j) = -RHSu(i, j) +
                           B.W(i, j) * us(i - 1, j) + B.E(i, j) * us(i + 1, j) +
                           B.S(i, j) * us(i, j - 1) + B.N(i, j) * us(i, j + 1);
                us(i, j) = us(i, j) * B.Pinv(i, j);

                vs(i, j) = -RHSv(i, j) +
                           B.W(i, j) * vs(i - 1, j) + B.E(i, j) * vs(i + 1, j) +
                           B.S(i, j) * vs(i, j - 1) + B.N(i, j) * vs(i, j + 1);
                vs(i, j) = vs(i, j) * B.Pinv(i, j);
            }
        }
        VelBC(us, vs); // impose boundary conditions on us and vs

        Ru *= 0;
        Rv *= 0;

        for (j = j_s; j <= j_e; j++)
        {
            for (i = i_s; i <= i_e; i++)
            {
                Ru(i, j) = -RHSu(i, j) +
                           B.W(i, j) * us(i - 1, j) + B.E(i, j) * us(i + 1, j) +
                           B.S(i, j) * us(i, j - 1) + B.N(i, j) * us(i, j + 1) -
                           B.P(i, j) * us(i, j);
                Rv(i, j) = -RHSv(i, j) +
                           B.W(i, j) * vs(i - 1, j) + B.E(i, j) * vs(i + 1, j) +
                           B.S(i, j) * vs(i, j - 1) + B.N(i, j) * vs(i, j + 1) -
                           B.P(i, j) * vs(i, j);
            }
        }
        errormom = (Ru.norm() + Rv.norm()) / (Nx * Ny);
        if (errormom < eps_ad && kmom > 5)
            break;
    }
}