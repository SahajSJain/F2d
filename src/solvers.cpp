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

void PointSOR(Eigen::ArrayXXd &p,
              Eigen::ArrayXXd &f,
              CLASSES::SolverCoeffs &A,
              Eigen::MatrixXd &Res,
              int &kSOR, double &error)
{
    kSOR = 0;
    error = 1000000.0;
    Res *= 0;
    int kgarb = 0;
    double Ds = 0;
    double pt = 0;
    while (kSOR < kmax && error > eps_pres) // k is the number of total iterations
    {
        kSOR = kSOR + kiner;
        for (kgarb = 0; kgarb < kiner; kgarb++)
        {
            for (j = j_s; j <= j_e; j++)
            {
                for (i = i_s; i <= i_e; i++)
                {
                    Ds = A.W(i, j) * p(i - 1, j) + A.E(i, j) * p(i + 1, j) +
                         A.S(i, j) * p(i, j - 1) + A.N(i, j) * p(i, j + 1);
                    pt = (Ds - f(i, j)) * A.Pinv(i, j);
                    p(i, j) = (1 - overrelax) * p(i, j) + overrelax * pt;
                }
            }
        }
        Res *= 0;
        for (j = j_s; j <= j_e; j++)
        {
            for (i = i_s; i <= i_e; i++)
            {
                Ds = A.W(i, j) * p(i - 1, j) + A.E(i, j) * p(i + 1, j) +
                     A.S(i, j) * p(i, j - 1) + A.N(i, j) * p(i, j + 1);
                Res(i, j) = p(i, j) * A.P(i, j) - (Ds - f(i, j));
            }
        }
        error = Res.norm() / (Nx * Ny);
    }
}
