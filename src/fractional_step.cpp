#include <eigen3/Eigen/Dense>
#include <iostream>
// #include <string> //Wheres Wallace????make
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/COEFFVARS.h"
#include "./headers/CLASSES.h"
#include "./headers/FLOWVARS.h"

// #include <iostream>
// #include<iostream>
using namespace INPUTDATA;
using namespace FIELDVARS;
using namespace Eigen;
using namespace COEFFVARS;
using namespace CLASSES;
using namespace std;
using namespace FLOWVARS;

void fractional_step()
{
    // start iterations, non VK
    Older = Old;
    Old = New;
    New.us = New.u; /// initializa
    New.vs = New.v; /// initialize

    // Assemble RHS for us
    DiffCalc(FSA.uD1, Old.u, ADC); // Calculate D(u_n)
    ConvCalc(FSA.uC1, Old.u, Old.u_i,
             Old.u_j, Old.ci, Old.cj); // Calculates NL(u_n)
    ConvCalc(FSA.uC0, Older.u, Older.u_i,
             Older.u_j, Older.ci, Older.cj); // Calculates NL(u_n-1)

    FSA.RHSu = -1.5 * FSA.uC1 +
               0.5 * FSA.uC0 +
               0.5 * Reinv * FSA.uD1 +
               Old.u * dtinv;
    FSA.RHSu += force_x; // apply force in x firection
    // Assemble RHS for vs
    DiffCalc(FSA.vD1, Old.v, ADC); // Calculate D(u_n)
    ConvCalc(FSA.vC1, Old.v, Old.u_i,
             Old.u_j, Old.ci, Old.cj); // Calculates NL(u_n)
    ConvCalc(FSA.vC0, Older.v, Older.u_i,
             Older.u_j, Older.ci, Older.cj); // Calculates NL(u_n-1)

    FSA.RHSv = -1.5 * FSA.vC1 +
               0.5 * FSA.vC0 +
               0.5 * Reinv * FSA.vD1 +
               Old.v * dtinv;
    FSA.RHSu += force_y; // apply force in x firection

    // cout<<Old.v;

    // Solve for us and vs

    int kmom = 0;
    double errormom = 0;
    VelocityPredictor(New.us, New.vs,
                      FSA.RHSu, FSA.RHSv,
                      FSA.Res1, FSA.Res2,
                      ADC, kmom, errormom);
    // cout << endl<<New.ci<<endl<<endl<<endl;
    HalfCalc(New.us, New.cis, New.u_j); // Calculate cis
    HalfCalc(New.vs, New.u_i, New.cjs); // Calculate cjs
    // apply bc on intermediate face fluxes
    FaceFlux(New.cis, New.cjs);
    // calculate divergence of (cis,cjs) field
    MassCalc(FSA.div_s, New.cis, New.cjs); // div_s is the intermeediate divergence
    FSA.f = FSA.div_s * dtinv;
    // Solve for pressure
    int kPsolver = 0;
    double Perror = 0;
    PointSOR(New.p_pert, FSA.f, DEq,
             FSA.Res1, kPsolver, Perror);
    // Calculate pressure gradient
    GradPFace(New.p_pert, FSA.gradpx_f, FSA.gradpy_f);
    GradPCenter(New.p_pert, FSA.gradpx, FSA.gradpy);
    // Correct velocity and flux field
    New.u = New.us - dt * FSA.gradpx;
    New.v = New.vs - dt * FSA.gradpy;
    New.ci = New.cis - dt * FSA.gradpx_f;
    New.cj = New.cjs - dt * FSA.gradpy_f;
    // apply boundary conditions for corrected velocity field
    VelBC(New.u, New.v);
    // FaceFlux(New.ci, New.cj);
    MassCalc(FSA.div_new, New.ci, New.cj); // div_new is the corrected divergence
    // cout<<endl<<endl<<"Intermediate divergence field"<<endl<<endl<<FSA.div_s<<endl;
    // cout<<endl<<endl<<"Final divergence field"<<endl<<endl<<FSA.div_new<<endl;
    double massold = FSA.div_s.matrix().norm();
    double massnew = FSA.div_new.matrix().norm();
    cout << "Old Mass total: " << massold << endl;
    cout << "New Mass total: " << massnew << endl;
    cout << "Mom Iteration count:" << kmom << endl;
    cout << "Pressure Iteration count:" << kPsolver << endl;

    cout << "Mom Error:" << errormom << endl;
    cout << "Pressure Error:" << Perror << endl;
}