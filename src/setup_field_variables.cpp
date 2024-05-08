// sets up field variables
#include <iostream>
#include <fstream>
#include <string> //Wheres Wallace????
#include <sys/types.h>
#include <sys/stat.h>
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/COEFFVARS.h"
#include <eigen3/Eigen/Dense>
#include "./headers/CLASSES.h"
#include "./headers/FLOWVARS.h"
#include "omp.h"

using namespace std;
using namespace INPUTDATA;
using namespace CLASSES;
using namespace FIELDVARS;
using namespace COEFFVARS;
using namespace Eigen;
using namespace FLOWVARS;

void setup_field_variables() //setup things like dx,dy,x,y,dx_inv... you get the idea
{
    // resize x,xc,y,ycv
    x.resize(Nx_f);
    xc.resize(Nx_c);
    y.resize(Ny_f);
    yc.resize(Ny_c);
    // resize dx, dy, delx dely
    dx.resize(Nx_f);
    dx_inv.resize(Nx_f);
    delx.resize(Nx_c);
    delx_inv.resize(Nx_c);

    dy.resize(Ny_f);
    dy_inv.resize(Ny_f);
    dely.resize(Ny_c);
    dely_inv.resize(Ny_c);

    // read x and y (face grid points)
    if (x_unif == 1) {
        x = ArrayXd::LinSpaced(Nx_f, x_start, x_end);
    } else {
        ifstream xgrid("xgrid.dat");
        for (i = if_s; i <= if_e; i++) {
            xgrid >> garbint >> x(i);
        }
        x_start = x.minCoeff(); // can use x(0). I just wanted to test minCoeff
        x_end = x.maxCoeff();
    }
    if (y_unif == 1) {
        y = ArrayXd::LinSpaced(Ny_f, y_start, y_end);
    } else {
        ifstream ygrid("ygrid.dat");
        for (j = jf_s; j <= jf_e; j++) {
            ygrid >> garbint >> y(j);
        }
        y_start = y.minCoeff();
        y_end = y.maxCoeff();
    }
    // create xc and yc (grid cell center points)
    for (i = i_s; i <= i_e; i++) {
        xc(i) = 0.5 * (x(i) + x(i - 1));
    }
    xc(0) = x(0) - (xc(1) - x(0));
    xc(Nx_c - 1) = x(last) + x(last) - xc(last - 1);

    // create xc and yc (grid cell center points)
    for (i = i_s; i <= i_e; i++) {
        xc(i) = 0.5 * (x(i) + x(i - 1));
    }
    xc(0) = x(0) - (xc(1) - x(0));
    xc.tail<1>()[0] = x(last) + x(last) - xc(last - 1);

    for (j = j_s; j <= j_e; j++) {
        yc(j) = 0.5 * (y(j) + y(j - 1));
    }
    yc(0) = y(0) - (yc(1) - y(0));
    yc.tail<1>()[0] = y(last) + y(last) - yc(last - 1);

    for (i = if_s; i <= if_e; i++) {
        dx(i) = xc(i + 1) - xc(i);
    }

    for (j = jf_s; j <= jf_e; j++) {
        dy(j) = yc(j + 1) - yc(j);
    }

    for (i = i_s; i <= i_e; i++)
        delx(i) = x(i) - x(i - 1);

    for (j = j_s; j <= j_e; j++)
        dely(j) = y(j) - y(j - 1);

    delx(0) = delx(1);
    delx.tail<1>()[0] = delx(last - 1);
    dely(0) = dely(1);
    dely.tail<1>()[0] = dely(last - 1);

    dx_inv = dx.cwiseInverse();
    dy_inv = dy.cwiseInverse();
    delx_inv = delx.cwiseInverse();
    dely_inv = dely.cwiseInverse();

    // cout<<"\n"<<"\n"<<"dxinv is:"<<"\n"<<"\n"<<delx_inv;
    // cout<<"\n"<<"\n"<<"dx is:"<<"\n"<<"\n"<<delx;

    // create dx and dy
    //  dx=0;
    // const int c_Nx_c=Nx_c;
    // const int c_Ny_c=Ny_c;

    // cout << "\n"<< dely << "\n";
    // blank.resize(Nx_c,Ny_c);
    // blank = Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(Nx_c, Ny_c, 0);
    // color_n = Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(Nx_c, Ny_c, 1);
    // color_s = Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(Nx_c, Ny_c, 1);
    // color_e = Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(Nx_c, Ny_c, 1);
    // color_w = Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(Nx_c, Ny_c, 1);
    // cout<<color_w<<"\n";
    // cout <<"\n"<< blank<<"\n";

    delx_min = delx.minCoeff();
    dely_min = dely.minCoeff();
    del_min = (delx_min < dely_min) ? delx_min : dely_min;
    CFL = dt / del_min;
    outputout << "Smallest grid size is:" << del_min << "\n";
    outputout << "CFL is:" << CFL << "\n";
}

void setup_coeff_field() {
    CoE.resize(Nx_f); // CDS half node interpolation coefficients//face values
    CoW.resize(Nx_f);
    CoN.resize(Ny_f);
    CoS.resize(Ny_f);

    for (i = if_s; i <= if_e; i++) {
        CoE(i) = delx(i) / (delx(i + 1) + delx(i));
        CoW(i) = delx(i + 1) / (delx(i + 1) + delx(i));
    }
    for (j = jf_s; j <= jf_e; j++) {
        CoN(j) = dely(j) / (dely(j + 1) + dely(j));
        CoS(j) = dely(j + 1) / (dely(j + 1) + dely(j));
    }


    GCOx = VectorXd::Zero(Nx_c); // Coefficients for gradient calculations
    GCOy = VectorXd::Zero(Ny_c); // Coefficients for gradient calculations

    for (i = i_s; i <= i_e; i++)
        GCOx(i) = 1 / (dx(i) + dx(i - 1));

    for (j = j_s; j <= j_e; j++)
        GCOy(j) = 1 / (dy(j) + dy(j - 1));

    // cout << GCOy;


    ADC.assign(Nx_c, Ny_c); // Advection Diffusion Coefficients
    DEq.assign(Nx_c, Ny_c); // Poisson Coefficients
    for (j = j_s; j <= j_e; j++) {
        for (i = i_s; i <= i_e; i++) {
            DEq.N(i, j) = 1 / (dely(j) * dy(j));
            DEq.S(i, j) = 1 / (dely(j) * dy(j - 1));
            DEq.E(i, j) = 1 / (delx(i) * dx(i));
            DEq.W(i, j) = 1 / (delx(i) * dx(i - 1));
            DEq.P(i, j) = DEq.N(i, j) + DEq.S(i, j) + DEq.E(i, j) + DEq.W(i, j);
            DEq.Pinv(i, j) = 1 / DEq.P(i, j);

            ADC.N(i, j) = -DEq.N(i, j) * 0.5 * Reinv;
            ADC.S(i, j) = -DEq.S(i, j) * 0.5 * Reinv;
            ADC.W(i, j) = -DEq.W(i, j) * 0.5 * Reinv;
            ADC.E(i, j) = -DEq.E(i, j) * 0.5 * Reinv;
            ADC.P(i, j) = -(dtinv + DEq.P(i, j) * 0.5 * Reinv);
            ADC.Pinv(i, j) = 1 / ADC.P(i, j);
        }
    }
    // cout<<delx;
    // cout<<ADC.W;
}

void setup_boundaryconditions() {

  outputout << "\n" << "\n" << " EAST FACE" << "\n";
  switch (bctype_e) { // East boundary conditions //same as west
  case 2:             // symmetric
    outputout << "\n" << "  Symmetric BC" << "\n";
    ubctype_e = 1; // Dir //u=0
    vbctype_e = 0; // Neumann //dv/dn=0
    break;
  case 3: // Inflow
    outputout << "\n" << "  Inflow BC" << "\n";
    ubctype_e = 1; // Dir //u=0
    vbctype_e = 1; // Dir
    break;
  case 4: // Outflow
    outputout << "\n" << "  Outflow BC" << "\n";
    ubctype_e = 0; // Neumann //du/dn=0
    vbctype_e = 0; // Neumann //dv/dn=0
    break;
  default: // wall
    outputout << "\n" << "  Wall BC" << "\n";
    ubctype_e = 1; // Dir
    vbctype_e = 1; // Dir
    break;
  }
  outputout << "\n" << " ubctype = " << ubctype_e << "\n";
  outputout << "\n" << " vbctype = " << vbctype_e << "\n";

  outputout << "\n" << "\n" << " WEST FACE" << "\n";
  switch (bctype_w) { // west boundary conditions
  case 2:             // symmetric
    outputout << "\n" << "  Symmetric BC" << "\n";
    ubctype_w = 1; // Dir //u=0
    vbctype_w = 0; // Neumann //dv/dn=0
    break;
  case 3: // Inflow
    outputout << "\n" << "  Inflow BC" << "\n";
    ubctype_w = 1; // Dir //u=0
    vbctype_w = 1; // Dir
    break;
  case 4: // Outflow
    outputout << "\n" << "  Outflow BC" << "\n";
    ubctype_w = 0; // Neumann //du/dn=0
    vbctype_w = 0; // Neumann //dv/dn=0
    break;
  default: // wall
    outputout << "\n" << "  Wall BC" << "\n";
    ubctype_w = 1; // Dir
    vbctype_w = 1; // Dir
    break;
  }
  outputout << "\n" << " ubctype = " << ubctype_w << "\n";
  outputout << "\n" << " vbctype = " << vbctype_w << "\n";

  outputout << "\n" << "\n" << " NORTH FACE" << "\n";
  switch (bctype_n) { // South boundary conditions
  case 2:             // symmetric //only symmetric gets changed
    outputout << "\n" << "  Symmetric BC" << "\n";
    ubctype_n = 0; // Neumann //du/dn=0
    vbctype_n = 1; // Dir //v=0
    break;
  case 3: // Inflow
    outputout << "\n" << "  Inflow BC" << "\n";
    ubctype_n = 1; // Dir //u=0
    vbctype_n = 1; // Dir
    break;
  case 4: // Outflow
    outputout << "\n" << "  Outflow BC" << "\n";
    ubctype_n = 0; // Neumann //du/dn=0
    vbctype_n = 0; // Neumann //dv/dn=0
    break;
  default: // wall
    outputout << "\n" << "  Wall BC" << "\n";
    ubctype_n = 1; // Dir
    vbctype_n = 1; // Dir
    break;
  }
  outputout << "\n" << " ubctype = " << ubctype_n << "\n";
  outputout << "\n" << " vbctype = " << vbctype_n << "\n";

  outputout << "\n" << "\n" << " SOUTH FACE" << "\n";
  switch (bctype_s) { // South boundary conditions
  case 2:             // symmetric //only symmetric gets changed
    outputout << "\n" << "  Symmetric BC" << "\n";
    ubctype_s = 0; // Neumann //du/dn=0
    vbctype_s = 1; // Dir //v=0
    break;
  case 3: // Inflow
    outputout << "\n" << "  Inflow BC" << "\n";
    ubctype_s = 1; // Dir //u=0
    vbctype_s = 1; // Dir
    break;
  case 4: // Outflow
    outputout << "\n" << "  Outflow BC" << "\n";
    ubctype_s = 0; // Neumann //du/dn=0
    vbctype_s = 0; // Neumann //dv/dn=0
    break;
  default: // wall
    outputout << "\n" << "  Wall BC" << "\n";
    ubctype_s = 1; // Dir
    vbctype_s = 1; // Dir
    break;
  }
  outputout << "\n" << " ubctype = " << ubctype_s << "\n";
  outputout << "\n" << " vbctype = " << vbctype_s << "\n";
}


void setup_primitives() {
    //Initialize fields associated with frac step
    FSA.InitializeCenters(Nx_c, Ny_c);
    FSA.InitializeFaces(Nx_f, Ny_f);

    //Initialize flow fields
    New.InitializeCenters(Nx_c, Ny_c, uinit, vinit, perturb);
    New.InitializeFaces(Nx_f, Ny_f, uinit, vinit);
    Old.InitializeCenters(Nx_c, Ny_c, uinit, vinit, perturb);
    Old.InitializeFaces(Nx_f, Ny_f, uinit, vinit);
    Older.InitializeCenters(Nx_c, Ny_c, uinit, vinit, perturb);
    Older.InitializeFaces(Nx_f, Ny_f, uinit, vinit);
    VelBC(New.u, New.v);
    FaceFluxBC(New.ci, New.cj);
    PresBC(New.p);
    VelBC(Old.u, Old.v);
    FaceFluxBC(Old.ci, Old.cj);
    PresBC(Old.p);
    VelBC(Older.u, Older.v);
    FaceFluxBC(Older.ci, Older.cj);
    PresBC(Older.p);
    // cout<<New.u_j;
    HalfCalc(New.u, New.u_i, New.u_j);
    // cout<<New.u_j;
}