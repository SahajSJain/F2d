#include "eigen3/Eigen/Dense"
#include <fstream>
#include <iostream>
// #include <string> //Wheres Wallace????make
#include "./headers/FUNCTIONS.h"
#include "./headers/INPUT.h"
#include <fstream>
// #include "./headers/FIELDVARS.h"
#include "./headers/CLASSES.h"
#include "./headers/COEFFVARS.h"
#include "./headers/FLOWVARS.h"
#include "./headers/FIELDVARS.h"

#include "headers/INPUT.h"
#include "./headers/IBMVARS.h"
// #include <iostream>
// #include<iostream>
using namespace INPUTDATA;
// using namespace FIELDVARS;
using namespace Eigen;
using namespace COEFFVARS;
using namespace CLASSES;
using namespace std;
using namespace FLOWVARS;
using namespace IBMVARS;
using namespace FIELDVARS;
void fractional_step_looper() {
  dumpstep = 0;
  outputout << "\n"
            << "\n"
            << "Prining tecplot files at n = " << timestep << "\n";
  tecplot_printer();
  Average.InitializeCenters(Nx_c, Ny_c, uinit, vinit, perturb);
  Average.u *= 0;
  Average.v *= 0;
  Average.p *= 0;
  Average.om *= 0;
  int numdump = 0;
  for (timestep = 0; timestep < Tn; timestep++) {
    soltime = soltime + dt;
    printflag = 0;
    if (timestep % nprint == 0)
      printflag = 1;
    //           if(timestep%nprint==0) printflag=1;
    //           if(printflag==0)
    if (printflag == 0) {
      outputout << "\n"
                << "\n"
                << "========================================" << "\n"
                << "TIMESTEP NO. = " << timestep << "\n"
                << "TIME = " << soltime << "\n"
                << "========================================" << "\n";
    }

    // Start simulatiom
    fractional_step(); // fractional step
    // Should we print tecplot files?
    if (timestep % ndump == 0) {
      dumpstep = dumpstep + 1;
      outputout << "\n"
                << "\n"
                << "Prining tecplot files at n = " << timestep << "\n";
      tecplot_printer();
    }
    if (printflag == 0)
      outputout << "========================================" << "\n";

// averaging
#define flagdump (timestep > ndumpinit)
    Average.u += New.u * flagdump;
    Average.v += New.v * flagdump;
    Average.p += New.p * flagdump;
    Average.om += New.u * flagdump;
    numdump += flagdump;
  }
  // print average files
  std::ofstream tecoutave("./out/dump.average.dat"); // read from input.dat
  tecoutave << "Title = \"F2d\" \n\n";
  tecoutave << "VARIABLES = \"x\" \"y\" \"z\" \"u\" \"v\" \"p\" \"blank\" "
            "\"ibnode\" \"Vorticity\"  \n\n ";
  tecoutave << "ZONE I=" << (Nx_c) << " J=" << (Ny_c) << " K=1, F=POINT \n\n";
  Average.u/=numdump;
  Average.v/=numdump;
  Average.p/=numdump;
  Average.om/=numdump;

  for (j = 0; j < Ny_c; j++) {
    for (i = 0; i < Nx_c; i++) {
      tecoutave << xc(i) << "\t" << yc(j) << "\t 0 \t" << Average.u(i, j) << "\t"
             << Average.v(i, j) << "\t" << Average.p(i, j) << "\t" << blank.blank(i, j)
             << "\t" << blank.ibnode(i, j) << "\t" << Average.om(i, j) << "\n\n";
    }
  }
}