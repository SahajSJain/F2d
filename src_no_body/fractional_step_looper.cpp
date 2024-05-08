#include <eigen3/Eigen/Dense>
#include <iostream>
// #include <string> //Wheres Wallace????make
#include "./headers/FUNCTIONS.h"
#include "./headers/INPUT.h"

// #include "./headers/FIELDVARS.h"
#include "./headers/CLASSES.h"
#include "./headers/COEFFVARS.h"
#include "./headers/FLOWVARS.h"
#include "headers/INPUT.h"


// #include <iostream>
// #include<iostream>
using namespace INPUTDATA;
// using namespace FIELDVARS;
using namespace Eigen;
using namespace COEFFVARS;
using namespace CLASSES;
using namespace std;
using namespace FLOWVARS;

void fractional_step_looper() {
  dumpstep = 0;
  outputout << "\n"
            << "\n"
            << "Prining tecplot files at n = " << timestep << "\n";
  tecplot_printer();
  for (timestep = 0; timestep < Tn; timestep++) {
    soltime = soltime + dt;
    printflag = 0;
    if (timestep % nprint == 0)
      printflag = 1;
    //           if(timestep%nprint==0) printflag=1;
    //           if(printflag==0)
    if (printflag == 0)
    {
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
  }
}