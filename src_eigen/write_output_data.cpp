#include <iostream>
#include "eigen3/Eigen/Dense"
#include <vector>
#include <fstream>
#include <string> //Wheres Wallace????
#include <sys/types.h>
#include <sys/stat.h>
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"

using namespace std;
using namespace INPUTDATA;

void write_output_data() {
    //     outputdir = "./out.Nx" + to_string(Nx) + ".Ny" + to_string(Ny) + "." + to_string(x_unif) + "." + to_string(y_unif) + "." + vkflagstring + "/";
    // // Creating a directory
    // string dirName = "./out/";
    if (mkdir(outputdir.c_str(), 0777) == 0)
        cout << " Directory created:" << outputdir << "\n";
    else
        cout << " Directory already present:" << outputdir << "\n";
    
    if (mkdir(outtecdir.c_str(), 0777) == 0)
        cout << " Directory created:" << outputdir << "\n";
    else
        cout << " Directory already present:" << outputdir << "\n";

//    if (mkdir(tecoutdir.c_str(), 0777) == 0)
//        cout << " Directory created:" << outputdir << "\n";
//    else
//        cout << " Directory already present:" << outputdir << "\n";


//    outputfname = outputdir + outputfilename;
    outputout << "Readin in input.dat \n";
    outputout << "Nx = " << Nx << " , Ny = " << Ny << "\n";
    outputout << " Uniform in x direction:" << x_unif << "\n";
    outputout << " Uniform in y direction:" << y_unif << "\n";
    outputout << "(Only aplicable for Unif) Bounds for x direction: (" << x_start << "," << x_end << ")" << "\n";
    outputout << "(Only aplicable for Unif) Bounds for y direction: (" << y_start << "," << y_end << ")" << "\n";

    outputout << "\n\n\n Field Details" << "\n";
    outputout << "dt=" << dt << "\n";
    outputout << "Re=" << Re << "\n";
    outputout << " Simulation will be run for " << Tn << " timesteps" << "\n";
    outputout << " Print data at every  " << nprint << " timesteps" << "\n";
    outputout << " Dump tecplot files every " << ndump << " timesteps starting from n=" << ndumpinit << "\n";
    outputout << " Fractional Step Mode :" << vkflagstring << "\n";
    outputout << "Ad Diff Solver: kmax = " << adkmax << " till eps_ad=" << eps_ad << "\n";
    outputout << " Pressure Solver Mode :" << PSolver << "\n";
    outputout << "Pressure Solver: kmax = " << kmax << " till eps_pres=" << eps_pres << "\n";
    outputout << "Residual reported every = " << kiner << " iterations and overrelaxation factor is = " << overrelax << "\n";

    outputout << "Left Face" << "\n";
    outputout << "Type? " << bctype_w << "\n";
    outputout << " u vel=" << ubc_w << " v vel=" << vbc_w << "\n"
              << "\n";

    outputout << "Right Face" << "\n";
    outputout << "Type? " << bctype_e << "\n";
    outputout << " u vel=" << ubc_e << " v vel=" << vbc_e << "\n"
              << "\n";

    outputout << "North Face" << "\n";
    outputout << "Type? " << bctype_n << "\n";
    outputout << " u vel=" << ubc_n << " v vel=" << vbc_n << "\n"
              << "\n";

    outputout << "South Face" << "\n";
    outputout << "Type? " << bctype_s << "\n";
    outputout << " u vel=" << ubc_s << " v vel=" << vbc_s << "\n"
              << "\n";

    outputout << "Field Initialization: u=" << uinit << " v=" << vinit << "\n";
    outputout << "Perturbation=" << perturb << "\n"
              << "\n";
    outputout << " Body forces: X dir: " << force_x << " Y dir: " << force_y << "\n";

    outputout << "\n\t\t~~~~~ STARTING SIMULATION~~~~\n\n";
}
