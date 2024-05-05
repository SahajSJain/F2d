#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <string> //Wheres Wallace????
#include <bits/stdc++.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"

using namespace std;
using namespace INPUTDATA;

void write_output_data()
{
    // outputdir = "./out.Nx" + to_string(Nx) + ".Ny" + to_string(Ny) + "." + to_string(x_unif) + "." + to_string(y_unif) + "." + vkflagstring + "/";
    // // Creating a directory
    // string dirName = "./out/";
    // if (mkdir(outputdir.c_str(), 0777) == 0)
    //     cout << " Directory created:" << outputdir << "\n";
    // else
    //     cout << " Directory already present:" << outputdir << "\n";

    outputfname = outputdir + outputfilename;
    outputout << "Readin in input.dat \n";
    outputout << "Nx = " << Nx << " , Ny = " << Ny << endl;
    outputout << " Uniform in x direction:" << x_unif << endl;
    outputout << " Uniform in y direction:" << y_unif << endl;
    outputout << "(Only aplicable for Unif) Bounds for x direction: (" << x_start << "," << x_end << ")" << endl;
    outputout << "(Only aplicable for Unif) Bounds for y direction: (" << y_start << "," << y_end << ")" << endl;

    outputout << "\n\n\n Field Details" << endl;
    outputout << "dt=" << dt << endl;
    outputout << "Re=" << Re << endl;
    outputout << " Simulation will be run for " << Tn << " timesteps" << endl;
    outputout << " Print data at every  " << nprint << " timesteps" << endl;
    outputout << " Dump tecplot files every " << ndump << " timesteps starting from n=" << ndumpinit << endl;
    outputout << " Fractional Step Mode :" << vkflagstring << endl;
    outputout << "Ad Diff Solver: kmax = " << adkmax << " till eps_ad=" << eps_ad << endl;
    outputout << " Pressure Solver Mode :" << PSolver << endl;
    outputout << "Pressure Solver: kmax = " << kmax << " till eps_pres=" << eps_pres << endl;
    outputout << "Residual reported every = " << kiner << " iterations and overrelaxation factor is = " << overrelax << endl;

    outputout << "Left Face" << endl;
    outputout << "Type? " << bctype_w << endl;
    outputout << " u vel=" << ubc_w << " v vel=" << vbc_w << endl
              << endl;

    outputout << "Right Face" << endl;
    outputout << "Type? " << bctype_e << endl;
    outputout << " u vel=" << ubc_e << " v vel=" << vbc_e << endl
              << endl;

    outputout << "North Face" << endl;
    outputout << "Type? " << bctype_n << endl;
    outputout << " u vel=" << ubc_n << " v vel=" << vbc_n << endl
              << endl;

    outputout << "South Face" << endl;
    outputout << "Type? " << bctype_s << endl;
    outputout << " u vel=" << ubc_s << " v vel=" << vbc_s << endl
              << endl;

    outputout << "Field Initialization: u=" << uinit << " v=" << vinit << endl;
    outputout << "Perturbation="<<perturb<<endl<<endl;
    outputout << " Body forces: X dir: " << force_x << " Y dir: " << force_y << endl;

    outputout << "\n\t\t~~~~~ STARTING SIMULATION~~~~\n\n";
}
