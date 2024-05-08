#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
// #include <string> //Wheres Wallace????make
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/COEFFVARS.h"
#include "./headers/CLASSES.h"
#include "./headers/FLOWVARS.h"
#include "headers/INPUT.h"
// #include <iostream>
// #include<iostream>
using namespace INPUTDATA;
using namespace FIELDVARS;
using namespace Eigen;
using namespace COEFFVARS;
using namespace CLASSES;
using namespace std;
using namespace FLOWVARS;

void tecplot_printer() {
    //    outfile
    std::stringstream ss;
    ss << std::setw(6) << std::setfill('0') << dumpstep;
    std::string s = ss.str();
    std::ofstream tecout(outtecdir + "dump." + s + ".dat"); // read from input.dat
    tecout << "Title = \"F2d\" \n\n";
    tecout << "VARIABLES = \"x\" \"y\" \"z\" \"u\" \"v\" \"p\" \n\n";
    tecout << "ZONE I=" <<(Nx_c) << " J=" <<(Ny_c) <<" K=1, F=POINT \n\n";
    for (j = 0; j < Ny_c; j++) {
        for (i = 0; i < Nx_c; i++) {
            tecout << xc(i) << "\t" << yc(j) << "\t 0 \t" << New.u(i, j) << "\t" << New.v(i, j) << "\t" << New.p(i, j) << "\n\n";
        }
    }
}