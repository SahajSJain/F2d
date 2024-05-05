// my_program.cpp
// libraries to include
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <fstream>
#include <string> //Wheres Wallace????
#include <bits/stdc++.h>
#include <sys/types.h>
#include <sys/stat.h>
// header files
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/COEFFVARS.h"
#include "./headers/FLOWVARS.h"

#include "./headers/CLASSES.h"

#include "./TECIO/TECIO.h"
#include "./TECIO/MASTER.h" /* for defintion of NULL */
#include <sstream>
using namespace std;
using namespace INPUTDATA;
using namespace CLASSES;
using namespace FIELDVARS;
using namespace FLOWVARS;
using namespace COEFFVARS;

// using namespace Eigen;
void funky(CLASSES::SolverCoeffs &A)
{
    cout << A.E;
}

int main()
{
    read_input_data();
    write_output_data();
    setup_field_variables();
    setup_coeff_field();
    setup_primitives();
    // calculate half node fluxes
    fractional_step();
    // funky(ADC.W);
    return 0;
    // cout<<bool;
}
