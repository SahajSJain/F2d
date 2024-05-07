// my_program.cpp
// libraries to include
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>
#include <string> //Wheres Wallace????
#include <vector>

// #include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>

// #include<mdspan>
// header files

#include "./headers/CLASSES.h"
#include "./headers/COEFFVARS.h"
#include "./headers/FIELDVARS.h"
#include "./headers/FLOWVARS.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/INPUT.h"
#include "headers/INPUT.h"
#include <iomanip>
#include <sstream>

using namespace std;
using namespace INPUTDATA;
using namespace CLASSES;
using namespace FIELDVARS;
using namespace FLOWVARS;
using namespace COEFFVARS;

// using namespace Eigen;
void funky(CLASSES::SolverCoeffs &A) { cout << A.E; }

int main() {
  // ios_base::sync_with_stdio(false);
  // cin.tie(NULL);

  read_input_data();
  write_output_data();
  setup_field_variables();
  setup_coeff_field();
  setup_boundaryconditions();
  setup_primitives();
  fractional_step_looper();
  if (debugexportflag == 1)
    debug_export();

  return 0;
}
