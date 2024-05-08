#include <iostream>
#include <fstream>
#include <string> //Wheres Wallace????
#include <sys/types.h>
#include <sys/stat.h>
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"

using namespace std;
using namespace INPUTDATA;

void read_input_data()
{
    char comma;
    // First read in everything from input.dat
    ifstream inputin(inputfilename); // read from input.dat
    getline(inputin, garbage);       // ========= GRID DETAILS ==================================
    inputin >> Nx >> comma >> Ny;
    getline(inputin, garbage); // Nx       Ny
    inputin >> x_unif >> comma >> x_start >> comma >> x_end;
    getline(inputin, garbage); // x_unif  Nx_start Nx_end
    // cout << garbage << "\n";
    // cout << x_unif << " " << Nx_start << " " << Nx_end << "\n";
    inputin >> y_unif >> comma >> y_start >> comma >> y_end;
    getline(inputin, garbage); // y_unif  Ny_start Ny_end
    // cout << x_unif << " " << y_unif << "\n";
    // getline(inputin, garbage); // y_unif  Ny_start Ny_end
    getline(inputin, garbage); //======== SOLVER DETAILS =================================
    inputin >> dt >> comma >> Re;
    getline(inputin, garbage); // dt (dt=CFL*min(dx))     Re

    dtinv = 1.0 / dt;
    Reinv = 1.0 / Re;
    inputin >> Tn >> comma >> nprint >> comma >> ndump >> comma >> ndumpinit;
    getline(inputin, garbage); // Tn (No. of timesteps)   nprint  ndump   ndumpinit

    inputin >> vkflag;
    getline(inputin, garbage); // VonKan (vkflag)

    if (vkflag == 1)
        vkflagstring = "VanKan";
    else
    {
        vkflagstring = "NonVanKan";
        vkflag=0;
    }

    inputin >> debugflag>> comma >> debugexportflag;
    getline(inputin, garbage); // VonKan (vkflag)

    getline(inputin, garbage); // ~~~Advection Diffusion Solver~~~
    inputin >> adkmax >> comma >> eps_ad;
    getline(inputin, garbage); // adkmax  eps_ad

    getline(inputin, garbage); // ~~~Pressure Solver Solver~~~
    inputin >> PSolver >> comma >> kiner >> comma >> kmax >> comma >> eps_pres >> comma >> overrelax;
    getline(inputin, garbage); // PSolver kiner   kmax    eps_pres        overrelax

    getline(inputin, garbage); //~~~~ SSM STUFF ~~~~~~
    inputin>>ssmflag>>comma>>forceflag;
    getline(inputin, garbage);  //SSMFlag, forceflag
    inputin>>ssm_bodytype>>comma>>ssm_x1>>comma>>ssm_y1>>comma>>ssm_x2>>comma>>ssm_y2>>comma>>ssm_rad;
    getline(inputin, garbage);   //ssm_bodytype, ssm_x1, ssm_y1, ssm_x2, ssm_y2, ssm_rad

    getline(inputin, garbage); //== Boundary Conditions and Field == (only zero p grad supported for now. Solver written for robin bc)
    getline(inputin, garbage); //~~~ Left BC ~~~
    inputin >> bctype_w >> comma >> ubc_w >> comma >> vbc_w;
    getline(inputin, garbage); //...>
    // cout<<garbage;

    getline(inputin, garbage); //~~~ right BC ~~~
    // cout<<garbage;
    inputin >> bctype_e >> comma >> ubc_e >> comma >> vbc_e;
    getline(inputin, garbage); //...>
    // cout<<garbage;

    getline(inputin, garbage); //~~~ north BC ~~~
    inputin >> bctype_n >> comma >> ubc_n >> comma >> vbc_n;
    getline(inputin, garbage); //...>

    getline(inputin, garbage); //~~~ south BC ~~~
    inputin >> bctype_s >> comma >> ubc_s >> comma >> vbc_s;
    getline(inputin, garbage); //...>

    getline(inputin, garbage); //~~ Field initialization~~
    inputin >> uinit >> comma >> vinit >> comma >> perturb >> comma >> force_x >> comma >> force_y;
    getline(inputin, garbage); // uinit           vinit           force_x         force_y

    // cout << Nx << " " << Ny << " \n ";
        // Initialize number of cells and faces
    Nx_c = Nx + 2;
    Ny_c = Ny + 2;
    Nx_f = Nx + 1;
    Ny_f = Ny + 1;
    // determine start for for loops
    i_s = 1;
    j_s = 1;
    if_s = 0;
    jf_s = 0;
    // determine the end for for loops
    i_e = Nx;
    j_e = Ny;
    if_e = Nx;
    jf_e = Ny;
    if(ssmflag!=1) 
    {
        forceflag=0;
        ssmflag=0;
    }
}
