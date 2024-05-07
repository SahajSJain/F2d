#pragma once


#include <string> //Wheres Wallace????
#include <fstream> //Wheres Wallace????

// using namespace std;
namespace INPUTDATA
{
    inline std::string inputfilename = "input.dat";
    inline int Nx{1}, Ny{1};
    inline bool x_unif{1}, y_unif{1};
    inline double x_start{0.}, y_start{0.};
    inline double x_end{1.}, y_end{1.};
    inline double CFL;
    inline std::string garbage = " ";
    inline std::string outputdir = "./out/";
//    inline std::string outputfilename = "./out/output.dat";
    inline std::string outtecdir="./out/tecplot/";
    inline std::string outputfname = "./out/output.dat";;
    inline double dt, Re, Reinv, dtinv;
    inline int Tn, nprint, ndump, ndumpinit,printstep{0},dumpstep{0};
    inline int velkiner{5};
    inline bool printflag=0;
    inline int timestep{0};
    inline double soltime{0.};
    inline bool vkflag{0};
    inline std::string vkflagstring{"vk"};
    inline int adkmax;
    inline double eps_ad;
    inline int PSolver, kiner, kmax;
    inline double eps_pres, overrelax;
    inline int bctype_e, bctype_w, bctype_n, bctype_s;
    inline double ubc_w, vbc_w;
    inline double ubc_n, vbc_n;
    inline double ubc_e, vbc_e;
    inline double ubc_s, vbc_s;
    inline bool ubctype_w, vbctype_w;
    inline bool ubctype_e, vbctype_e;
    inline bool ubctype_n, vbctype_n;
    inline bool ubctype_s, vbctype_s;
    
    inline double uinit, vinit,perturb, force_x, force_y;
    inline  std::ofstream outputout(outputfname); // read from input.dat
    inline int i_s, i_e, j_s, j_e; //Starting and ending index for centers
    inline int if_s, if_e, jf_s, jf_e; //starting and ending index for faces
    inline int Nx_c, Ny_c; //No. of cell centers
    inline int Nx_f, Ny_f; //No. of cell faces
    inline int i{0}, j{0}, n{0};
    inline int garbint{0};
    inline bool debugflag=0, debugexportflag=0;
    inline double fluxintegral=0, flowintegralpost=0, outflowarea=0;
    inline double inflowintegral=0, outflowintegral=0;
    inline double inflowintegralpost=0, outflowintegralpost=0;

}

