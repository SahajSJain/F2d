#ifndef INPUT_H // include guard
#define INPUT_H

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
    inline std::string outputfilename = "output.dat";
    inline std::string outputfname = outputfilename;
    inline double dt, Re, Reinv, dtinv;
    inline int Tn, nprint, ndump, ndumpinit;
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
    inline double uinit, vinit,perturb, force_x, force_y;
    inline  std::ofstream outputout("output.dat"); // read from input.dat
    inline int i_s, i_e, j_s, j_e; //Starting and ending index for centers
    inline int if_s, if_e, jf_s, jf_e; //starting and ending index for faces
    inline int Nx_c, Ny_c; //No. of cell centers
    inline int Nx_f, Ny_f; //No. of cell faces
    inline int i{0}, j{0}, n{0};
    inline int garbint{0};
    
}
#endif
