  #include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>

// #include <string> //Wheres Wallace????make
// #include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/FIELDVARS.h"
// #include "./headers/COEFFVARS.h"
// #include "./headers/CLASSES.h"
#include "./headers/FLOWVARS.h"
// #include "headers/INPUT.h"

// #include <iostream>
// #include<iostream>
// using namespace INPUTDATA;
using namespace FIELDVARS;
using namespace Eigen;
// using namespace COEFFVARS;
// using namespace CLASSES;
using namespace std;
using namespace FLOWVARS;

void debug_export()
{
  ofstream outxAA("out/x.dat");
  outxAA  << x;

  ofstream outy("out/y.dat");
  outy <<  y;

  ofstream outxc("out/xc.dat");
  outxc <<  xc;

  ofstream outyc("out/yc.dat");
  outyc  << yc;

  ofstream outdx("out/dx.dat");
  outdx  << dx;

  ofstream outdy("out/dy.dat");
  outdy  << dy;

  ofstream outdelx("out/delx.dat");
  outdelx  << delx;
  ofstream outdely("out/dely.dat");
  outdely << dely;


  ofstream uout("out/u.dat");
  uout << New.u;
  ofstream usout("out/us.dat");
  usout << New.us;

  ofstream u0out("out/u0.dat");
  u0out << Old.u;

  ofstream vout("out/v.dat");
  vout << New.v;
  ofstream vsout("out/vs.dat");
  vsout << New.vs;

  ofstream v0out("out/v0.dat");
  v0out << Old.v;

  ofstream ciout("out/ci.dat");
  ciout << New.ci;
  ofstream cjout("out/cj.dat");
  cjout << New.cj;

  ofstream cisout("out/cis.dat");
  cisout << New.cis;
  ofstream cjsout("out/cjs.dat");
  cjsout << New.cjs;

  ofstream RHSuout("out/RHSu.dat");
  RHSuout << FSA.RHSu;
  ofstream RHSvout("out/RHSv.dat");
  RHSvout << FSA.RHSv;

  ofstream massdout("out/massd.dat");
  massdout << FSA.div_s;
  ofstream massdnewout("out/massdnew.dat");
  massdnewout << FSA.div_new;

  ofstream Duout("out/Du.dat");
  Duout << FSA.uD1;
  ofstream Dvout("out/Dv.dat");
  Dvout << FSA.vD1;

  ofstream xgpout("out/xgradp.dat");
  xgpout << FSA.gradpx;
  ofstream ygpout("out/ygradp.dat");
  ygpout << FSA.gradpy;

  ofstream fout("out/f.dat");
  fout << FSA.f;
}