#include <eigen3/Eigen/Dense>
// #include <string> //Wheres Wallace????make
#include "./headers/INPUT.h"
#include "./headers/FUNCTIONS.h"
#include "./headers/FIELDVARS.h"
// #include<iostream>
using namespace INPUTDATA;
using namespace FIELDVARS;
using namespace Eigen;
// using namespace std;
void VelBC(Eigen::ArrayXXd& u,Eigen::ArrayXXd& v)
{
    u.row(0)=-u.row(1);//bottom row(i)
    v.row(0)=-v.row(1);//bottom row(i)
    u.row(u.rows()-1)=2-u.row(u.rows()-2);//top row(i)
    v.row(v.rows()-1)=-v.row(v.rows()-2);//top row(i)
    
    u.col(0)=-u.col(1);//left row(i)
    v.col(0)=-v.col(1);//left row(i)
    u.col(u.cols()-1)=-u.col(u.cols()-2);//right row(i)
    v.col(v.cols()-1)=-v.col(v.cols()-2);//right row(i)
    // v.row(last)=-v.row(last-1);//bottom row(i)
    // cout<<endl<<endl<<u(u.rows()-1,0)<<endl<<endl;
}

void FaceFlux(Eigen::ArrayXXd& ci,Eigen::ArrayXXd& cj)
{
    ci.row(0)=0;//bottom row(i)
    ci.row(cj.rows()-1)=0;//top row(i)
    
    cj.col(0)=0;//left row(i)
    cj.col(cj.cols()-1)=0;//right row(i)
    // v.row(last)=-v.row(last-1);//bottom row(i)
    // cout<<endl<<endl<<u(u.rows()-1,0)<<endl<<endl;
}


void PresBC(Eigen::ArrayXXd& p)
{
    p.row(0)=-p.row(1);//bottom row(i)
    p.row(p.rows()-1)=2-p.row(p.rows()-2);//top row(i)
    p.col(0)=-p.col(1);//left row(i)
    p.col(p.cols()-1)=-p.col(p.cols()-2);//right row(i)
    // v.row(last)=-v.row(last-1);//bottom row(i)
    // cout<<endl<<endl<<u(u.rows()-1,0)<<endl<<endl;
}

