#pragma once
#include "./CLASSES.h"
//#include <eigen3/Eigen/Dense>
#include "../eigen3/Eigen/Dense"

namespace FIELDVARS{
    inline Eigen::VectorXd x(1), y(1), xc(1), yc(1);
    inline Eigen::VectorXd delx(1), dely(1), dx(1), dy(1); //dels stored at center, ds stores at faces
    inline double delx_min, dely_min, del_min;
    inline Eigen::VectorXd delx_inv(1), dely_inv(1); //Inverse of dels
    inline Eigen::VectorXd dx_inv(1), dy_inv(1); //Inverse of dels
} 
