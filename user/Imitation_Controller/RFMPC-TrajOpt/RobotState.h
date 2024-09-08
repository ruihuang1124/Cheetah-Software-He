#ifndef _RobotState
#define _RobotState

#include <eigen3/Eigen/Dense>

#include "common_types.h"

using Eigen::Matrix;
using Eigen::Quaternionf;

#include "common_types.h"
class RobotState {
 public:
  RobotState();
  void set(flt* p, flt* v, flt* q, flt* w, flt* r, flt yaw);
  // void compute_rotations();
  void print();
  Matrix<fpt, 3, 1> p, v, w;
  Matrix<fpt, 3, 4> r_feet;
  Matrix<fpt, 3, 3> R;
  Matrix<fpt, 3, 3> R_yaw;
  Matrix<fpt, 3, 3> I_body_mini, I_body_milab, I_body_cheetah3, I_body, I_body_arcdog, I_body_arcdog_mini;
  Quaternionf q;
  fpt yaw;
  // Total mass of Quadruped
  fpt m_mini = 9;
  fpt m_milab = 23;      // 25.7;
  fpt m_arcdog = 8.098;  // TODO. Copy the value from Arcdog.h
  fpt m_arcdog_mini = 6.847;
  fpt m_cheetah3 = 43;
  fpt m;
  // fpt m = 50.236; //DH
  // private:
};
#endif
