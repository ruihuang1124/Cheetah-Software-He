#include "RobotState.h"

#include <math.h>

#include <eigen3/Eigen/Dense>
#include <iostream>


using std::cout;
using std::endl;

void RobotState::set(flt* p_, flt* v_, flt* q_, flt* w_, flt* r_, flt yaw_) {
  for (u8 i = 0; i < 3; i++) {
    this->p(i) = p_[i];
    this->v(i) = v_[i];
    this->w(i) = w_[i];
  }
  this->q.w() = q_[0];
  this->q.x() = q_[1];
  this->q.y() = q_[2];
  this->q.z() = q_[3];
  this->yaw = yaw_;

  // for(u8 i = 0; i < 12; i++)
  //     this->r_feet(i) = r[i];
  for (u8 rs = 0; rs < 3; rs++)
    for (u8 c = 0; c < 4; c++) this->r_feet(rs, c) = r_[rs * 4 + c];

  R = this->q.toRotationMatrix();
  fpt yc = cos(yaw_);
  fpt ys = sin(yaw_);

  R_yaw << yc, -ys, 0, ys, yc, 0, 0, 0, 1;
  // Spatial inertias from Dynamic
  Matrix<fpt, 3, 1> Id_mini, Id_milab, Id_cheetah3, Id_arcdog, Id_arcdog_mini;
  Id_mini << .07f, 0.26f, 0.242f;
  Id_cheetah3 << 0.41f, 2.1f, 2.1f;
  //    Id_milab << 0.0996f, 0.765f, 0.765f;//25.7kg
  Id_milab << 0.0891f, 0.685f, 0.685f;  // 23kg
  //    Id_milab << 0.1084f, 0.834f, 0.834f;//28kg
  Id_arcdog << 0.029679f, 0.13969f,
      0.15473f;  // TODO. Copy the value from Arcdog.h
                  // 29679.983, -2.510, -267.014, -2.510, 139698.558, 11.971, -267.014, 11.971, 154733.982
  Id_arcdog_mini << 0.0572f, 0.02918f, 0.07133f;
  //  21000, 1000, 2000, 1000, 55000, 0, 2000, 0, 61000;
  // 57200, 160.8, 2.76725, 160.8, 29180, 1050, 2.76725, 1050, 71330;
  I_body_mini.diagonal() = Id_mini;
  I_body_milab.diagonal() = Id_milab;
  I_body_cheetah3.diagonal() = Id_cheetah3;
  I_body_arcdog.diagonal() = Id_arcdog;
  I_body_arcdog_mini.diagonal() = Id_arcdog_mini;
}

void RobotState::print() {
  cout << "Robot State:" << endl
       << "Position\n"
       << p.transpose() << "\nVelocity\n"
       << v.transpose() << "\nAngular Veloctiy\n"
       << w.transpose() << "\nRotation\n"
       << R << "\nYaw Rotation\n"
       << R_yaw << "\nFoot Locations\n"
       << r_feet << "\nInertia\n"
       << I_body << "\nMass\n"
       << m << endl;
}
RobotState::RobotState() {
  I_body.setZero();
  I_body_mini.setZero();
  I_body_milab.setZero();
  I_body_arcdog.setZero();
  I_body_arcdog_mini.setZero();
  I_body_cheetah3.setZero();
  Matrix<fpt, 3, 1> Id_mini, Id_milab, Id_cheetah3, Id_arcdog, Id_arcdog_mini;
  Id_mini << .07f, 0.26f, 0.242f;
  Id_cheetah3 << 0.41f, 2.1f, 2.1f;
  //    Id_milab << 0.0996f, 0.765f, 0.765f;//25.7kg
  Id_milab << 0.0891f, 0.685f, 0.685f;  // 23kg
  //    Id_milab << 0.1084f, 0.834f, 0.834f;//28kg
  Id_arcdog << 0.029679f, 0.13969f,
      0.15473f;  // TODO. Copy the value from Arcdog.h
                 // 29679.983, -2.510, -267.014, -2.510, 139698.558, 11.971, -267.014, 11.971, 154733.982
  Id_arcdog_mini << 0.0572f, 0.02918f, 0.07133f;
  //  21000, 1000, 2000, 1000, 55000, 0, 2000, 0, 61000;
  // 57200, 160.8, 2.76725, 160.8, 29180, 1050, 2.76725, 1050, 71330;
  I_body_mini.diagonal() = Id_mini;
  I_body_milab.diagonal() = Id_milab;
  I_body_cheetah3.diagonal() = Id_cheetah3;
  I_body_arcdog.diagonal() = Id_arcdog;
  I_body_arcdog_mini.diagonal() = Id_arcdog_mini;
}
