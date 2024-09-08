//
// Created by ray on 7/10/24.
//

#ifndef ARCDOG_SOFTWARE_MPCSOLVER_H
#define ARCDOG_SOFTWARE_MPCSOLVER_H
#include <stdio.h>

#include <eigen3/Eigen/Dense>
#include <iostream>

#include "common_types.h"
#include "MPC_interface.h"
#include "MathUtilities.h"

using Eigen::Matrix;
using Eigen::Quaterniond;
using Eigen::Quaternionf;

template <class T>
void vbmpc_print_array(T* array, u16 rows, u16 cols) {
  for (u16 r = 0; r < rows; r++) {
    for (u16 c = 0; c < cols; c++) std::cout << (fpt)array[c + r * cols] << " ";
    printf("\n");
  }
}

template <class T>
void vbmpc_print_named_array(const char* name, T* array, u16 rows, u16 cols) {
  printf("%s:\n", name);
  vbmpc_print_array(array, rows, cols);
}

// print named variable
template <class T>
void vbmpc_pnv(const char* name, T v) {
  printf("%s: ", name);
  std::cout << v << std::endl;
}

template <class T>
T vbmpc_t_min(T a, T b) {
  if (a < b) return a;
  return b;
}

template <class T>
T vbmpc_sq(T a) {
  return a * a;
}

void solve_cmpc(mpc_update_data_t* update, mpc_problem_setup* setup, bool arcdog_mini, bool arcdog);

void solve_vbmpc(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd,
                         Eigen::MatrixXf& Ud, VBMPCParameters& vbmpc_parameters);

void cmpc_quat_to_rpy(Quaternionf q, Matrix<fpt, 3, 1>& rpy);
void vbmpc_ct_ss_mats(Matrix<fpt, 3, 3> I_world, fpt m, Matrix<fpt, 3, 4> r_feet, Matrix<fpt, 3, 3> R_yaw, Matrix<fpt, 13, 13>& A,
                Matrix<fpt, 13, 12>& B);
void cmpc_resize_qp_mats(s16 horizon);
void vbmpc_resize_qp_mats(s16 horizon);
void cmpc_c2qp(Matrix<fpt, 13, 13> Ac, Matrix<fpt, 13, 12> Bc, fpt dt, s16 horizon);
mfp* get_cmpc_q_soln();
mfp* get_vbmpc_q_soln();

#endif  // ARCDOG_SOFTWARE_MPCSOLVER_H
