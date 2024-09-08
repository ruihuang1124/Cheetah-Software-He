//
// Created by ray on 7/10/24.
//

#ifndef ARCDOG_SOFTWARE_MPC_INTERFACE_H
#define ARCDOG_SOFTWARE_MPC_INTERFACE_H

#define MPC_MAX_GAIT_SEGMENTS 72  // horizon for gait should match value in here

// #include "common_types.h"

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#include <pthread.h>
#include <stdio.h>
#include <string.h>

#include <eigen3/Eigen/Dense>

using namespace Eigen;

struct mpc_problem_setup {
  float dt;
  float mu;
  float f_max;
  int horizon;
};

struct mpc_update_data_t {
  float p[3];
  float v[3];
  float q[4];
  float w[3];
  float r[12];
  float yaw;
  float weights[12];
  float traj[12 * MPC_MAX_GAIT_SEGMENTS];
  float alpha;
  unsigned char gait[MPC_MAX_GAIT_SEGMENTS];
  unsigned char hack_pad[1000];
  int max_iterations;
  double rho, sigma, solver_alpha, terminate;
  float x_drag;
};

EXTERNC void mpc_setup_problem(double dt, int horizon, double mu, double f_max);
EXTERNC double get_cmpc_solution(int index);
EXTERNC double get_vbmpc_solution(int index);
EXTERNC void update_cmpc_solver_settings(int max_iter, double rho, double sigma, double solver_alpha, double terminate);
EXTERNC void update_vbmpc_solver_settings(int max_iter, double rho, double sigma, double solver_alpha, double terminate);
EXTERNC void update_cmpc_problem_data_floats(float* p, float* v, float* q, float* w, float* r, float yaw, float* weights,
                                              float* state_trajectory, float alpha, int* gait, bool arcdog_mini, bool arcdog);
EXTERNC void update_vbmpc_problem_data_float(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd,
                                       Eigen::MatrixXf& Ud, float* weights, float alpha, int* gait, bool arcdog_mini,
                                       bool arcdog);

EXTERNC void update_vbmpc_problem_data_float_new(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd,
                                             Eigen::MatrixXf& Ud, float* weights_Q, float* weights_Qf, float* weights_R, float alpha, int* gait, bool arcdog_mini,
                                             bool arcdog);
EXTERNC void update_vbmpc_problem_data_float_latest(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd,
                                                 Eigen::MatrixXf& Ud, float* weights_Q, float* weights_Qf, float* weights_R, float alpha);

#endif  // ARCDOG_SOFTWARE_MPC_INTERFACE_H
