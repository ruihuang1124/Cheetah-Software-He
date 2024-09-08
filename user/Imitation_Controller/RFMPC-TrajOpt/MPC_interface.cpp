//
// Created by ray on 7/10/24.
//
#include "MPC_interface.h"

#include "RobotState.h"
#include "common_types.h"
#include "RFMPCSolver.h"

#define K_NUM_LEGS 4

mpc_problem_setup mpc_problem_configuration_;
VBMPCParameters vbmpc_parameters_;
pthread_mutex_t problem_cfg_mt_;
pthread_mutex_t update_mt_;
mpc_update_data_t mpc_update_data_;
pthread_t solve_thread_;
u8 first_run_ = 1;
int cmpc_has_solved_ = 0;
int vbmpc_has_solved_ = 0;
RobotState robot_state_;

// inline to motivate gcc to unroll the loop in here.
inline void vbmpc_mfp_to_flt(flt* dst, mfp* src, s32 n_items) {
  for (s32 i = 0; i < n_items; i++) *dst++ = *src++;
}

inline void vbmpc_mint_to_u8(u8* dst, mint* src, s32 n_items) {
  for (s32 i = 0; i < n_items; i++) *dst++ = *src++;
}

void initialize_mpc_thread() {
  // printf("Initializing MPC!\n");
  if (pthread_mutex_init(&problem_cfg_mt_, NULL) != 0) printf("[MPC ERROR] Failed to initialize problem configuration mutex.\n");

  if (pthread_mutex_init(&update_mt_, NULL) != 0) printf("[MPC ERROR] Failed to initialize update data mutex.\n");
}

void mpc_setup_problem(double dt, int horizon, double mu, double f_max) {
  // mu = 0.6;
  if (first_run_) {
    first_run_ = false;
    initialize_mpc_thread();
  }
  // pthread_mutex_lock(&problem_cfg_mt);

  mpc_problem_configuration_.horizon = horizon;
  mpc_problem_configuration_.f_max = f_max;
  mpc_problem_configuration_.mu = mu;
  mpc_problem_configuration_.dt = dt;

  // pthread_mutex_unlock(&problem_cfg_mt);
  cmpc_resize_qp_mats(horizon);
  vbmpc_resize_qp_mats(horizon);
}

void update_cmpc_solver_settings(int max_iter, double rho, double sigma, double solver_alpha, double terminate) {
  mpc_update_data_.max_iterations = max_iter;
  mpc_update_data_.rho = rho;
  mpc_update_data_.sigma = sigma;
  mpc_update_data_.solver_alpha = solver_alpha;
  mpc_update_data_.terminate = terminate;
}

void update_vbmpc_solver_settings(int max_iter, double rho, double sigma, double solver_alpha, double terminate) {
  mpc_update_data_.max_iterations = max_iter;
  mpc_update_data_.rho = rho;
  mpc_update_data_.sigma = sigma;
  mpc_update_data_.solver_alpha = solver_alpha;
  mpc_update_data_.terminate = terminate;
}

void update_cmpc_problem_data_floats(float* p, float* v, float* q, float* w, float* r, float yaw, float* weights,
                                      float* state_trajectory, float alpha, int* gait, bool arcdog_mini, bool arcdog) {
  mpc_update_data_.alpha = alpha;
  mpc_update_data_.yaw = yaw;
  vbmpc_mint_to_u8(mpc_update_data_.gait, gait, 4 * mpc_problem_configuration_.horizon);
  memcpy((void*)mpc_update_data_.p, (void*)p, sizeof(float) * 3);
  memcpy((void*)mpc_update_data_.v, (void*)v, sizeof(float) * 3);
  memcpy((void*)mpc_update_data_.q, (void*)q, sizeof(float) * 4);
  memcpy((void*)mpc_update_data_.w, (void*)w, sizeof(float) * 3);
  memcpy((void*)mpc_update_data_.r, (void*)r, sizeof(float) * 12);
  memcpy((void*)mpc_update_data_.weights, (void*)weights, sizeof(float) * 12);
  memcpy((void*)mpc_update_data_.traj, (void*)state_trajectory, sizeof(float) * 12 * mpc_problem_configuration_.horizon);
  solve_cmpc(&mpc_update_data_, &mpc_problem_configuration_, arcdog_mini, arcdog);  // TODO, use the corresponding solver for vbmpc.
  cmpc_has_solved_ = 1;
}

void update_vbmpc_problem_data_float(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd,
                               Eigen::MatrixXf& Ud, float* weights, float alpha, int* gait, bool arcdog_mini, bool arcdog) {
  memcpy((void*)mpc_update_data_.weights, (void*)weights, sizeof(float) * 12);
  vbmpc_mint_to_u8(mpc_update_data_.gait, gait, 4 * mpc_problem_configuration_.horizon);
  mpc_update_data_.alpha = alpha;
  vbmpc_parameters_.Tmpc = mpc_problem_configuration_.dt;
  vbmpc_parameters_.mu = mpc_problem_configuration_.mu;
  vbmpc_parameters_.Umax = mpc_problem_configuration_.f_max;
  vbmpc_parameters_.predHorizon = mpc_problem_configuration_.horizon;
  vbmpc_parameters_.Q.setZero();
  vbmpc_parameters_.R.setZero();
  vbmpc_parameters_.Qf.setZero();
  for (int i = 0; i < 12; ++i) {
    vbmpc_parameters_.Q(i, i) = mpc_update_data_.weights[i];
    vbmpc_parameters_.Qf(i,i) = weights[i];
    vbmpc_parameters_.R(i, i) = mpc_update_data_.alpha;
  }
  if (arcdog_mini) {
    robot_state_.m = robot_state_.m_arcdog_mini;
    robot_state_.I_body = robot_state_.I_body_arcdog_mini;
  } else if (arcdog) {
    robot_state_.m = robot_state_.m_arcdog;
    robot_state_.I_body = robot_state_.I_body_arcdog;
  } else {
    robot_state_.m = robot_state_.m_mini;
    robot_state_.I_body = robot_state_.I_body_mini;
  }
  vbmpc_parameters_.mass = robot_state_.m;
  vbmpc_parameters_.inertia_body = robot_state_.I_body;
  solve_vbmpc(Xt, Ut, Xd, Ud, vbmpc_parameters_);
  vbmpc_has_solved_ = 1;
}


void update_vbmpc_problem_data_float_new(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd,
                                         Eigen::MatrixXf& Ud, float* weights_Q, float* weights_Qf, float* weights_R, float alpha, int* gait, bool arcdog_mini,
                                         bool arcdog) {
  memcpy((void*)mpc_update_data_.weights, (void*)weights_Q, sizeof(float) * 12);
  vbmpc_mint_to_u8(mpc_update_data_.gait, gait, 4 * mpc_problem_configuration_.horizon);
  mpc_update_data_.alpha = alpha;
  vbmpc_parameters_.Tmpc = mpc_problem_configuration_.dt;
  vbmpc_parameters_.mu = mpc_problem_configuration_.mu;
  vbmpc_parameters_.Umax = mpc_problem_configuration_.f_max;
  vbmpc_parameters_.predHorizon = mpc_problem_configuration_.horizon;
  vbmpc_parameters_.Q.setZero();
  vbmpc_parameters_.R.setZero();
  vbmpc_parameters_.Qf.setZero();
  for (int i = 0; i < 12; ++i) {
    vbmpc_parameters_.Q(i, i) = weights_Q[i];
    vbmpc_parameters_.Qf(i,i) = weights_Qf[i];
    vbmpc_parameters_.R(i, i) = weights_R[i];
  }
  if (arcdog_mini) {
    robot_state_.m = robot_state_.m_arcdog_mini;
    robot_state_.I_body = robot_state_.I_body_arcdog_mini;
  } else if (arcdog) {
    robot_state_.m = robot_state_.m_arcdog;
    robot_state_.I_body = robot_state_.I_body_arcdog;
  } else {
    robot_state_.m = robot_state_.m_mini;
    robot_state_.I_body = robot_state_.I_body_mini;
  }
  vbmpc_parameters_.mass = robot_state_.m;
  vbmpc_parameters_.inertia_body = robot_state_.I_body;
  solve_vbmpc(Xt, Ut, Xd, Ud, vbmpc_parameters_);
  vbmpc_has_solved_ = 1;
}


void update_vbmpc_problem_data_float_latest(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd,
                                         Eigen::MatrixXf& Ud, float* weights_Q, float* weights_Qf, float* weights_R, float alpha) {
    memcpy((void*)mpc_update_data_.weights, (void*)weights_Q, sizeof(float) * 12);
    mpc_update_data_.alpha = alpha;
    vbmpc_parameters_.Tmpc = mpc_problem_configuration_.dt;
    vbmpc_parameters_.mu = mpc_problem_configuration_.mu;
    vbmpc_parameters_.Umax = mpc_problem_configuration_.f_max;
    vbmpc_parameters_.predHorizon = mpc_problem_configuration_.horizon;
    vbmpc_parameters_.Q.setZero();
    vbmpc_parameters_.R.setZero();
    vbmpc_parameters_.Qf.setZero();
    for (int i = 0; i < 12; ++i) {
        vbmpc_parameters_.Q(i, i) = weights_Q[i];
        vbmpc_parameters_.Qf(i,i) = weights_Qf[i];
        vbmpc_parameters_.R(i, i) = weights_R[i];
    }

    robot_state_.m = robot_state_.m_mini;
    robot_state_.I_body = robot_state_.I_body_mini;
    vbmpc_parameters_.mass = robot_state_.m;
    vbmpc_parameters_.inertia_body = robot_state_.I_body;
    solve_vbmpc(Xt, Ut, Xd, Ud, vbmpc_parameters_);
    vbmpc_has_solved_ = 1;
}

double get_cmpc_solution(int index) {
  if (!cmpc_has_solved_) return 0.f;
  mfp* qs = get_cmpc_q_soln();
  return qs[index];
}
double get_vbmpc_solution(int index) {
  if (!vbmpc_has_solved_) return 0.f;
  mfp* qs = get_vbmpc_q_soln();
  return qs[index];
}
