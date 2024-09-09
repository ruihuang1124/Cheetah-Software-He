//
// Created by ray on 7/10/24.
//
#include "RFMPCSolver.h"

#include <Utilities/Timer.h>
#include <qpSWIFT.h>
#include <stdio.h>
#include <sys/time.h>

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <qpOASES.hpp>

#include "RobotState.h"
#include "common_types.h"
#include "MPC_interface.h"

// #define K_PRINT_EVERYTHING
#define MPC_BIG_NUMBER 5e10

RobotState rs_;
using Eigen::Dynamic;
using std::cout;
using std::endl;

// qpOASES::real_t a;
qp_real* vbmpc_q_soln_;
Matrix<fpt, Dynamic, 13> cmpc_A_qp_;
Matrix<fpt, Dynamic, Dynamic> cmpc_B_qp_;
Matrix<fpt, 13, 12> cmpc_Bdt_;
Matrix<fpt, 13, 13> cmpc_Adt_;
Matrix<fpt, 25, 25> cmpc_ABc_, cmpc_expmm_;
Matrix<fpt, Dynamic, Dynamic> cmpc_S_;
Matrix<fpt, Dynamic, 1> cmpc_X_d_;
Matrix<fpt, Dynamic, 1> cmpc_U_b_;
Matrix<fpt, Dynamic, Dynamic> cmpc_fmat_;

Matrix<fpt, Dynamic, Dynamic> cmpc_qH_;
Matrix<fpt, Dynamic, 1> cmpc_qg_;

Matrix<fpt, Dynamic, Dynamic> cmpc_eye_12h_;

qpOASES::real_t* cmpc_H_qpoases_;
qpOASES::real_t* cmpc_g_qpoases_;
qpOASES::real_t* cmpc_A_qpoases_;
qpOASES::real_t* cmpc_lb_qpoases_;
qpOASES::real_t* cmpc_ub_qpoases_;
qpOASES::real_t* cmpc_q_soln_;

qpOASES::real_t* cmpc_H_red_;
qpOASES::real_t* cmpc_g_red_;
qpOASES::real_t* cmpc_A_red_;
qpOASES::real_t* cmpc_lb_red_;
qpOASES::real_t* cmpc_ub_red_;
qpOASES::real_t* cmpc_q_red_;
u8 cmpc_real_allocated_ = 0;
u8 vbmpc_real_allocated_ = 0;

char cmpc_var_elim_[2000];
char cmpc_con_elim_[2000];

Matrix<fpt, 13, 1> cmpc_x_0_;
Matrix<fpt, 3, 3> cmpc_I_world_;
Matrix<fpt, 13, 13> cmpc_A_ct_;
Matrix<fpt, 13, 12> cmpc_B_ct_r_;

mfp* get_cmpc_q_soln() { return cmpc_q_soln_; }

mfp* get_vbmpc_q_soln() { return vbmpc_q_soln_; }

s8 cmpc_near_zero(fpt a) { return (a < 0.01 && a > -.01); }

s8 cmpc_near_one(fpt a) { return cmpc_near_zero(a - 1); }
void cmpc_matrix_to_real(qpOASES::real_t* dst, Matrix<fpt, Dynamic, Dynamic> src, s16 rows, s16 cols) {
  s32 a = 0;
  for (s16 r = 0; r < rows; r++) {
    for (s16 c = 0; c < cols; c++) {
      dst[a] = src(r, c);
      a++;
    }
  }
}

void cmpc_c2qp(Matrix<fpt, 13, 13> Ac, Matrix<fpt, 13, 12> Bc, fpt dt, s16 horizon) {
  cmpc_ABc_.setZero();
  cmpc_ABc_.block(0, 0, 13, 13) = Ac;
  cmpc_ABc_.block(0, 13, 13, 12) = Bc;
  cmpc_ABc_ = dt * cmpc_ABc_;
  cmpc_expmm_ = cmpc_ABc_.exp();
  cmpc_Adt_ = cmpc_expmm_.block(0, 0, 13, 13);
  cmpc_Bdt_ = cmpc_expmm_.block(0, 13, 13, 12);
#ifdef K_PRINT_EVERYTHING
  cout << "Adt: \n" << Adt << "\nBdt:\n" << Bdt << endl;
#endif
  if (horizon > 30) {  // horizon for gait should match value in here
    throw std::runtime_error("horizon is too long!");
  }

  Matrix<fpt, 13, 13> powerMats[30];  // horizon length should match the value in here
  powerMats[0].setIdentity();
  for (int i = 1; i < horizon + 1; i++) {
    powerMats[i] = cmpc_Adt_ * powerMats[i - 1];
  }

  for (s16 r = 0; r < horizon; r++) {
    cmpc_A_qp_.block(13 * r, 0, 13, 13) = powerMats[r + 1];  // Adt.pow(r+1);
    for (s16 c = 0; c < horizon; c++) {
      if (r >= c) {
        s16 a_num = r - c;
        cmpc_B_qp_.block(13 * r, 12 * c, 13, 12) = powerMats[a_num] /*Adt.pow(a_num)*/ * cmpc_Bdt_;
      }
    }
  }

#ifdef K_PRINT_EVERYTHING
  cout << "AQP:\n" << A_qp << "\nBQP:\n" << B_qp << endl;
#endif
}

void cmpc_resize_qp_mats(s16 horizon) {
  int mcount = 0;
  int h2 = horizon * horizon;

  cmpc_A_qp_.resize(13 * horizon, Eigen::NoChange);
  mcount += 13 * horizon * 1;

  cmpc_B_qp_.resize(13 * horizon, 12 * horizon);
  mcount += 13 * h2 * 12;

  cmpc_S_.resize(13 * horizon, 13 * horizon);
  mcount += 13 * 13 * h2;

  cmpc_X_d_.resize(13 * horizon, Eigen::NoChange);
  mcount += 13 * horizon;

  cmpc_U_b_.resize(20 * horizon, Eigen::NoChange);
  mcount += 20 * horizon;

  cmpc_fmat_.resize(20 * horizon, 12 * horizon);
  mcount += 20 * 12 * h2;

  cmpc_qH_.resize(12 * horizon, 12 * horizon);
  mcount += 12 * 12 * h2;

  cmpc_qg_.resize(12 * horizon, Eigen::NoChange);
  mcount += 12 * horizon;

  cmpc_eye_12h_.resize(12 * horizon, 12 * horizon);
  mcount += 12 * 12 * horizon;

  // printf("realloc'd %d floating point numbers.\n",mcount);
  mcount = 0;

  cmpc_A_qp_.setZero();
  cmpc_B_qp_.setZero();
  cmpc_S_.setZero();
  cmpc_X_d_.setZero();
  cmpc_U_b_.setZero();
  cmpc_fmat_.setZero();
  cmpc_qH_.setZero();
  cmpc_eye_12h_.setIdentity();

  // TODO: use realloc instead of free/malloc on size changes

  if (cmpc_real_allocated_) {
    free(cmpc_H_qpoases_);
    free(cmpc_g_qpoases_);
    free(cmpc_A_qpoases_);
    free(cmpc_lb_qpoases_);
    free(cmpc_ub_qpoases_);
    free(cmpc_q_soln_);
    free(cmpc_H_red_);
    free(cmpc_g_red_);
    free(cmpc_A_red_);
    free(cmpc_lb_red_);
    free(cmpc_ub_red_);
    free(cmpc_q_red_);
  }

  cmpc_H_qpoases_ = (qpOASES::real_t*)malloc(12 * 12 * horizon * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * 12 * h2;
  cmpc_g_qpoases_ = (qpOASES::real_t*)malloc(12 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * horizon;
  cmpc_A_qpoases_ = (qpOASES::real_t*)malloc(12 * 20 * horizon * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * 20 * h2;
  cmpc_lb_qpoases_ = (qpOASES::real_t*)malloc(20 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 20 * horizon;
  cmpc_ub_qpoases_ = (qpOASES::real_t*)malloc(20 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 20 * horizon;
  cmpc_q_soln_ = (qpOASES::real_t*)malloc(12 * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * horizon;

  cmpc_H_red_ = (qpOASES::real_t*)malloc(12 * 12 * horizon * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * 12 * h2;
  cmpc_g_red_ = (qpOASES::real_t*)malloc(12 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * horizon;
  cmpc_A_red_ = (qpOASES::real_t*)malloc(12 * 20 * horizon * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * 20 * h2;
  cmpc_lb_red_ = (qpOASES::real_t*)malloc(20 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 20 * horizon;
  cmpc_ub_red_ = (qpOASES::real_t*)malloc(20 * 1 * horizon * sizeof(qpOASES::real_t));
  mcount += 20 * horizon;
  cmpc_q_red_ = (qpOASES::real_t*)malloc(12 * horizon * sizeof(qpOASES::real_t));
  mcount += 12 * horizon;
  cmpc_real_allocated_ = 1;

  // printf("malloc'd %d floating point numbers.\n",mcount);

#ifdef K_DEBUG
  printf("RESIZED MATRICES FOR HORIZON: %d\n", horizon);
#endif
}

void vbmpc_resize_qp_mats(s16 horizon) {
  // TODO: use realloc instead of free/malloc on size changes
  if (vbmpc_real_allocated_) {
    free(vbmpc_q_soln_);
  }
  vbmpc_q_soln_ = (qp_real*)malloc(24 * horizon * sizeof(qp_real));
  vbmpc_real_allocated_ = 1;
  // printf("malloc'd %d floating point numbers.\n",mcount);
}

inline Matrix<fpt, 3, 3> cmpc_cross_mat(Matrix<fpt, 3, 3> I_inv, Matrix<fpt, 3, 1> r) {
  Matrix<fpt, 3, 3> cm;
  cm << 0.f, -r(2), r(1), r(2), 0.f, -r(0), -r(1), r(0), 0.f;
  return I_inv * cm;
}
// continuous time state space matrices.
void cmpc_ct_ss_mats(Matrix<fpt, 3, 3> I_world, fpt m, Matrix<fpt, 3, 4> r_feet, Matrix<fpt, 3, 3> R_yaw, Matrix<fpt, 13, 13>& A,
                     Matrix<fpt, 13, 12>& B, float x_drag) {
  A.setZero();
  A(3, 9) = 1.f;
  A(11, 9) = x_drag;
  A(4, 10) = 1.f;
  A(5, 11) = 1.f;

  A(11, 12) = 1.f;
  A.block(0, 6, 3, 3) = R_yaw.transpose();

  B.setZero();
  Matrix<fpt, 3, 3> I_inv = I_world.inverse();

  for (s16 b = 0; b < 4; b++) {
    B.block(6, b * 3, 3, 3) = cmpc_cross_mat(I_inv, r_feet.col(b));
    B.block(9, b * 3, 3, 3) = Matrix<fpt, 3, 3>::Identity() / m;
  }
}

void cmpc_quat_to_rpy(Quaternionf q, Matrix<fpt, 3, 1>& rpy) {
  // from my MATLAB implementation

  // edge case!
  fpt as = vbmpc_t_min(-2. * (q.x() * q.z() - q.w() * q.y()), .99999);
  rpy(0) = atan2(2.f * (q.x() * q.y() + q.w() * q.z()), vbmpc_sq(q.w()) + vbmpc_sq(q.x()) - vbmpc_sq(q.y()) - vbmpc_sq(q.z()));
  rpy(1) = asin(as);
  rpy(2) = atan2(2.f * (q.y() * q.z() + q.w() * q.x()), vbmpc_sq(q.w()) - vbmpc_sq(q.x()) - vbmpc_sq(q.y()) + vbmpc_sq(q.z()));
}
void print_cmpc_problem_setup(mpc_problem_setup* setup) {
  printf("DT: %.3f\n", setup->dt);
  printf("Mu: %.3f\n", setup->mu);
  printf("F_Max: %.3f\n", setup->f_max);
  printf("Horizon: %d\n", setup->horizon);
}

void print_cmpc_update_data(mpc_update_data_t* update, s16 horizon) {
  vbmpc_print_named_array("p", update->p, 1, 3);
  vbmpc_print_named_array("v", update->v, 1, 3);
  vbmpc_print_named_array("q", update->q, 1, 4);
  vbmpc_print_named_array("w", update->r, 3, 4);
  vbmpc_pnv("Yaw", update->yaw);
  vbmpc_print_named_array("weights", update->weights, 1, 12);
  vbmpc_print_named_array("trajectory", update->traj, horizon, 12);
  vbmpc_pnv("Alpha", update->alpha);
  vbmpc_print_named_array("gait", update->gait, horizon, 4);
}

void solve_cmpc(mpc_update_data_t* update, mpc_problem_setup* setup, bool arcdog_mini, bool arcdog) {
  if (arcdog_mini) {
    rs_.m = rs_.m_arcdog_mini;
    rs_.I_body = rs_.I_body_arcdog_mini;
  } else if (arcdog) {
    rs_.m = rs_.m_arcdog;
    rs_.I_body = rs_.I_body_arcdog;
  } else {
    rs_.m = rs_.m_mini;
    rs_.I_body = rs_.I_body_mini;
  }
  rs_.set(update->p, update->v, update->q, update->w, update->r, update->yaw);
#ifdef K_PRINT_EVERYTHING

  printf("-----------------\n");
  printf("   PROBLEM DATA  \n");
  printf("-----------------\n");
  print_problem_setup(setup);

  printf("-----------------\n");
  printf("    ROBOT DATA   \n");
  printf("-----------------\n");
  rs.print();
  print_update_data(update, setup->horizon);
#endif

  // roll pitch yaw
  Matrix<fpt, 3, 1> rpy;
  cmpc_quat_to_rpy(rs_.q, rpy);

  // initial state (13 state representation)
  cmpc_x_0_ << rpy(2), rpy(1), rpy(0), rs_.p, rs_.w, rs_.v, -9.8f;
  cmpc_I_world_ = rs_.R_yaw * rs_.I_body * rs_.R_yaw.transpose();  // original
  // I_world = rs.R_yaw.transpose() * rs.I_body * rs.R_yaw;
  // cout<<rs.R_yaw<<endl;
  cmpc_ct_ss_mats(cmpc_I_world_, rs_.m, rs_.r_feet, rs_.R_yaw, cmpc_A_ct_, cmpc_B_ct_r_, update->x_drag);

#ifdef K_PRINT_EVERYTHING
  cout << "Initial state: \n" << x_0 << endl;
  cout << "World Inertia: \n" << I_world << endl;
  cout << "A CT: \n" << A_ct << endl;
  cout << "B CT (simplified): \n" << B_ct_r << endl;
#endif
  // QP matrices
  cmpc_c2qp(cmpc_A_ct_, cmpc_B_ct_r_, setup->dt, setup->horizon);

  // weights
  Matrix<fpt, 13, 1> full_weight;
  for (u8 i = 0; i < 12; i++) full_weight(i) = update->weights[i];
  full_weight(12) = 0.f;
  cmpc_S_.diagonal() = full_weight.replicate(setup->horizon, 1);

  // trajectory
  for (s16 i = 0; i < setup->horizon; i++) {
    for (s16 j = 0; j < 12; j++) cmpc_X_d_(13 * i + j, 0) = update->traj[12 * i + j];
  }
  // cout<<"XD:\n"<<X_d<<endl;

  // note - I'm not doing the shifting here.
  s16 k = 0;
  for (s16 i = 0; i < setup->horizon; i++) {
    for (s16 j = 0; j < 4; j++) {
      cmpc_U_b_(5 * k + 0) = MPC_BIG_NUMBER;
      cmpc_U_b_(5 * k + 1) = MPC_BIG_NUMBER;
      cmpc_U_b_(5 * k + 2) = MPC_BIG_NUMBER;
      cmpc_U_b_(5 * k + 3) = MPC_BIG_NUMBER;
      cmpc_U_b_(5 * k + 4) = update->gait[i * 4 + j] * setup->f_max;
      k++;
    }
  }

  fpt mu = 1.f / setup->mu;
  Matrix<fpt, 5, 3> f_block;

  f_block << mu, 0, 1.f, -mu, 0, 1.f, 0, mu, 1.f, 0, -mu, 1.f, 0, 0, 1.f;

  for (s16 i = 0; i < setup->horizon * 4; i++) {
    cmpc_fmat_.block(i * 5, i * 3, 5, 3) = f_block;
  }

  cmpc_qH_ = 2 * (cmpc_B_qp_.transpose() * cmpc_S_ * cmpc_B_qp_ + update->alpha * cmpc_eye_12h_);
  cmpc_qg_ = 2 * cmpc_B_qp_.transpose() * cmpc_S_ * (cmpc_A_qp_ * cmpc_x_0_ - cmpc_X_d_);

  cmpc_matrix_to_real(cmpc_H_qpoases_, cmpc_qH_, setup->horizon * 12, setup->horizon * 12);
  cmpc_matrix_to_real(cmpc_g_qpoases_, cmpc_qg_, setup->horizon * 12, 1);
  cmpc_matrix_to_real(cmpc_A_qpoases_, cmpc_fmat_, setup->horizon * 20, setup->horizon * 12);
  cmpc_matrix_to_real(cmpc_ub_qpoases_, cmpc_U_b_, setup->horizon * 20, 1);

  for (s16 i = 0; i < 20 * setup->horizon; i++) cmpc_lb_qpoases_[i] = 0.0f;

  s16 num_constraints = 20 * setup->horizon;
  s16 num_variables = 12 * setup->horizon;

  qpOASES::int_t nWSR = 100;

  int new_vars = num_variables;
  int new_cons = num_constraints;

  for (int i = 0; i < num_constraints; i++) cmpc_con_elim_[i] = 0;

  for (int i = 0; i < num_variables; i++) cmpc_var_elim_[i] = 0;

  for (int i = 0; i < num_constraints; i++) {
    if (!(cmpc_near_zero(cmpc_lb_qpoases_[i]) && cmpc_near_zero(cmpc_ub_qpoases_[i]))) continue;
    double* c_row = &cmpc_A_qpoases_[i * num_variables];
    for (int j = 0; j < num_variables; j++) {
      if (cmpc_near_one(c_row[j])) {
        new_vars -= 3;
        new_cons -= 5;
        int cs = (j * 5) / 3 - 3;
        cmpc_var_elim_[j - 2] = 1;
        cmpc_var_elim_[j - 1] = 1;
        cmpc_var_elim_[j] = 1;
        cmpc_con_elim_[cs] = 1;
        cmpc_con_elim_[cs + 1] = 1;
        cmpc_con_elim_[cs + 2] = 1;
        cmpc_con_elim_[cs + 3] = 1;
        cmpc_con_elim_[cs + 4] = 1;
      }
    }
  }
  // if(new_vars != num_variables)
  if (1 == 1) {
    int var_ind[new_vars];
    int con_ind[new_cons];
    int vc = 0;
    for (int i = 0; i < num_variables; i++) {
      if (!cmpc_var_elim_[i]) {
        if (!(vc < new_vars)) {
          printf("BAD ERROR 1\n");
        }
        var_ind[vc] = i;
        vc++;
      }
    }
    vc = 0;
    for (int i = 0; i < num_constraints; i++) {
      if (!cmpc_con_elim_[i]) {
        if (!(vc < new_cons)) {
          printf("BAD ERROR 1\n");
        }
        con_ind[vc] = i;
        vc++;
      }
    }
    for (int i = 0; i < new_vars; i++) {
      int olda = var_ind[i];
      cmpc_g_red_[i] = cmpc_g_qpoases_[olda];
      for (int j = 0; j < new_vars; j++) {
        int oldb = var_ind[j];
        cmpc_H_red_[i * new_vars + j] = cmpc_H_qpoases_[olda * num_variables + oldb];
      }
    }

    for (int con = 0; con < new_cons; con++) {
      for (int st = 0; st < new_vars; st++) {
        float cval = cmpc_A_qpoases_[(num_variables * con_ind[con]) + var_ind[st]];
        cmpc_A_red_[con * new_vars + st] = cval;
      }
    }
    for (int i = 0; i < new_cons; i++) {
      int old = con_ind[i];
      cmpc_ub_red_[i] = cmpc_ub_qpoases_[old];
      cmpc_lb_red_[i] = cmpc_lb_qpoases_[old];
    }

    Timer solve_timer;
    qpOASES::QProblem problem_red(new_vars, new_cons);
    qpOASES::Options op;
    op.setToMPC();
    op.printLevel = qpOASES::PL_NONE;
    problem_red.setOptions(op);
    // int_t nWSR = 50000;

    int rval = problem_red.init(cmpc_H_red_, cmpc_g_red_, cmpc_A_red_, NULL, NULL, cmpc_lb_red_, cmpc_ub_red_, nWSR);
    (void)rval;
    int rval2 = problem_red.getPrimalSolution(cmpc_q_red_);
    if (rval2 != qpOASES::SUCCESSFUL_RETURN) printf("failed to solve!\n");

    // printf("solve time: %.3f ms, size %d, %d\n", solve_timer.getMs(), new_vars, new_cons);

    vc = 0;
    for (int i = 0; i < num_variables; i++) {
      if (cmpc_var_elim_[i]) {
        cmpc_q_soln_[i] = 0.0f;
      } else {
        cmpc_q_soln_[i] = cmpc_q_red_[vc];
        vc++;
      }
    }
  }
#ifdef K_PRINT_EVERYTHING
  // cout<<"fmat:\n"<<fmat<<endl;
#endif
}

void solve_vbmpc(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd, Eigen::MatrixXf& Ud,
                 VBMPCParameters& vbmpc_parameters) {
  Eigen::MatrixXf H, Aeq, Aineq;
  Eigen::VectorXf g, beq, bineq;
  //    Timer t2;
  fcn_get_QP_form_eta(Xt, Ut, Xd, Ud, vbmpc_parameters, H, g, Aineq, bineq, Aeq, beq);
  //    printf("generate QP for vbmpc use %f ms\n", t2.getMs());

  //    std::cout << "The matrix H is:\n" << H << "\n\n";
  //    std::cout << "The matrix g is:\n" << g << "\n\n";
  //    std::cout << "The matrix Aineq is:\n" << Aineq << "\n\n";
  //    std::cout << "The matrix bineq is:\n" << bineq << "\n\n";
  //    std::cout << "The matrix Aeq is:\n" << Aeq << "\n\n";
  //    std::cout << "The matrix beq is:\n" << beq << "\n\n";

  MatrixXd P = H.cast<double>();
  MatrixXd A = Aeq.cast<double>();
  MatrixXd G = Aineq.cast<double>();
  VectorXd c = g.cast<double>();
  VectorXd b = beq.cast<double>();
  VectorXd h = bineq.cast<double>();
  QP* myQP;
  //  std::cout << "number of Decision Variables is:" << P.rows() << "\n";
  //  std::cout << "number of Inequality Constraints is:" << G.rows() << "\n";
  //  std::cout << "number of equality Constraints is:" << A.rows() << "\n";
  myQP = QP_SETUP_dense(P.rows(), G.rows(), A.rows(), P.data(), A.data(), G.data(), c.data(), h.data(), b.data(), NULL,
                        COLUMN_MAJOR_ORDERING);
  myQP->options->maxit = 10000;
  myQP->options->reltol = 1e-3;
  myQP->options->abstol = 1e-3;
  myQP->options->verbose = 0;
  qp_int ExitCode = QP_SOLVE(myQP);
  if (myQP != NULL) printf("Setup Time     : %f ms\n", myQP->stats->tsetup * 1000.0);
  if (ExitCode == QP_OPTIMAL) {
    printf("Solve Time     : %f ms\n", (myQP->stats->tsolve + myQP->stats->tsetup) * 1000.0);
    printf("KKT_Solve Time : %f ms\n", myQP->stats->kkt_time * 1000.0);
    printf("LDL Time       : %f ms\n", myQP->stats->ldl_numeric * 1000.0);
    printf("Diff	       : %f ms\n", (myQP->stats->kkt_time - myQP->stats->ldl_numeric) * 1000.0);
    printf("Iterations     : %ld\n", myQP->stats->IterationCount);
    printf("Optimal Solution Found\n");
  }
  if (ExitCode == QP_MAXIT) {
    printf("Solve Time     : %f ms\n", myQP->stats->tsolve * 1000.0);
    printf("KKT_Solve Time : %f ms\n", myQP->stats->kkt_time * 1000.0);
    printf("LDL Time       : %f ms\n", myQP->stats->ldl_numeric * 1000.0);
    printf("Diff	       : %f ms\n", (myQP->stats->kkt_time - myQP->stats->ldl_numeric) * 1000.0);
    printf("Iterations     : %ld\n", myQP->stats->IterationCount);
    printf("Maximum Iterations reached\n");
  }
  if (ExitCode == QP_FATAL) {
    printf("Unknown Error Detected\n");
  }
  if (ExitCode == QP_KKTFAIL) {
    printf("LDL Factorization fail\n");
  }
  /*! The Solution can be found as real pointer in myQP->x;It is an array of Dimension n*/
  //  std::cout << "First 12 Solution is:" << std::endl;
  //  for (int i = 0; i < 12; ++i) {
  //    std::cout << "x[" << i << "]: " << myQP->x[i] << std::endl;
  //  }
  //  s16 num_constraints = 20 * vbmpc_parameters.predHorizon;
  s16 num_optimized_variables = 24 * vbmpc_parameters.predHorizon;
  for (int i = 0; i < num_optimized_variables; i++) {
    vbmpc_q_soln_[i] = myQP->x[i];
  }
  std::cout << "finish setting vbmpc_q_soln_:" << std::endl;

  QP_CLEANUP_dense(myQP);
}