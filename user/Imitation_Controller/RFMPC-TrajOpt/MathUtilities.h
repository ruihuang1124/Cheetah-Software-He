//
// Created by ray on 7/11/24.
//

#ifndef ARCDOG_SOFTWARE_MATHUTILITIES_H
#define ARCDOG_SOFTWARE_MATHUTILITIES_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

#include <Eigen/Dense>
#include <cmath>
#include <eigen3/unsupported/Eigen/KroneckerProduct>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <vector>
#include "iostream"

using namespace Eigen;

struct VBMPCParameters {
  float Tmpc = 0.2f;
  float mass = 5.5f;
  Matrix3f inertia_body = (Matrix3f() << 0.026f, 0.f, 0.f, 0.f, 0.112f, 0.f, 0.f, 0.f, 0.075f).finished();
  float g = 9.81f;
  float mu = 1.0f;
  Matrix<float, 3, 4> pf34 =
      (Matrix<float, 3, 4>() << 0.15f, 0.15f, -0.15f, -0.15f, 0.094f, -0.094f, 0.094f, -0.094f, 0.f, 0.f, 0.f, 0.f).finished();
  int predHorizon = 10;
  float Umax = 180;
  float decayRate = 1;
  Eigen::Matrix<float, 12, 12> R;
  Eigen::Matrix<float, 12, 12> Q;
  Eigen::Matrix<float, 12, 12> Qf;
};

EXTERNC Matrix<float, 9, 3> fcn_get_N();
EXTERNC Matrix<float, 3, 9> fcn_get_N_inv();
EXTERNC Matrix<float, 9, 3> fcn_get_D(const Vector3f& in);
EXTERNC Matrix<float, 3, 9> fcn_get_F(const Vector3f& k);
EXTERNC Matrix3f hatMap(const Vector3f& a);
EXTERNC void eta_co_xv(const Matrix<float, 12, 1>& force_operation, float dt, float mass, float g, Eigen::Matrix3f& Cx_x,
                       Eigen::Matrix3f& Cx_v, Eigen::Matrix3f& Cv_v, Matrix<float, 3, 12>& Cv_u, Vector3f& Cv_c);
EXTERNC Eigen::VectorXf vec(const Eigen::MatrixXd& m);
EXTERNC void eta_co_R(const Matrix3f& Rop, const Vector3f& wop, float dt, Matrix3f& CE_eta, Matrix3f& CE_w, Vector3f& CE_c);
EXTERNC void eta_co_w(Eigen::Vector3f& xop, const Matrix3f& Rop, const Vector3f& wop, const Matrix<float, 12, 1>& force_operation,
                      float dt, const Eigen::MatrixXf& J, const Eigen::MatrixXf& pf, Matrix3f& CW_x, Matrix3f& CW_eta, Matrix3f& CW_w,
                      Matrix<float, 3, 12>& CW_u, Vector3f& CW_c);
EXTERNC void fcn_get_ABD_eta(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, const VBMPCParameters& p, MatrixXf& A,
                             MatrixXf& B, MatrixXf& D);
EXTERNC void fcn_get_QP_form_eta(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd, Eigen::MatrixXf& Ud,
                                 const VBMPCParameters& p, Eigen::MatrixXf& H, Eigen::VectorXf& g, Eigen::MatrixXf& Aineq,
                                 Eigen::VectorXf& bineq, Eigen::MatrixXf& Aeq, Eigen::VectorXf& beq);
EXTERNC Eigen::Vector3f veeMap(const Eigen::Matrix3f& in);

#endif  // ARCDOG_SOFTWARE_MATHUTILITIES_H
