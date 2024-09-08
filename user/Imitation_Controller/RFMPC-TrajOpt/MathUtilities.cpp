//
// Created by ray on 7/11/24.
//
#include "MathUtilities.h"
Matrix<float, 9, 3> fcn_get_N() {
  Matrix<float, 9, 3> N;
  N << 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0;
  return N;
}

Matrix<float, 3, 9> fcn_get_N_inv() {
  Matrix<float, 3, 9> N_inv;
  N_inv << -0, -0, -0, -0, -0, 0.5, -0, -0.5, -0, -0, -0, -0.5, -0, -0, -0, 0.5, -0, -0, -0, 0.5, 0, -0.5, -0, -0, -0, -0, -0;
  return N_inv;
}

Matrix<float, 9, 3> fcn_get_D(const Vector3f& in) {
  float d = in(0);
  float e = in(1);
  float f = in(2);
  Matrix<float, 9, 3> D;
  D << 0, 0, 0, e, -d, 0, f, 0, -d, -e, d, 0, 0, 0, 0, 0, f, -e, -f, 0, d, 0, -f, e, 0, 0, 0;
  return D;
}

Matrix<float, 3, 9> fcn_get_F(const Vector3f& k) {
  Matrix<float, 3, 9> F;
  F << k(0), k(1), k(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, k(0), k(1), k(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, k(0), k(1), k(2);
  return F;
}

Matrix3f hatMap(const Vector3f& a) {
  Matrix3f H;
  H << 0, -a(2), a(1), a(2), 0, -a(0), -a(1), a(0), 0;
  return H;
}

Eigen::Vector3f veeMap(const Eigen::Matrix3f& in) {
  Eigen::Vector3f out;
  out(0) = -in(1, 2);
  out(1) = in(0, 2);
  out(2) = -in(0, 1);
  return out;
}

void eta_co_xv(const Matrix<float, 12, 1>& force_operation, float dt, float mass, float g, Eigen::Matrix3f& Cx_x,
               Eigen::Matrix3f& Cx_v, Eigen::Matrix3f& Cv_v, Matrix<float, 3, 12>& Cv_u, Vector3f& Cv_c) {
  // Initialize identity matrices
  Cx_x = Eigen::Matrix3f::Identity();
  Cx_v = Eigen::Matrix3f::Identity() * dt;
  Cv_v = Eigen::Matrix3f::Identity();

  // Calculate Cv_u
  Eigen::Matrix<float, 3, 12> temp = Eigen::Matrix<float, 3, 12>::Zero();
  temp << Eigen::Matrix3f::Identity(), Eigen::Matrix3f::Identity(), Eigen::Matrix3f::Identity(), Eigen::Matrix3f::Identity();
  Cv_u = (temp / mass) * dt;

  // Calculate Cv_c
  Cv_c = Cv_u * force_operation + Eigen::Vector3f(0, 0, -g) * dt;
  //  std::cout << "The matrix Cx_x is:\n" << Cx_x << "\n\n";
  //  std::cout << "The matrix Cx_v is:\n" << Cx_v << "\n\n";
  //  std::cout << "The matrix Cv_v is:\n" << Cv_v << "\n\n";
  //  std::cout << "The matrix Cv_u is:\n" << Cv_u << "\n\n";
  //  std::cout << "The matrix Cv_c is:\n" << Cv_c << "\n\n";
}

Eigen::VectorXf vboc_vec(const Eigen::MatrixXf& m) { return Eigen::Map<const Eigen::VectorXf>(m.data(), m.rows() * m.cols()); }
void eta_co_R(const Matrix3f& Rop, const Vector3f& wop, float dt, Matrix3f& CE_eta, Matrix3f& CE_w, Vector3f& CE_c) {
  Matrix<float, 9, 3> N = fcn_get_N();
  Matrix<float, 9, 3> D = fcn_get_D(wop);
  Matrix<float, 9, 3> C_eta = (Eigen::kroneckerProduct(Matrix3f::Identity(), Rop * hatMap(wop)) * N +
                               Eigen::kroneckerProduct(Matrix3f::Identity(), Rop) * D)
                                  .eval();
  Matrix<float, 9, 3> C_w = Eigen::kroneckerProduct(Matrix3f::Identity(), Rop) * N;
  VectorXf C_c = vboc_vec(Rop * hatMap(wop)) - Eigen::kroneckerProduct(Matrix3f::Identity(), Rop) * N * wop;
  MatrixXf invN = fcn_get_N_inv();
  Matrix<float, 3, 9> invN_kron_RopT = invN * Eigen::kroneckerProduct(Matrix3f::Identity(), Rop.transpose());

  CE_eta = Matrix3f::Identity() + invN * dt * Eigen::kroneckerProduct(Matrix3f::Identity(), Rop.transpose()) * C_eta;
  CE_w = invN_kron_RopT * dt * C_w;
  CE_c = invN_kron_RopT * dt * C_c;
  //    std::cout << "The matrix N is:\n" << N << "\n\n";
  //    std::cout << "The matrix invN is:\n" << invN << "\n\n";
  //  std::cout << "The matrix CE_eta is:\n" << CE_eta << "\n\n";
  //  std::cout << "The matrix CE_w is:\n" << CE_w << "\n\n";
  //  std::cout << "The matrix CE_c is:\n" << CE_c << "\n\n";
}

void eta_co_w(Eigen::Vector3f& xop, const Matrix3f& Rop, const Vector3f& wop, const Matrix<float, 12, 1>& force_operation,
              float dt, const Eigen::MatrixXf& J, const Eigen::MatrixXf& pf, Matrix3f& CW_x, Matrix3f& CW_eta, Matrix3f& CW_w,
              Matrix<float, 3, 12>& CW_u, Vector3f& CW_c) {
  Matrix<float, 9, 3> N = fcn_get_N();

  Eigen::Vector3f r1 = pf.col(0) - xop;
  Eigen::Vector3f r2 = pf.col(1) - xop;
  Eigen::Vector3f r3 = pf.col(2) - xop;
  Eigen::Vector3f r4 = pf.col(3) - xop;
  Eigen::Matrix<float, 3, 12> temp1 = Eigen::Matrix<float, 3, 12>::Zero();
  Eigen::Matrix<float, 3, 12> temp2 = Eigen::Matrix<float, 3, 12>::Zero();
  temp1 << hatMap(r1), hatMap(r2), hatMap(r3), hatMap(r4);
  //    std::cout << "The matrix hatMap(r1) is:\n" << hatMap(r1) << "\n\n";
  //    std::cout << "The matrix temp is:\n" << temp << "\n\n";
  Eigen::MatrixXf Mop = temp1 * force_operation;
  Eigen::MatrixXf temp_J_w = hatMap(J * wop) - hatMap(wop) * J;
  temp2 << Eigen::Matrix3f::Identity(), Eigen::Matrix3f::Identity(), Eigen::Matrix3f::Identity(), Eigen::Matrix3f::Identity();
  //    std::cout << "The matrix temp is:\n" << temp << "\n\n";
  Eigen::MatrixXf sum_fop = temp2 * force_operation;
  Eigen::MatrixXf Cx = Rop.transpose() * hatMap(sum_fop);
  Eigen::MatrixXf Ceta = fcn_get_F(Rop.transpose() * Mop) * N - temp_J_w * hatMap(wop);
  Eigen::MatrixXf Cw = temp_J_w;
  Eigen::MatrixXf Cu = Rop.transpose() * temp1;
  Eigen::MatrixXf Cc = -hatMap(wop) * J * wop + Rop.transpose() * Mop - temp_J_w * wop - Cx * xop;
  //    std::cout << "The matrix Cx is:\n" << Cx << "\n\n";
  //    std::cout << "The matrix Ceta is:\n" << Ceta << "\n\n";
  //    std::cout << "The matrix Cw is:\n" << Cw << "\n\n";
  //    std::cout << "The matrix Cu is:\n" << Cu << "\n\n";
  //    std::cout << "The matrix Cc is:\n" << Cc << "\n\n";

  //
  //    // Placeholder for calculating coefficients
  //    // Replace with actual calculations based on the system dynamics
  CW_x = dt * (J.inverse() * Cx);
  CW_eta = dt * (J.inverse() * Ceta);
  CW_w = dt * (J.inverse() * Cw) + Eigen::Matrix3f::Identity();
  CW_u = dt * (J.inverse() * Cu);
  CW_c = dt * (J.inverse() * Cc);

  //    std::cout << "The matrix CW_x is:\n" << CW_x << "\n\n";
  //    std::cout << "The matrix CW_eta is:\n" << CW_eta << "\n\n";
  //    std::cout << "The matrix CW_w is:\n" << CW_w << "\n\n";
  //    std::cout << "The matrix CW_u is:\n" << CW_u << "\n\n";
  //    std::cout << "The matrix CW_c is:\n" << CW_c << "\n\n";
}

void fcn_get_ABD_eta(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, const VBMPCParameters& p, MatrixXf& A,
                     MatrixXf& B, MatrixXf& D) {
  float dt = p.Tmpc;
  Vector3f xop = Xt.block<3, 1>(0, 0);
  Vector3f vop = Xt.block<3, 1>(3, 0);
  Eigen::MatrixXf Rt(3, 3);
  Rt << Xt(6), Xt(7), Xt(8), Xt(9), Xt(10), Xt(11), Xt(12), Xt(13), Xt(14);
  Rt.transposeInPlace();
  Eigen::MatrixXf Rop = Rt;
  Vector3f wop = Xt.block<3, 1>(15, 0);
  Eigen::MatrixXf pf34(3, 4);
  pf34 << Xt(18), Xt(21), Xt(24), Xt(27), Xt(19), Xt(22), Xt(25), Xt(28), Xt(20), Xt(23), Xt(26), Xt(29);
  Matrix<float, 3, 3> Cx_x = Matrix<float, 3, 3>::Zero();
  Matrix<float, 3, 3> Cx_v = Matrix<float, 3, 3>::Zero();
  Matrix<float, 3, 3> Cv_v = Matrix<float, 3, 3>::Zero();
  Matrix<float, 3, 12> Cv_u = Matrix<float, 3, 12>::Zero();
  Vector3f Cv_c = Vector3f::Zero();
  Matrix<float, 3, 3> CE_eta = Matrix<float, 3, 3>::Zero();
  Matrix<float, 3, 3> CE_w = Matrix<float, 3, 3>::Zero();
  Vector3f CE_c = Vector3f::Zero();
  //    Matrix<float, 3, 3> Cw_x, Cw_eta, Cw_w;
  Matrix<float, 3, 3> Cw_x = Matrix<float, 3, 3>::Zero();
  Matrix<float, 3, 3> Cw_eta = Matrix<float, 3, 3>::Zero();
  Matrix<float, 3, 3> Cw_w = Matrix<float, 3, 3>::Zero();
  //    Matrix<float, 3, 12> Cw_u;
  Matrix<float, 3, 12> Cw_u = Matrix<float, 3, 12>::Zero();
  //    Vector3f Cw_c;
  Vector3f Cw_c = Vector3f::Zero();
  eta_co_xv(Ut, dt, p.mass, p.g, Cx_x, Cx_v, Cv_v, Cv_u, Cv_c);
  eta_co_R(Rop, wop, dt, CE_eta, CE_w, CE_c);
  eta_co_w(xop, Rop, wop, Ut, dt, p.inertia_body, pf34, Cw_x, Cw_eta, Cw_w, Cw_u, Cw_c);
  A =  MatrixXf(12,12).setZero();
  B = MatrixXf(12,12).setZero();
  D = MatrixXf(12,1).setZero();
  A << Cx_x, Cx_v, MatrixXf::Zero(3, 6), MatrixXf::Zero(3, 3), Cv_v, MatrixXf::Zero(3, 6),
      MatrixXf::Zero(3, 6), CE_eta, CE_w, Cw_x, MatrixXf::Zero(3, 3), Cw_eta, Cw_w;
  B  << MatrixXf::Zero(3, 12), Cv_u, MatrixXf::Zero(3, 12), Cw_u;
  D  << Vector3f::Zero(), Cv_c, CE_c, Cw_c;
  //  std::cout << "The matrix A is:\n" << A << "\n\n";
  //  std::cout << "The matrix B is:\n" << B << "\n\n";
  //  std::cout << "The matrix D is:\n" << D << "\n\n";
}

void fcn_get_QP_form_eta(const Matrix<float, 30, 1>& Xt, const Matrix<float, 12, 1>& Ut, Eigen::MatrixXf& Xd, Eigen::MatrixXf& Ud,
                         const VBMPCParameters& p, Eigen::MatrixXf& H, Eigen::VectorXf& g, Eigen::MatrixXf& Aineq,
                         Eigen::VectorXf& bineq, Eigen::MatrixXf& Aeq, Eigen::VectorXf& beq) {
  // min. 0.5 * x' * H *x + g' * x
  // s.t. Aineq *x <= bineq
  // Aeq * x <= beq
  // X = [pc dpc vR wb pf]': [30,1]
  // q = [pc dpc eta wb]: [12 1]
  // lb/ub - [4,n_hor]

  float mu = p.mu;
  int n_hor = p.predHorizon;
  float Umax = p.Umax;
  float decayRate = p.decayRate;

  Eigen::MatrixXf R = p.R;
  Eigen::MatrixXf Q = p.Q;
  Eigen::MatrixXf Qf = p.Qf;
  //    std::cout << "Qf is: \n" << Qf << "\n";


  Eigen::MatrixXf Qx = Q.block<3, 3>(0, 0);
  Eigen::MatrixXf Qv = Q.block<3, 3>(3, 3);
  Eigen::MatrixXf Qeta = Q.block<3, 3>(6, 6);
  Eigen::MatrixXf Qw = Q.block<3, 3>(9, 9);

  Eigen::MatrixXf Qxf = Qf.block<3, 3>(0, 0);
  Eigen::MatrixXf Qvf = Qf.block<3, 3>(3, 3);
  Eigen::MatrixXf Qetaf = Qf.block<3, 3>(6, 6);
  Eigen::MatrixXf Qwf = Qf.block<3, 3>(9, 9);
  //    std::cout << "Qxf is: \n" << Qxf << "\n";

  int nX = 12;
  int nU = 12;

  // A, B, d matrices for linear dynamics
  Eigen::MatrixXf A, B, d;
  // Assuming fcn_get_ABD_eta is implemented elsewhere
  fcn_get_ABD_eta(Xt, Ut, p, A, B, d);

  // Decompose
  Eigen::MatrixXf Rt(3, 3);
  Rt << Xt(6), Xt(7), Xt(8), Xt(9), Xt(10), Xt(11), Xt(12), Xt(13), Xt(14);
  Rt.transposeInPlace();
  //    Rt << Xt.block<9, 1>(6, 0);
  //    std::cout << "Rt is: \n" << Rt << "\n";
  //    std::cout << "Xt.block<9, 1>(6, 0) is: \n" << Xt.block<9, 1>(6, 0) << "\n";
  Eigen::VectorXf qt(nX);
  qt << Xt.topRows(6), 0, 0, 0, Xt(15), Xt(16), Xt(17);
  //    std::cout << "qt is: \n" << qt << "\n";
  Eigen::MatrixXf Fzd(4, Ud.cols());      // create a matrix of size 4 x Ud.cols() to hold the
  std::vector<int> rows = {2, 5, 8, 11};  // indices of the rows to be selected
  for (int i = 0; i < 4; i++) {
    Fzd.row(i) = Ud.row(rows[i]);
  }
  Eigen::MatrixXf lb = -1 * Fzd;
  Eigen::MatrixXf ub = 2 * Fzd;
  //    std::cout << "lb is: \n" << lb << "\n";
  //    std::cout << "ub is: \n" << ub << "\n";

  // Initialize matrices for QP
  H.resize((nX + nU) * n_hor, (nX + nU) * n_hor);
  g.resize(H.rows());
  Aeq.resize(nX * n_hor, (nX + nU) * n_hor);
  beq.resize(Aeq.rows());
  int nAineq_unit = 6;
  Eigen::MatrixXf Aineq_unit(nAineq_unit, 3);
  Aineq_unit << 1, 0, -mu, -1, 0, -mu, 0, 1, -mu, 0, -1, -mu, 0, 0, 1, 0, 0, -1;
  Aineq.resize(4 * nAineq_unit * n_hor, (nX + nU) * n_hor);
  bineq = Eigen::VectorXf::Zero(Aineq.rows());
  H.setZero();
  g.setZero();
  Aineq.setZero();
  bineq.setZero();
  Aeq.setZero();
  beq.setZero();

  // Loop over the prediction horizon
  for (int i_hor = 0; i_hor < n_hor; ++i_hor) {
    Eigen::Vector3f xd = Xd.block<3, 1>(0, i_hor);
    Eigen::Vector3f vd = Xd.block<3, 1>(3, i_hor);
    Eigen::Matrix3f Rd;
    Rd << Xd(6, i_hor), Xd(7, i_hor), Xd(8, i_hor), Xd(9, i_hor), Xd(10, i_hor), Xd(11, i_hor), Xd(12, i_hor), Xd(13, i_hor),
        Xd(14, i_hor);
    Rd.transposeInPlace();
    Eigen::Vector3f wd = Xd.block<3, 1>(15, i_hor);

    //    std::cout << "xd is: \n" << xd << "\n";
    //    std::cout << "vd is: \n" << vd << "\n";
    //        std::cout << "Rd is: \n" << Rd << "\n";
    //    std::cout << "wd is: \n" << wd << "\n";

    //        std::cout << "a is: \n" << veeMap((Rd.transpose() * Rt).log()) << "\n";
    //        std::cout << "a is: \n" << (Rd.transpose() * Rt) << "\n";
    //        std::cout << "a1 is: \n" << (Rd.transpose() * Rt).log() << "\n";
    //        std::cout << "a2 is: \n" << veeMap((Rd.transpose() * Rt).log()) << "\n";

    // Objective function
    int idx_u = (i_hor) * (nX + nU);
    int idx_x = (i_hor) * (nX + nU) + nU;

    if (i_hor == n_hor - 1) {
      H.block(idx_x, idx_x, nX, nX) = Qf * std::pow(decayRate, i_hor);
      g.segment(idx_x, nX) << -Qxf * xd, -Qvf * vd, Qetaf * veeMap((Rd.transpose() * Rt).log()), -Qwf * wd;
      g.segment(idx_x, nX) = g.segment(idx_x, nX) * std::pow(decayRate, i_hor - 1);
      //            std::cout << "g final is: \n" << g.segment(idx_x, nX) << "\n";
    } else {
      H.block(idx_x, idx_x, nX, nX) = Q * std::pow(decayRate, i_hor);
      g.segment(idx_x, nX) << -Qx * xd, -Qv * vd, Qeta * veeMap((Rd.transpose() * Rt).log()), -Qw * wd;
      g.segment(idx_x, nX) = g.segment(idx_x, nX) * std::pow(decayRate, i_hor);
      //            std::cout << "b is: \n" << Qeta * veeMap((Rd.transpose() * Rt).log()) << "\n";
    }
    //        std::cout << "H.block(idx_x, idx_x, nX, nX) : " << i_hor << " is: \n"<<H.block(idx_x, idx_x, nX, nX)<<"\n";
    //        std::cout << "g.segment(idx_x, nX) : " << i_hor << " is: \n" << g.segment(idx_x, nX) << "\n";
    H.block(idx_u, idx_u, nU, nU) = R * std::pow(decayRate, i_hor);
    g.segment(idx_u, nU) = R.transpose() * (Ut - Ud.col(i_hor)) * std::pow(decayRate, i_hor);
    // Ut only current control input vector.
    //        std::cout << "The H.block(idx_u, idx_u, nU, nU) : " << i_hor << " is:
    //        \n"<<H.block(idx_u, idx_u, nU, nU)<<"\n"; std::cout << "The g.segment(idx_u, nU) :
    //        " << i_hor << " is: \n"<<g.segment(idx_u, nU)<<"\n";

    // Equality constraints
    if (i_hor == 0) {
      Aeq.block(0, 0, nX, nU + nX) << -B, Eigen::MatrixXf::Identity(nX, nX);
      beq.head(nX) = A * qt + d;
    } else {
      Aeq.block((i_hor)*nX, (i_hor - 1) * (nX + nU) + nU, nX, 2 * nX + nU) << -A, -B, Eigen::MatrixXf::Identity(nX, nX);
      beq.segment((i_hor)*nX, nX) = d;
    }

    // Inequality constraints
    Eigen::Matrix<float, 24, 12> Fi;
    Fi.setZero();
    Eigen::VectorXf hi = Eigen::VectorXf::Zero(4 * nAineq_unit);
    for (int i_leg = 0; i_leg < 4; ++i_leg) {
      int idx_F = i_leg * nAineq_unit;
      int idx_constraint_u = i_leg * 3;
      Fi.block(idx_F, idx_constraint_u, nAineq_unit, 3) = Aineq_unit;
      hi.segment(idx_F, 6) << mu * Ut(2 + i_leg * 3) - Ut(i_leg * 3), mu * Ut(2 + i_leg * 3) + Ut(i_leg * 3),
          mu * Ut(2 + i_leg * 3) - Ut(1 + i_leg * 3), mu * Ut(2 + i_leg * 3) + Ut(1 + i_leg * 3),
          ub(i_leg, i_hor) - Ut(2 + i_leg * 3) + Ud(2 + i_leg * 3, i_hor),
          -lb(i_leg, i_hor) + Ut(2 + i_leg * 3) - Ud(2 + i_leg * 3, i_hor);
    }
    //        std::cout << "The matrix Fi is:\n" << Fi << "\n\n";
    int idx_A = (i_hor)*4 * nAineq_unit;
    int idx_z = (i_hor) * (nX + nU);
    Aineq.block(idx_A, idx_z, 4 * nAineq_unit, nU) = Fi;
    bineq.segment(idx_A, 4 * nAineq_unit) = hi;
  }
  //  std::cout << "The matrix H is:\n" << H << "\n\n";
  //  std::cout << "The matrix g is:\n" << g << "\n\n";
  //  std::cout << "The matrix Aeq is:\n" << Aeq << "\n\n";
  //  std::cout << "The matrix beq is:\n" << beq << "\n\n";
  //  std::cout << "The matrix Aineq is:\n" << Aineq << "\n\n";
  //  std::cout << "The matrix bineq is:\n" << bineq << "\n\n";
}
