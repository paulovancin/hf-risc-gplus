#include "EKF.h"

struct ekf_out EKF(typ_var X[7], typ_var acc[3], typ_var gyro[3], typ_var mag[3], struct Matrix Q, struct Matrix R, struct Matrix P_in, typ_var Ts)
{

  typ_var quat_measured[4], gx, gy, gz, wx, wy, wz, w[3], q[4], q2[4], n, e[3], x_norm[7];

  struct Matrix O34, O43, O33, I, W, quats, eps, Se, Sw, omega, E, F, A, O;
  struct Matrix Neg_transp_W, Neg_transp_E, Aux_E, Aux_F, XMat, FTs, A_Trans, H;
  struct Matrix H_Trans, K1, K2, K, eye4, Aux_X, Aux_P, Aux_P1, Aux_P2, eye7, Aux_X2, es;
  struct Matrix P, P1, P2, XMat2, K2_inv, Pout;
  struct ekf_out OUT;

  // Transform Sensor Data Into Quaternions

  imu2q(quat_measured, acc, mag);

  // Set Necessary Matrices

  O34 = zeros(3, 4);
  O43 = zeros(4, 3);
  O33 = zeros(3, 3);
  O = zeros(1, 1);
  I = setEye(3);

  w[0] = gyro[0] - X[4];
  w[1] = gyro[1] - X[5];
  w[2] = gyro[2] - X[6];

  W = set_values(3, 1, w);

  q[0] = X[0];
  q[1] = X[1];
  q[2] = X[2];
  q[3] = X[3];

  n = q[0];
  e[0] = q[1];
  e[1] = q[2];
  e[2] = q[3];

  quats = set_values(4, 1, q);
  es = set_values(3, 1, e);

  Se = crossProduct3(e[0], e[1], e[2]);
  Sw = crossProduct3(w[0], w[1], w[2]);
  Neg_transp_W = transposed(W);
  Neg_transp_W = multE(Neg_transp_W, val(-1));
  Sw = multE(Sw, val(-1));

  int conf_omega[4] = {2,2,4,4};
  omega = customMat(4, conf_omega, O, Neg_transp_W, W, Sw);
  omega = multE(omega, val(0.5));

  Neg_transp_E = transposed(es);
  Neg_transp_E = multE(Neg_transp_E, val(-1));
  Aux_E = multE(I, n);
  Aux_E = sum(Aux_E, Se);
  int conf_E[4] = {2,1,4,3};
  E = customMat(2, conf_E, Neg_transp_E, Aux_E);
  E = multE(E, val(0.5));

  Aux_F = multiplication(omega, quats);

  int conf_F[4] = {4,1,7,1};
  F = customMat(4, conf_F, Aux_F, O, O, O);

  int conf_A[4] = {2,2,7,7};
  A = customMat(4, conf_A, omega, E, O34, O33);


  // // ------------------ Prediction --------------------------------

  XMat = set_values(7, 1, X);
  FTs = multE(F, Ts);
  XMat = sum(XMat, FTs);
  A_Trans = transposed(A);
  P1 = multiplication(A, P_in);
  P2 = multiplication(P1, A_Trans);
  P = sum(P2, Q);


  // ------------------ Kalman Gain --------------------------------

  eye4 = setEye(4);
  eye7 = setEye(7);

  int conf_H[4] = {1,2,4,7};
  H = customMat(2, conf_H, eye4, O43);
  H_Trans = transposed(H);
  K1 = multiplication(P, H_Trans);
  K2 = multiplication(H, K1);
  K2 = sum(K2, R);
  K2_inv = Invert(K2);

  K = multiplication(K1, K2_inv);

  // ------------------ Correction --------------------------------

  for(int i = 0; i < 4; i++){
    q[i] = get_value(XMat, i, 0);
  }

  q2[0] = quat_measured[0] - q[0];
  q2[1] = quat_measured[1] - q[1];
  q2[2] = quat_measured[2] - q[2];
  q2[3] = quat_measured[3] - q[3];

  quats = set_values(4, 1, q2);
  Aux_X2 = multiplication(K, quats);
  XMat2 = sum(XMat, Aux_X2);

  Aux_P1 = multiplication(K, H);
  Aux_P2 = subtraction(eye7, Aux_P1);
  Pout = multiplication(Aux_P2, P);

  for(int i = 0; i < 4; i++){
    q[i] = get_value(XMat2, i, 0);
  }

  for(int i = 0; i < 3; i++){
    x_norm[i] = get_value(XMat2, (i+4), 0);
  }

  qnormalize(q);

  typ_var Xo[7] = {q[0], q[1], q[2], q[3], x_norm[0], x_norm[1], x_norm[2]};
  XMat2 = set_values(7, 1, Xo);

  OUT.states = copy(XMat2);
  OUT.MatP = copy(Pout);

  SHOW(q[0]);
  printf("\n");

  return OUT;

}

// ---------------------------------------------------------------------------
