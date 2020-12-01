#include "EKF.h"
#include "out.h"


void main(int argc, char* argv[])
{

  int samples;
  samples = 431;
  typ_var Ts, q_init[4];
  typ_var ACC[3], GYR[3], MAG[3];
  typ_var ACCINIT[3], GYRINIT[3], MAGINIT[3];
  Ts = val(0.05);

  struct Matrix Q, R, Pin;
  struct ekf_out ST;


  // -------------------Setting Covariance Matrices-----------------------------


  R = setEye(4);
  R = multE(R, val(2));
  Q = setEye(7);
  Q = multE(Q, val(0.1));


  // -------------------Setting Initial States----------------------------------

  Pin = setEye(7);

  // Get Sensor Data For Initial Values - For Now Reading Data From "out.h"

  ACCINIT[1] = mul(val(-1),accx[0]);
  ACCINIT[0] = mul(val(-1),accy[0]);
  ACCINIT[2] = accz[0];
  GYRINIT[1] = mul(val(-1),gyrx[0]);
  GYRINIT[0] = mul(val(-1),gyry[0]);
  GYRINIT[2] = gyrz[0];
  MAGINIT[1] = mul(val(-1),magx[0]);
  MAGINIT[0] = mul(val(-1),magy[0]);
  MAGINIT[2] = magz[0];

  imu2q(q_init, ACCINIT, MAGINIT);
  typ_var Xin[7] = {q_init[0], q_init[1], q_init[2], q_init[3], val(0), val(0), val(0)};
  ST = EKF(Xin, ACCINIT, GYRINIT, MAGINIT, Q, R, Pin, Ts);
  //print_matrix(ST.states);


  // -----------------------Running The Filter----------------------------------


  for(int i = 1; i < 431; i++){

      // Get Sensor Data - For Now Reading Data From "out.h"

      ACC[1] = mul(val(-1),accx[i]);
      ACC[0] = mul(val(-1),accy[i]);
      ACC[2] = accz[i];
      GYR[1] = mul(val(-1),gyrx[i]);
      GYR[0] = mul(val(-1),gyry[i]);
      GYR[2] = gyrz[i];
      MAG[1] = mul(val(-1),magx[i]);
      MAG[0] = mul(val(-1),magy[i]);
      MAG[2] = magz[i];

      // Get Previous States and P Matrix

      for(int j = 0; j < 7; j++){
        Xin[j] = get_value(ST.states, j, 0);
      }

      Pin = copy(ST.MatP);

      // Call EKF

      ST = EKF(Xin, ACC, GYR, MAG, Q, R, Pin, Ts);

      // Printing The Quaternions and Gyro Bias - Later Put Quaternions Into A Vector For PID
      // Also Need to Send Timestamp for PID
      // Maybe Pass Barometer Data Together For Sync

      //print_matrix(ST.states);

    }

}
