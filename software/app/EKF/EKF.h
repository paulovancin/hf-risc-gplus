#include "setType.h"
#include "matrix.h"
#include "quaternions.h"


struct ekf_out EKF(typ_var X[7], typ_var acc[3], typ_var gyro[3], typ_var mag[3], struct Matrix Q, struct Matrix R, struct Matrix P_in, typ_var Ts);
