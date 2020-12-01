#include "setType.h"

void qnormalize(typ_var *quat);
void qtimes(typ_var *quat, typ_var a[4], typ_var b[4]);
void imu2q(typ_var *quat, typ_var acc[3], typ_var mag[3]);
typ_var q2ea(typ_var quat[4]);
