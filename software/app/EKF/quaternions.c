#include "quaternions.h"

//------------------------------------------------------------------------------

void qnormalize(typ_var *quat)
{

	 typ_var ssum;

	 ssum = mul(quat[0], quat[0]) + mul(quat[1], quat[1]) + mul(quat[2], quat[2]) + mul(quat[3], quat[3]);
	 ssum = square_root(ssum);


	 quat[0] = div(quat[0], ssum);
	 quat[1] = div(quat[1], ssum);
	 quat[2] = div(quat[2], ssum);
	 quat[3] = div(quat[3], ssum);


}

void qtimes(typ_var *quat, typ_var a[4], typ_var b[4])
{

		quat[0] = mul(a[0],b[0]) - mul(a[1],b[1]) - mul(a[2],b[2]) - mul(a[3],b[3]);
		quat[1] = mul(a[1],b[0]) + mul(a[0],b[1]) - mul(a[3],b[2]) + mul(a[2],b[3]);
		quat[2] = mul(a[2],b[0]) + mul(a[3],b[1]) + mul(a[0],b[2]) - mul(a[1],b[3]);
		quat[3] = mul(a[3],b[0]) - mul(a[2],b[1]) + mul(a[1],b[2]) + mul(a[0],b[3]);

}

void imu2q(typ_var *quat, typ_var acc[3], typ_var mag[3])
{

		typ_var tau, lx, ly, lz, ax, ay, az;
		typ_var qacc[4], qmag[4], qacc_aux[4], qmag_aux[4], q[4];

		lx = mag[0];
		ly = mag[1];
		lz = mag[2];

		ax = acc[0];
		ay = acc[1];
		az = acc[2];


		tau = (mul(lx,lx) + mul(ly,ly));

		if (az >= 0){
			qacc_aux[0] = square_root((div((az+val(1)),val(2))));
			qacc_aux[1] = mul(val(-1),(div(ay,(mul(val(2),(square_root((az+val(1)))))))));
			qacc_aux[2] = div(ax, (mul(val(2),(square_root((az+val(1)))))));
			qacc_aux[3] = val(0);
			qnormalize(qacc_aux);
		}

		if (az < 0){
			qacc_aux[0] = -(div(ay,(mul(val(2),square_root((az-val(1)))))));
			qacc_aux[1] = square_root(div((az-val(1)), val(2)));
			qacc_aux[2] = val(0);
			qacc_aux[3] = div(ax, (mul(val(2),(square_root((az-val(1)))))));
			qnormalize(qacc_aux);
		}

		if (lx >= 0){
			qmag_aux[0] = div((square_root(tau + mul(lx, square_root(tau)))), (square_root(mul(val(2),tau))));
			qmag_aux[1] = val(0);
			qmag_aux[2] = val(0);
			qmag_aux[3] = div(ly, (mul(square_root(val(2)), square_root(tau + mul(lx, square_root(tau))))));
			qnormalize(qmag_aux);
		}

		if (lx < 0){
			qmag_aux[0] = div(ly, (mul(square_root(val(2)),(square_root(tau - mul(lx,(square_root(tau))))))));
			qmag_aux[1] = val(0);
			qmag_aux[2] = val(0);
			qmag_aux[3] = div((square_root(tau - mul((lx),(square_root(tau))))),(square_root(mul(val(2),tau))));
			qnormalize(qmag_aux);
		}


		qtimes(q, qacc_aux, qmag_aux);


		quat[0] = q[0];
		quat[1] = q[1];
		quat[2] = q[2];
		quat[3] = q[3];

}


typ_var q2ea(typ_var quat[4])
{

typ_var out[3];

  // Roll (x-axis rotation)
  typ_var sinr_cosp = mul(val(2), (mul(quat[0], quat[1]) + mul(quat[2], quat[3])));
  typ_var cosr_cosp = val(1) + (mul(val(2), (mul(quat[1], quat[1]) + mul(quat[2], quat[2]))));
  out[0] = arc_tan2(sinr_cosp, cosr_cosp);
  out[0] = mul(out[0], div(val(180), PI));

  // Pitch (y-axis rotation)
  typ_var sinp = mul(val(2), (mul(quat[0], quat[2]) - mul(quat[3], quat[1])));
  if (fabs(sinp) >= val(1)){
    out[1] = cpsign(HALF_PI, sinp); // If GIMBAL LOCK
    out[1] = mul(out[1], div(val(180), PI));
  }
  else{
    out[1] = arc_sin(sinp);
    out[1] = mul(out[1], div(val(180), PI));
  }

  // Yaw (z-axis rotation)
  typ_var siny_cosp = mul(val(2), (mul(quat[0], quat[3]) + mul(quat[1], quat[2])));
  typ_var cosy_cosp = val(1) - (mul(val(2), (mul(quat[2], quat[2]) + mul(quat[3], quat[3]))));
  out[2] = arc_tan2(siny_cosp, cosy_cosp);
  out[2] = mul(out[2], div(val(180), PI));


return out;

}
