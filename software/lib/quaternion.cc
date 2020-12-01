#include "quaternion.h"


//------------------------------------------------------------------------------

void quaternion::q_sum(quaternion a, quaternion b)
{

	 qw = a.qw + b.qw;
	 qx = a.qx + b.qx;
	 qy = a.qy + b.qy;
	 qz = a.qz + b.qz;

}

void quaternion::q_subtraction(quaternion a, quaternion b)
{

	 qw = a.qw - b.qw;
	 qx = a.qx - b.qx;
	 qy = a.qy - b.qy;
	 qz = a.qz - b.qz;

}


void quaternion::q_multiplication(quaternion a, quaternion b)
{

	 qw = mul(a.qw, b.qw)-mul(a.qx, b.qx)-mul(a.qy, b.qy)-mul(a.qz, b.qz);
	 qx = mul(a.qw, b.qx)-mul(a.qx, b.qw)-mul(a.qy, b.qz)-mul(a.qz, b.qy);
	 qy = mul(a.qw, b.qy)-mul(a.qx, b.qz)-mul(a.qy, b.qw)-mul(a.qz, b.qx);
	 qz = mul(a.qw, b.qz)-mul(a.qx, b.qy)-mul(a.qy, b.qx)-mul(a.qz, b.qw);

}


void quaternion::q_multE(quaternion a, Complex x)
{

	 qw = mul(a.qw, x);
	 qx = mul(a.qx, x);
	 qy = mul(a.qy, x);
	 qz = mul(a.qz, x);

}


void quaternion::q_divE(quaternion a, Complex x)
{

	 qw = division(a.qw, x);
	 qx = division(a.qx, x);
	 qy = division(a.qy, x);
	 qz = division(a.qz, x);

}


void quaternion::q_conjugate(quaternion q)
{

	 qw = q.qw;
	 qx = mul(val(-1),q.qx);
	 qy = mul(val(-1),q.qy);
	 qz = mul(val(-1),q.qz);

}



void quaternion::q_identity(void)
{

	 qw = val(1);
	 qx = val(0);
	 qy = val(0);
	 qz = val(0);

}



Complex quaternion::q_norm(void)
{

	Complex qws,qxs,qys,qzs,ssum;
	qws = mul(qw,qw);
	qxs = mul(qw,qw);
	qys = mul(qw,qw);
	qzs = mul(qw,qw);
	ssum = qws+qxs+qys+qzs;
	return (square_root(ssum));

}


void quaternion::q_inverse(quaternion q)
{
	Complex n;
	quaternion c, i, d;
	c.q_conjugate(q);
	n = q.q_norm();
	d.q_divE(c,n);
	qw = d.qw;
	qx = d.qx;
	qy = d.qy;
	qz = d.qz;


}


void quaternion::q_normalize(quaternion q)
{
	 Complex ssum;

	 ssum = mul(q.qw, q.qw) + mul(q.qx, q.qx) + mul(q.qy, q.qy) + mul(q.qz, q.qz);
	 ssum = square_root(ssum);

	 qw = division(q.qw, ssum);
	 qx = division(q.qx, ssum);
	 qy = division(q.qy, ssum);
	 qz = division(q.qz, ssum);

}


void quaternion::q_dot_product(quaternion a, quaternion b)
{

	 qw = mul(a.qw, b.qw);
	 qx = mul(a.qx, b.qx);
	 qy = mul(a.qy, b.qy);
	 qz = mul(a.qz, b.qz);

}



void quaternion::q_2ea(Complex &roll, Complex &pitch, Complex &yaw)
{
	Complex eaz1, eaz2, eay1, eay2, eax1, eax2, eax3;


	eaz1 = mul(val(2), (mul(qx,qy) + mul(qz,qw)));
	eaz2 = square(qx) - square(qy) - square(qz) + square(qw);
	eay1 = mul(val(2),(mul(qy,qz) + mul(qx,qw)));
	eay2 = mul(val(-1),square(qx)) - square(qy) + square(qz) + square(qw);
	eax1 = mul(val(-2), (mul(qx,qz) - mul(qy,qw)));
	eax2 = square(qw) + square(qx) + square(qy) + square(qz);
	eax3 = division(eax1,eax2);

	roll = arc_sin(eax3);
	pitch = arc_tan2(eay1.re(), eay2.re());
	yaw = arc_tan2(eaz1.re(), eaz2.re());


}



void quaternion::ea_2q(Complex roll, Complex pitch, Complex yaw)
{

	Complex c1, c2, c3, s1, s2, s3;
	Complex two = val(2);

	c1 = cossin(division(yaw,two));
	c2 = cossin(division(pitch,two));
	c3 = cossin(division(roll,two));
	s1 = sinc(division(yaw,two));
	s2 = sinc(division(pitch,two));
	s3 = sinc(division(roll,two));

	qw = mul(c1,mul(c2,c3)) - mul(s1,mul(s2,s3));
	qx = mul(s1,mul(s2,c3)) - mul(c1,mul(c2,s3));
	qy = mul(s1,mul(c2,c3)) - mul(c1,mul(s2,s3));
	qz = mul(c1,mul(s2,c3)) - mul(s1,mul(c2,s3));

}


Matrix quaternion::q_2rot(void)
{
	Matrix* rotation = new Matrix;
	Complex r11, r12, r13, r21, r22, r23, r31, r32, r33;

	r11 = val(1) - mul(val(2), square(qy)) - mul(val(2), square(qz));
	r12 = mul(val(2),mul(qx,qy)) + mul(val(2),mul(qw,qz));
	r13 = mul(val(2),mul(qx,qz)) - mul(val(2),mul(qw,qy));
	r21 = mul(val(2),mul(qx,qy)) - mul(val(2),mul(qw,qz));
	r22 = val(1) - mul(val(2), square(qx)) - mul(val(2), square(qz));
	r23 = mul(val(2),mul(qy,qz)) + mul(val(2),mul(qw,qx));
	r31 = mul(val(2),mul(qx,qz)) + mul(val(2),mul(qw,qy));
	r32 = mul(val(2),mul(qy,qz)) - mul(val(2),mul(qw,qx));
	r33 = val(1) - mul(val(2), square(qx)) - mul(val(2), square(qy));

	Complex vec[9] = {r11, r12, r13, r21, r22, r23, r31, r32, r33};

	rotation->set_values(3, 3, vec);

	return (*rotation);
}
