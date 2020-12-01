#ifndef __QUATERNION_H
#define __QUATERNION_H

#include "setType.h"
#include "Matrix.h"

class quaternion {
  private:
    Complex qw, qx, qy, qz;
  public:
  	void q_sum(quaternion a, quaternion b);
  	void q_subtraction(quaternion a, quaternion b);
  	void q_multiplication(quaternion a, quaternion b);
  	void q_multE(quaternion a, Complex x);
  	void q_divE(quaternion a, Complex x);
  	void q_conjugate(quaternion q);
  	void q_identity(void);
  	Complex q_norm(void);
  	void q_inverse(quaternion q);
    void q_normalize(quaternion q);
    void q_dot_product(quaternion a, quaternion b);
    void q_2ea(Complex& roll, Complex& picth, Complex& yaw);
    void ea_2q(Complex roll, Complex picth, Complex yaw);
    Matrix q_2rot(void);
};

#endif
