#ifndef _COMPLEX_H_
#define _COMPLEX_H_

#include "setType.h"

class Complex {
private:
	typ_var real, img;
public:
	Complex(void);
	Complex(typ_var r);
	Complex(typ_var r, typ_var i);
	typ_var re(void);
	typ_var im(void);
	Complex conj(void);
	typ_var abs(void);
	typ_var arg(void);
	typ_var norm(void);
	Complex operator+(Complex z);
	Complex operator+(typ_var a);
	friend Complex operator+(typ_var a, Complex z);
	Complex operator-(Complex z);
	Complex operator-(typ_var a);
	friend Complex operator-(typ_var a, Complex z);
	Complex operator*(Complex z);
	Complex operator*(typ_var a);
	friend Complex operator*(typ_var a, Complex z);
	Complex operator/(Complex z);
	Complex operator/(typ_var a);
	friend Complex operator/(typ_var a, Complex z);
	const Complex& operator+=(const Complex& z);
  const Complex& operator-=(const Complex& z);
  int operator==(Complex z);
  Complex sqrt_C(void);
	Complex log_C(void);
	Complex exp_C(void);
	Complex pow_C(typ_var c);
	Complex sin_C(void);
	Complex cos_C(void);
	Complex tan_C(void);
	Complex sec_C(void);
	Complex csc_C(void);
	Complex cot_C(void);
	Complex sinh_C(void);
	Complex cosh_C(void);
	Complex tanh_C(void);
	Complex sech_C(void);
	Complex csch_C(void);
	Complex coth_C(void);
	Complex asin_C(void);
	Complex acos_C(void);
	Complex atan_C(void);
	Complex asinh_C(void);
	Complex acosh_C(void);
	Complex atanh_C(void);
	void print_C(void);

};


#endif
