#include "complex.h"
#include "fixed.h"

static Complex j = Complex(val(0),val(1));
static Complex i = Complex(val(0),val(1));


Complex::Complex(void) {real = val(0); img = val(0);}
Complex::Complex(typ_var r) {real = r; img = val(0);}
Complex::Complex(typ_var r,typ_var i) {real = r; img = i;}

typ_var Complex::re(void)
{
	return(this->real);
}

typ_var re(Complex z)
{
	return z.re();
}

typ_var Complex::im(void)
{
	return(this->img);
}

typ_var im(Complex z)
{
	return z.im();
}

Complex Complex::conj(void)
{
	return Complex(this->real,-this->img);
}

Complex conj(Complex z)
{
	return Complex(z.re(),-z.im());
}

typ_var Complex::abs(void)
{
	return fix_sqrt(fix_mul(real,real)+fix_mul(img,img));
}

typ_var abs(Complex z)
{
	return fix_sqrt(fix_mul(z.re(),z.re())+fix_mul(z.im(),z.im()));
}

typ_var Complex::arg(void)
{
	return fix_atan2(this->img,this->real);
}

typ_var arg(Complex z)
{
	return fix_atan2(z.im(),z.re());
}

typ_var Complex::norm(void)
{
	return fix_mul(real,real)+fix_mul(img,img);
}

typ_var norm(Complex z)
{
	return fix_mul(z.re(),z.re())+fix_mul(z.im(),z.im());
}


Complex Complex::operator+(Complex z)
{
  return Complex(this->real+z.real,this->img+z.img);
}


Complex Complex::operator+(typ_var a)
{
  return Complex(this->real+a,this->img);
}


Complex operator+(typ_var a, Complex z)
{
  return Complex(a+z.re(),z.im());
}


Complex Complex::operator-(Complex z)
{
  return Complex(this->real-z.real,this->img-z.img);
}


Complex Complex::operator-(typ_var a)
{
  return Complex(this->real-a,this->img);
}


Complex operator-(typ_var a, Complex z)
{
  return Complex(a-z.re(),-z.im());
}


Complex Complex::operator*(Complex z)
{
  return Complex((fix_mul(this->real,z.real) - fix_mul(this->img,z.img)),
		 (fix_mul(this->real,z.img)+fix_mul(this->img,z.real)));
}


Complex Complex::operator*(typ_var a)
{
  return Complex(fix_mul(real,a),fix_mul(img,a));
}

Complex operator*(typ_var a, Complex z)
{
  return Complex(fix_mul(a,z.re()),fix_mul(a,z.im()));
}


Complex Complex::operator/(Complex z)
{
  Complex top((*this)*z.conj());
  typ_var bottom(z.norm());
  Complex quo(top/bottom);
  return quo;
}

Complex Complex::operator/(typ_var a)
{
  return Complex(fix_div(this->real,a),fix_div(this->img,a));
}

Complex operator/(typ_var a, Complex z)
{
  Complex top((a)*z.conj());
  typ_var bottom(z.norm());
  Complex quo(top/bottom);
  return quo;
}


const Complex& Complex::operator+=(const Complex& z)
{
  this->real+=z.real;
  this->img+=z.img;
  return *this;
}


const Complex& Complex::operator-=(const Complex& z)
{
  this->real-=z.real;
  this->img-=z.img;
  return *this;
}


int Complex::operator==(Complex z)
{
  if (this->real == z.re() && this->img == z.im()) {
    return 1;
  }
  else {
    return 0;
  }
}

Complex Complex::sqrt_C(void)
{
  typ_var zsqre,zsqim;

  zsqre = fix_sqrt(fix_mul(val(0.5),(this->abs()+this->re())));
  zsqim = fix_sqrt(fix_mul(val(0.5),(this->abs()-this->re())));

  if (this->im() > val(0.000015)) {
    return Complex(zsqre,zsqim);
  }
  else {
    return Complex(zsqre,val(0));
  }
}


Complex Complex::log_C(void)
{
  if (this->re() < 0 && this->im() == 0.0) {
    return fix_ln(this->abs())+j*PI;
  }
  else {
    return fix_ln(this->abs())+j*this->arg();
  }
}


Complex Complex::exp_C(void)
{
  return fix_exp(this->re())*(fix_cos(this->im())+j*fix_sin(this->im()));
}


Complex Complex::pow_C(typ_var c)
{
  Complex x, y;
  x = this->log_C();
  y = c*x;
  return y.exp_C();
}


Complex Complex::sin_C(void)
{
  Complex x1, x2, x3, x4, y1, y2, y3, y4, z;
  z.real = this->re();
  z.img = this->im();
  x1 = j*z;
  x2 = x1.exp_C();
  x3 = (-1)*j;
  x4 = 0.5*x2*x3;
  y1 = (-1)*j*z;
  y2 = y1.exp_C();
  y3 = (-1)*j;
  y4 = 0.5*y2*y3;

  return (x4-y4);
}

Complex Complex::cos_C(void)
{
  Complex x1, x2, x3, y1, y2, y3, z;
  z.real = this->re();
  z.img = this->im();
  x1 = j*z;
  x2 = x1.exp_C();
  x3 = 0.5*x2;
  y1 = (-1)*j*z;
  y2 = y1.exp_C();
  y3 = 0.5*y2;

  return (x3+y3);
}

Complex Complex::tan_C(void)
{
	return this->sin_C()/this->cos_C();
}

Complex Complex::sec_C(void)
{
	return 1/this->cos_C();
}

Complex Complex::csc_C(void)
{
	return 1/this->sin_C();
}

Complex Complex::cot_C(void)
{
	return this->cos_C()/this->sin_C();
}

Complex Complex::sinh_C(void)
{
	Complex x, z;
	z.real = this->re();
    z.img = this->im();
    x = (-1)*z;
	return ((z.exp_C()-x.exp_C())/2);
}

Complex Complex::cosh_C(void)
{
	Complex x, z;
	z.real = this->re();
    z.img = this->im();
    x = (-1)*z;
	return ((z.exp_C()+x.exp_C())/2);
}

Complex Complex::tanh_C(void)
{
	return this->sinh_C()/this->cosh_C();
}

Complex Complex::sech_C(void)
{
	return 1/this->cosh_C();
}

Complex Complex::csch_C(void)
{
	return 1/this->sinh_C();
}

Complex Complex::coth_C(void)
{
	return this->cosh_C()/this->sinh_C();
}

Complex Complex::asin_C(void)
{
	Complex x1, x2, x3, x4, x5, z;
	z.real = this->re();
    z.img = this->im();
    x1 = (-1)*j;
    x2 = 1-z*z;
    x3 = x2.sqrt_C();
    x4 = j*z+x3;
    x5 = x4.log_C();
	return (x1*x5);
}

Complex Complex::acos_C(void)
{
	Complex x1, x2, x3, x4, x5, z;
	z.real = this->re();
    z.img = this->im();
    x1 = (1)*j;
    x2 = z*z-1;
    x3 = x2.sqrt_C();
    x4 = z+x3;
    x5 = x4.log_C();
	return (x1*x5);
}

Complex Complex::atan_C(void)
{
	Complex x1, x2, x3, x4, x5, z;
	z.real = this->re();
    z.img = this->im();
    x1 = j+z;
    x2 = j-z;
    x3 = x1/x2;
    x4 = x3.log_C();
    x5 = (0.5*j)*x4;
	return (x5);
}

Complex Complex::asinh_C(void)
{
	Complex x1, x2, x3, x4, z;
	z.real = this->re();
    z.img = this->im();
    x1 = z*z+1;
    x2 = x1.sqrt_C();
    x3 = z+x2;
    x4 = x3.log_C();
	return (x4);
}

Complex Complex::acosh_C(void)
{
	Complex x1, x2, x3, x4, z;
	z.real = this->re();
    z.img = this->im();
    x1 = z*z-1;
    x2 = x1.sqrt_C();
    x3 = z+x2;
    x4 = x3.log_C();
	return (x4);
}

Complex Complex::atanh_C(void)
{
	Complex x1, x2, x3, x4, x5, z;
	z.real = this->re();
    z.img = this->im();
    x1 = 1+z;
    x2 = 1-z;
    x3 = x1/x2;
    x4 = x3.log_C();
    x5 = 0.5*x4;
	return (x5);
}

void Complex::print_C(void)
{
	if (this->im() > 0){
		fix_print(this->re());
		printf("+");
		fix_print(this->im());
		printf("i");
	}

	else if (this->im() == 0){
		fix_print(this->re());
	}

	else {
		fix_print(this->re());
		fix_print(this->im());
		printf("i");
	}

}
