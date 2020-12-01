#ifndef __SETTYPE_H
#define __SETTYPE_H

#include "hf-risc.h"
#include "new.h"
#define NEW_SIZE 5000

#ifdef TFIXED //===> FIXED POINT !!!!

#include "fixed.h"
typedef fixed_t typ_var;
#include "complex.h"

#define show(A)         A.print_C()
#define mul(A,B)		    A*B
#define division(A,B) 	A/B
#define val(A)          fix_val(A)
#define square_root(A)  A.sqrt_C()
#define power(A,B)      A.pow_C(B)
#define fabs(A)         A.abs()
#define sinc(A)      	  A.sin_C()
#define cossin(A)      	A.cos_C()
#define arc_sin(A)      A.asin_C()
#define arc_tan2(A,B)   fix_atan2(A,B)
#define HALF_PI         fix_val(1.57079632679489661923)
#define PI              fix_val(3.14159265358979323846)
#define square(A)       A*A
#define swap(A,B)       {y=(A);(A)=(B);(B)=y;}


#endif


#endif
