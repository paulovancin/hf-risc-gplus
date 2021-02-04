#ifndef __FIXED_H
#define __FIXED_H


#ifdef __cplusplus
extern "C" {
#endif

#include "hf-risc.h"


 typedef unsigned char			  	uint8_t;
 typedef char				        int8_t;
 typedef unsigned short int			uint16_t;
 typedef short int			      	int16_t;
 typedef unsigned int			  	uint32_t;
 typedef int				        int32_t;
 typedef unsigned long long			uint64_t;
 typedef long long			      	int64_t;
 typedef unsigned long			  	size_t;
 typedef void				        (*funcptr)();
 typedef int32_t		          	fixed_t;


#ifndef FIX_IBITS
#define FIX_IBITS	16
#endif

#ifndef FIX_MULDIV_WIDTH
#define FIX_MULDIV_WIDTH	64
#endif

#if FIX_MULDIV_WIDTH != 32 && FIX_MULDIV_WIDTH != 64
#error "FIX_MULDIV_WIDTH should be either 32 or 64"
#endif

#if FIX_IBITS < 16 && FIX_MULDIV_WIDTH == 32
#error "FIX_IBITS must be greater or equal to 16 when FIX_MULDIV_WIDTH is 32"
#endif

#define FIX_FBITS		(32 - FIX_IBITS)
#define FIX_FMASK		(((fixed_t)1 << FIX_FBITS) - 1)
#define FIX_ONE			((fixed_t)((fixed_t)1 << FIX_FBITS))

#define fix_val(V)		((fixed_t)((V) * FIX_ONE + ((V) >= 0 ? 0.5 : -0.5)))
#define fix_int(F)		((F) >> FIX_FBITS)
#define fix_frac(A)		((fixed_t)(A) & FIX_FMASK)
#define fix_abs(A)		((A) < 0 ? -(A) : (A))
#define fix_add(A,B)	((A) + (B))
#define fix_sub(A,B)	((A) - (B))

#define FIX_HALF		(FIX_ONE >> 1)
#define FIX_TWO			(FIX_ONE + FIX_ONE)
#define FIX_PI			fix_val(3.14159265358979323846)
#define FIX_TWO_PI		fix_val(6.28318530717958647692)
#define FIX_HALF_PI		fix_val(1.57079632679489661923)
#define FIX_FOURTH_PI	fix_val(0.78539816339744830961)
#define FIX_E			fix_val(2.71828182845904523536)
#define FIX_LN2			fix_val(0.69314718055994530942)
#define FIX_LN2_INV		fix_val(1.44269504088896340736)

#if FIX_MULDIV_WIDTH == 64
typedef	int64_t	fixedd_t;

#define fix_mul(A,B)		((fixed_t)(((fixedd_t)(A) * (fixedd_t)(B)) >> FIX_FBITS))
#define fix_div(A,B)		((fixed_t)(((fixedd_t)(A) << FIX_FBITS) / (fixedd_t)(B)))

#else
fixed_t fix_mul(fixed_t x, fixed_t y);
fixed_t fix_div(fixed_t x, fixed_t y);
#endif

void fixtoa(fixed_t a, char *str, int32_t dec);
fixed_t fix_sqrt(fixed_t a);
fixed_t fix_exp(fixed_t fp);
fixed_t fix_ln(fixed_t fp);
fixed_t fix_log(fixed_t fp, fixed_t base);
fixed_t fix_pow(fixed_t fp, fixed_t exp);
fixed_t fix_rad(fixed_t deg);
fixed_t fix_sin(fixed_t rad);
fixed_t fix_cos(fixed_t rad);
fixed_t fix_tan(fixed_t rad);
static fixed_t xatan(fixed_t arg);
static fixed_t satan(fixed_t arg);
fixed_t fix_atan(fixed_t arg);
fixed_t fix_atan2(fixed_t arg1, fixed_t arg2);
fixed_t fix_asin(fixed_t arg);
fixed_t fix_acos(fixed_t arg);
fixed_t fix_sinh(fixed_t arg);
fixed_t fix_cosh(fixed_t arg);
fixed_t fix_tanh(fixed_t arg);
void fix_print(fixed_t input);
fixed_t fix_copysign(fixed_t arg1, fixed_t arg2);
fixed_t fix_sign(fixed_t arg1);
fixed_t fix_norm(int size, fixed_t vec[]);
fixed_t fix_square(fixed_t arg1);



#ifdef __cplusplus
}

#endif


#endif
