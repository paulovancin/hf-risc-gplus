#include <hf-risc.h>
#include"fixedptc.h"

void main(void){

float i;
fixedpt test_a, test_b, result;
char str, str2;

test_a = -2*256;
test_b = -1*256;

result = test_a + test_b;


i = fixedpt_tofloat(result);
fixedpt_str(result, &str, 10);
printf("Fixed Point:%s\n", &str);
fixedpt_str(FIXEDPT_ONE, &str2, 10);
printf("Fixed Point One:%s\n", &str2);
printf("float:%f\n", i);

}
