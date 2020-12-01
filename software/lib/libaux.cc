#include "libaux.h"

Complex cpsign(Complex arg1, Complex arg2){
	Complex out;
	if (arg2.re() >= 0){
		if (arg1.re() >= 0){
			out = arg1;
			return out;
		}
		else{
			out = (-1.0)*arg1;
			return out;
		}
	}

	else{
		if (arg1.re() >= 0){
			out = (-1.0)*arg1;
			return out;
		}
		else{
			out = arg1;
			return out;
		}
	}
}


Complex sign(Complex arg1)
{
	Complex out;
	if (arg1.re() > 0){
		out = val(1);
	}
	else if (arg1.re() == 0){
		out = val(0);
	}
	else{
		out = val(-1);
	}

	return out;

}

Complex norm_vector(int size, Complex vec[]){
	Complex aux, out;

	for(int i = 0; i < size; i++) {
		aux += square(vec[i]);
	}

	out = square_root(aux);
	return out;
}
