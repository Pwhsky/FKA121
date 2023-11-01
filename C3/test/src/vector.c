#include <stdio.h>
#include <math.h>
#include "vector.h"

double
norm(
	    double *v1,
	    unsigned int len
	   ){
		double sum = 0;
		for(int i = 0; i<len;i++){
			sum += pow(v1[i],2);
	   }
	   sum = sqrt(sum);
		return sum;
}

void
elementwise_addition(
		     double *res,
		     double *v1,
		     double *v2,
		     unsigned int len
	            )
{
    for(int i = 0; i < len; ++i){
	res[i] = v1[i] + v2[i];
    }
}

void
constant_multiplication(
			   double *res,
			   double *v1,
			   double a,
			   unsigned int len
	                  )
{
    for(int i = 0; i < len; ++i){
	res[i] = v1[i] * a;
    }
}

double
dot_product(
	    double *v1,
	    double *v2,
	    unsigned int len
	   )
{
    double result = 0;
    for(int i = 0; i < len; ++i){
	result += v1[i] * v2[i];
    }
    return result;
}
