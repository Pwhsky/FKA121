#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#define PI 3.14159265359
typedef struct {
    double integral;
    double error;
} result_t;
/*Generate N points x_i uniformly on the interval [0,1] and evaluate the
integral with the Monte Carlo method. Perform the calculations for
N = 1e1 , 1e2 , 1e3 , and 1e4 , and estimate the error. Compare your
results with the exact value.*/
double function(double x){
    return x*(1-x);
}

result_t MC_with_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    double x;
    double y;
    double sum = 0;
    double squareSum= 0;  
    for (int i = 0; i < N; i++){
        x = gsl_rng_uniform(k);
        y = (1.0/PI)*acos(1.0-2.0*x);
        
        sum += (2.0/PI)*(y*(1-y))/(sin(PI*y));
        squareSum += pow((2.0/PI)*(y*(1-y))/(sin(PI*y)),2);
    }  
    double sigma = sqrt(squareSum/(double)N - pow(sum/(double)N,2))/(sqrt(N));
    double average = sum/(double)N;

    result.integral = average; // Write the integral here
    result.error = sigma;// Write the error here

    return result;
}
result_t MC_without_importance_sampling(int N, gsl_rng *k)
{
    result_t result;
    double sum = 0;
    double squareSum = 0;

    for (int i = 0; i < N; i++){
        double x = gsl_rng_uniform(k);
        sum += function(x);
        squareSum += pow((function(x)),2);
    }  
    double average = sum/(double)N;
    double sigma = sqrt(squareSum/(double)N - pow(sum/(double)N,2))/sqrt(N);

    double exact   = 1.0/6.0;
    result.integral = average; // Write the integral here
    result.error =  sigma;// Write the error here

    return result;
}

int main() {
    const gsl_rng_type *T;
    gsl_rng *k;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    k = gsl_rng_alloc(T);
    gsl_rng_set(k, time(NULL));

    int N[4] = {10,100,1000,10000};
    for(int i = 0; i < 4; i++ ){
        result_t with_sampling    = MC_with_importance_sampling(N[i],k); 
        result_t without_sampling = MC_without_importance_sampling(N[i],k); 
        printf("N = %d:\n",N[i]);
        printf("Without Importance sampling: Integral = %lf, Error = %lf \n",without_sampling.integral,without_sampling.error);    
        printf("With Importance sampling: Integral = %lf, Error = %lf \n \n",with_sampling.integral,with_sampling.error);

    }


    return 0;
}