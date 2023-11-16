#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#define PI 3.14159265359
#define nPoints 1000000
/*Generate N points x_i uniformly on the interval [0,1] and evaluate the
integral with the Monte Carlo method. Perform the calculations for
N = 1e1 , 1e2 , 1e3 , and 1e4 , and estimate the error. Compare your
results with the exact value.*/
typedef struct{
    double  weight;
    double  function_value;
    int accepted;
} result_t;

double weight(double* x){
    return exp(-x[0]*x[0] - x[1]*x[1] - x[2]*x[2]);
}

double function(
                double *x)
{

 return x[0]*x[0] + x[0]*x[1]*x[0]*x[1] + x[0]*x[1]*x[2]*x[0]*x[1]*x[2];
}

result_t MCMC_step_displace_all(
                                 double *x,
                                 double delta,
                                 gsl_rng *k)
{

    result_t result; 
    double* rands = (double*)malloc(sizeof(double)*3);
    double* trials = (double*)malloc(sizeof(double)*3);

    rands[0] = gsl_rng_uniform(k); 
    rands[1] = gsl_rng_uniform(k); 
    rands[2] = gsl_rng_uniform(k);

    for (int i = 0; i <3; i++){
        trials[i] = x[i] + delta*(rands[i]-0.5);
    }

    double w       = weight(x);
    double w_trial = weight(trials);
    double decision = gsl_rng_uniform(k);

    if (w_trial>w ||(w_trial/w) > decision){

        for(int i = 0; i<3; i++){
            x[i] = trials[i];
        }
        result.accepted = 1;
    }else {
        result.accepted = 0;
    }
    result.weight         = weight(x);
    result.function_value = function(x);
    return result;
}


int main() {

    //set up random number generator
    const gsl_rng_type *T;
    gsl_rng *k;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    k = gsl_rng_alloc(T);
    gsl_rng_set(k, time(NULL));

    double* x     = (double*)malloc(sizeof(double)*3);
    int N = 1000000;
    x[0] = gsl_rng_uniform(k);
    x[1] = gsl_rng_uniform(k);
    x[2] = gsl_rng_uniform(k);
    double delta = 2.1;
    double sum = 0;

    //burn off starting points
    for (int t = 0; t < 20000; t++){
        result_t warmup = MCMC_step_displace_all(x,delta,k);
    }
    //Compute integral
    for (int t = 0; t < N; t++){
        result_t results = MCMC_step_displace_all(x,delta,k);
        sum += results.function_value;
    }
    double integral = sum/(double)N;


    //Convergence does not happen instantly, so discard the first 2000 points in a "burn"

    double exact   = 7.0/8.0;
    printf("Exact value = %lf \n",exact);
    printf("Approximated value = %lf \n",integral);


    return 0;
}