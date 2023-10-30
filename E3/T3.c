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
double function(double x,double y, double z){
    return pow(x,2) + pow(x,2)*pow(y,2) + pow(x,2)*pow(y,2)*pow(z,2);
}

double weight_function(double x,double y, double z){
    double exponent = -3.0/2.0;
    return exp(-pow(x,2) - pow(y,2) - pow(z,2))/pow(PI, (3.0/2.0));
}



void write_to_file(char *fname,double* x, double* y ) {
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "x,y\n");
    for (int i = 0; i < nPoints; ++i) {
        fprintf(fp, "%f,%f\n",x[i] ,y[i] );
    }
    fclose(fp);
}
int main() {

    //set up random number generator
    const gsl_rng_type *T;
    gsl_rng *B;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    B = gsl_rng_alloc(T);
    gsl_rng_set(B, time(NULL));

    double* x      = (double*)malloc(sizeof(double)*nPoints);
    double* y      = (double*)malloc(sizeof(double)*nPoints);
    double* z      = (double*)malloc(sizeof(double)*nPoints);
    double* trials = (double*)malloc(sizeof(double)*3);
    double* rands  = (double*)malloc(sizeof(double)*3);
    x[0] = gsl_rng_uniform(B);
    y[0] = gsl_rng_uniform(B);
    z[0] = gsl_rng_uniform(B);
    double delta = 2.1;
    double sum = 0.0;
    size_t counter =0;
    for (int i = 0; i < nPoints-1; i++){
        //Sample rands and propose to trial
        for(int j = 0; j<3; j++){
            rands[j] = gsl_rng_uniform(B);
        }
        trials[0] = x[i] + delta*(rands[0]-0.5);
        trials[1] = y[i] + delta*(rands[1]-0.5);
        trials[2] = z[i] + delta*(rands[2]-0.5);
        //get corresponding weights for trials and compare to the current points.
        double w       = weight_function(x[i],y[i],z[i]);
        double w_trial = weight_function(trials[0],trials[1],trials[2]);
        double decision = gsl_rng_uniform(B);
        if (w_trial > w || (w_trial/w) > decision){
            x[i+1] = trials[0];
            y[i+1] = trials[1];
            z[i+1] = trials[2];
            counter++;
        }else {
            x[i+1] =  x[i];
            y[i+1] =  y[i];
            z[i+1] =  z[i];
        }
        //Profit
    }  
    //Convergence does not happen instantly, so discard the first 2000 points in a "burn"
    int nBurnoff = 2000;
    for (int i = nBurnoff; i< nPoints; i++){
        sum += function(x[i],y[i],z[i]);
    }
    double average = sum/(double)(nPoints-nBurnoff);
    double accept_rate = (double)counter/nPoints;
    double exact   = 7.0/8.0;
    printf("Exact value = %lf \n",exact);
    printf("Approximated value = %lf \n",average);
    printf("Acceptance rate = %lf \n",accept_rate);

    return 0;
}