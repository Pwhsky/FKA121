#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>

/*Generate N points x_i uniformly on the interval [0,1] and evaluate the
integral with the Monte Carlo method. Perform the calculations for
N = 1e1 , 1e2 , 1e3 , and 1e4 , and estimate the error. Compare your
results with the exact value.*/
double function(double x){
    return x*(1-x);
}

int main() {
    srand(time(NULL));
    size_t nPoints = 10000;
    double stepSize = 1/nPoints;
    double* linspace = (double*)malloc(sizeof(double)*nPoints);
    double sum = 0;
    for (int i = 0; i < nPoints; i++){
        sum += function( (double)rand()/RAND_MAX );
    }  
    double average = sum/nPoints;
    double exact   = 1.0/6.0;
    double error   = fabs((exact - average)/exact)*100;
    printf("Exact value = %lf \n",exact);
    printf("Approximated value = %lf \n",average);
    printf("relative error = %lf%\n",error);

    return 0;
}