#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#define PI 3.14159265359
#define nPoints 10000
/*Generate N points x_i uniformly on the interval [0,1] and evaluate the
integral with the Monte Carlo method. Perform the calculations for
N = 1e1 , 1e2 , 1e3 , and 1e4 , and estimate the error. Compare your
results with the exact value.*/
double function(double x){
    return x*(1-x);
}
double weighted_function(double y){
    return (2.0/PI)*(y*(1-y))/(sin(PI*y));
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
    srand(time(NULL));
    double stepSize = 1/nPoints;
    double* x = (double*)malloc(sizeof(double)*nPoints);
    double* y = (double*)malloc(sizeof(double)*nPoints); 

    double sum = 0;
    for (int i = 0; i < nPoints; i++){
        x[i] = (double) rand()/RAND_MAX;
        y[i] = (1.0/PI)*acos(1.0-2.0*x[i]);

        sum += weighted_function(y[i]);
    }  
    double average = sum/nPoints;
    double exact   = 1.0/6.0;
    double error   = fabs((exact - average)/exact)*100;
    printf("Exact value = %lf \n",exact);
    printf("Approximated value = %lf \n",average);
    printf("relative error = %lf%\n",error);
    write_to_file("weighted.csv",x,y);
    return 0;
}