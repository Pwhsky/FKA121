#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#define PI 3.14159265359
#define N 100000

/*Generate N points x_i uniformly on the interval [0,1] and evaluate the
integral with the Monte Carlo method. Perform the calculations for
N = 1e1 , 1e2 , 1e3 , and 1e4 , and estimate the error. Compare your
results with the exact value.*/


void read_data(char *fname, double *array) {
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
        perror("error: ");
        exit(1);
    }
    int i = 0;
    while (i < N && fscanf(fp, "%lf", &array[i]) == 1) {
        i++;
    }
    fclose(fp);
}

void write_to_file(char *fname,double* x, double* y ) {
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "Block_size,estimate\n");
    for (int i = 0; i < N; ++i) {
        fprintf(fp, "%f,%f\n",x[i] ,y[i] );
    }
    fclose(fp);
}
double computeVariance(double* data,int nElements){
    double variance;
    double sumSquares = 0;
    double squareSum = 0;
    
    for(int i = 0; i < nElements; i++){
        sumSquares += pow(data[i],2);
        squareSum  += data[i];
    }
    sumSquares = sumSquares/nElements;
    squareSum  = pow(squareSum/nElements,2);
    variance = sumSquares - squareSum;
    return variance;
}


int main() {
    double* data      = (double*)malloc(sizeof(double)*N);
    read_data("MC.txt",data);
    //First, compute variance from dataset:
    double f = computeVariance(data,N);
    int Block_size = 100;
    size_t nElements  = N/Block_size;
    double* block_averages = (double*)malloc(sizeof(double)*nElements);
    int index = 0;
    for (int i = 0; i < nElements; i+=Block_size){
        double average = 0;
        
        for(int j = i; j< Block_size+i; j++){
            average+=data[j];
        }

        block_averages[index] = average/((double)Block_size);
        index++;
    }

    double F = computeVariance(block_averages,nElements);
    double n_s = Block_size*F/f;
    printf("%lf\n",f);

    return 0;
}