#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#define PI 3.14159265359
#define N 1000000

int nPoints = 100;
/*Generate N points x_i uniformly on the interval [0,1] and evaluate the
integral with the Monte Carlo method. Perform the calculations for
N = 1e1 , 1e2 , 1e3 , and 1e4 , and estimate the error. Compare your
results with the exact value.*/
double average(double* data, int nElements);
void read_data(char *fname, double *array); 
double variance(double *v1, unsigned int len);




double block_average(double *data,
	                 int block_size,
	                 int data_len
	                )
{
    if (block_size == 0){
        block_size++;
    }
    int lenPerBlock = data_len/(block_size);
    double F[lenPerBlock];
    for (int i=0; i<lenPerBlock; i++)
        F[i] = 0;
    
    for (int i=0; i<lenPerBlock; i++){
        for (int k=0; k<block_size; k++){
            F[i] += data[k+i*block_size]/((double)block_size);
        }
    }
    double variance_F = 0;
    double mean_F = 0;
    for (int i=0; i<lenPerBlock; i++){
        variance_F += F[i]*F[i]/((double)lenPerBlock);
        mean_F += F[i]/((double)lenPerBlock);
    }
    variance_F = variance_F - mean_F*mean_F;
    double var = variance(data,data_len);

    double s = (double)block_size * variance_F / var;
    return s;
}

double autocorrelation(double *data, int time_lag_ind, int max_shift_ind) {
    
    double correlation = 0.0;
    double phi = 0.0;
    double mean = average(data,N);

    for(int j = 0; j<max_shift_ind;j++){
        correlation += data[j+time_lag_ind]*data[j]/(double)(max_shift_ind);
    }
    
    phi = correlation/(variance(data,N));


    return phi;
}


void write_to_file(char *fname,double* x, double* y ) {
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "Block_size,estimate\n");
    for (int i = 0; i < nPoints; ++i) {
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
    double* data         = (double*)malloc(sizeof(double)*N);
    double* block_values = (double*)malloc(sizeof(double)*nPoints);
    double* block_sizes  = (double*)malloc(sizeof(double)*nPoints);

    double* correlation_values = (double*)malloc(sizeof(double)*nPoints);
    double* correlation_lags    = (double*)malloc(sizeof(double)*nPoints);

    read_data("MC.txt",data);
    double mean = average(data,N);
    //shift data to mean
    for(int i = 0; i< N; i++){
        data[i] -= mean;
    }
    //compute s for different block sizes:
    int block_size; int lag;
    for(int i = 0; i< nPoints; i++){
        block_size = i*50;
        lag = i*10;

        correlation_lags[i] = lag;
        block_sizes[i] = block_size;
        correlation_values[i] = autocorrelation(data,lag,N-lag);
        block_values[i]       = block_average(data,block_size,N);
    }   
    write_to_file("block_averaging.csv",block_sizes,block_values);
    write_to_file("correlation.csv",correlation_lags,correlation_values);

    //Determine statistical inefficency by finding the value for exp(-2)
    double s = exp(-2.0);
    int trial = 0;
    for(int i = 0; i<nPoints;i++){
        if ( sqrt((s-correlation_values[i])*(s-correlation_values[i])) < sqrt((s-correlation_values[trial])*(s-correlation_values[trial] ))){
            trial = i;
        }
    }
    double statistical_inefficency = correlation_lags[trial];
    printf("Statistical inefficency for autocorrelation function = %lf",statistical_inefficency);


    

    return 0;
}

double average(double* data, int nElements){
    double sum= 0.0;
    for(int i = 0; i< nElements; i++){
            sum+= data[i];
    }
    sum/=(double)nElements;
    return sum;
}

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


double variance(
                       double *v1,
                       unsigned int len
                  )
{
    double mean = average(v1,len);
    double sum = 0;
    for(int i = 0; i<len; i++){
        sum += pow((v1[i]-mean),2);
    }
    
    sum = sum/(double)len;
    return sum;
}
