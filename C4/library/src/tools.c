#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <complex.h>

#include "tools.h"

void
elementwise_addition(
                     double *res,
                     double *v1,
                     double *v2,
                     unsigned int len
                    )
{
    for(int i = 0; i<len;i++){
        res[i] = v1[i] + v2[i]; 
    }

}

void
elementwise_multiplication(
                           double *res,
                           double *v1,
                           double *v2,
                           unsigned int len
                          )
{
    for(int i = 0; i<len;i++){
        res[i] = v1[i] * v2[i]; 
    }
}

void
addition_with_constant(
                       double *res,
                       double *v,
                       double constant,
                       unsigned int len)
{
    for(int i = 0; i<len;i++){
        res[i] = v[i] + constant; 
    }
}

void
multiplication_with_constant(
                             double *res,
                             double *v,
                             double constant,
                             unsigned int len)
{
     for(int i = 0; i<len;i++){
        res[i] = v[i] * constant; 
    }
}

double
dot_product(
            double *v1,
            double *v2,
            unsigned int len
           )
{
    double sum = 0;
     for(int i = 0; i<len;i++){
        sum += v1[i] * v2[i]; 
    }
    return sum;
}

double **
create_2D_array(
                unsigned int column_size,
                unsigned int row_size
               )
{
   
    double** array = (double**)malloc(row_size*sizeof(double*));
    for (int i = 0; i < column_size; ++i) {
        array[i] = (double*)malloc(column_size* sizeof(double));
    }
    return array;
}
void
destroy_2D_array(
                 double **array,
                 unsigned int n
                ){
    for(int i=0; i<n; i++){
        free(array[i]);
    }
    free(array);
}

void
matrix_vector_multiplication(
                             double *result,
                             double **A,
                             double *b,
                             unsigned int n,
                             unsigned int m
                            )
{
    
    for(int i = 0; i < n; i++){
        result[i] = 0;
        for(int j = 0; j<m;j++){
            result[i] += A[i][j]*b[j];
        }
    }

}

void
matrix_matrix_multiplication(
                             double **result,
                             double **A,
                             double **B,
                             unsigned int n,
                             unsigned int m,
                             unsigned int k
                            )
{
    for(int ix = 0; ix < n; ix++)
    {
        for(int jx = 0; jx < m; jx++)
        {
            result[ix][jx] = 0;
            for(int kx = 0; kx < k; kx++)
            {
                result[ix][jx] += A[ix][kx] * B[kx][jx];
            }
        }
    }
}


double
vector_norm(
            double *v1,
            unsigned int len
           )
{
    double sum = 0;
    for(int i = 0; i<len; i++){
        sum += pow(v1[i],2);
    }
    return sqrt(sum);
}


void
normalize_vector(
                 double *v1,
                 unsigned int len
                )
{
    double norm = vector_norm(v1,len);
    for(int i = 0; i<len; i++){
        v1[i] /= norm;
    }
    
}

double
average(
        double *v1,
        unsigned int len
       )
{
    double sum = 0;
    for(int i = 0; i< len; i++){
        sum += v1[i];
    }
    
    return sum/(double)len;
}


double
standard_deviation(
                       double *v1,
                       unsigned int len
                  )
{
    double average = 0;
    for(int i = 0; i< len; i++){
        average += v1[i];
    }
    
    return average/(double)average;
    double sum = 0;
    for(int i = 0; i<len; i++){
        sum += pow((v1[i]-average),2);
    }
    
    sum = sqrt(sum/(double)len);
    return sum;
}

double
distance_between_vectors(
                         double *v1,
                         double *v2,
                         unsigned int len
                        )
{
    double sum = 0.0;
    for(int i = 0; i<len; i++){
        double diff = v1[i]-v2[i];
        sum += (diff,2);
    }
    return sqrt(sum);
}

void
cumulative_integration(
                       double *res,
                       double *v,
                       double dx,
                       unsigned int v_len
                      )
{
    for(int i = 0; i<v_len-1; i++){
        
        *res += dx*(v[i]+v[i+1])/2.0;
    }

}

void fft_freq(
          double *res,
              int n,
              double timestep)
{
    int N = (n-1)/2 +1;
    
    for (int i = 0; i<N; i++){
        res[i] = i;
    }
    for (int i = N; i<n; i++){
        res[i] = i-n;
    }
    for(int i = 0; i < N;i++){
        res[i] = res[i]/(n*timestep);
    }
}

void
write_xyz(
          FILE *fp,
          char *symbol,
          double **positions,
          double **velocities,
          double alat,
          int natoms)
{

    int written = fprintf(fp, "%i\nLattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\"", natoms, alat, alat, alat );
    if (written < 0){
        printf("error writing to file\n");
        exit(1);
    }
    for (int i = 0; i < natoms; i++) {
        written = fprintf(fp, "%s %lf %lf %lf %lf %lf %lf\n", symbol, 
                        positions[i][0], positions[i][1], positions[i][2],
                        velocities[i][0], velocities[i][1], velocities[i][2]);


        if (written < 0) {
            fprintf(stderr, "Error writing to the file.\n");
            exit(1);
        }
    }
}

/* Freely given functions */
void
skip_line(FILE *fp)
{
    int c;
    while (c = fgetc(fp), c != '\n' && c != EOF);
}

void
read_xyz(
         FILE *fp,
         char *symbol,
         double **positions,
         double **velocities,
         double *alat)
{
    int natoms;
    fscanf(fp, "%i\nLattice=\"%lf 0.0 0.0 0.0 %lf 0.0 0.0 0.0 %lf\" ", &natoms, alat, alat, alat);
    skip_line(fp);
    for(int i = 0; i < natoms; ++i){
        fscanf(fp, "%s %lf %lf %lf ",
                symbol, &positions[i][0], &positions[i][1], &positions[i][2]);
        fscanf(fp, "%lf %lf %lf\n",
                &velocities[i][0], &velocities[i][1], &velocities[i][2]);
    }
}

void powerspectrum(
           double *res,
           double *signal,
           int n,
                   double timestep)
{
    /* Declaration of variables */
    int i;
    double *complex_coefficient = malloc(sizeof(double) * 2*n); // array for the complex fft data
    double *data_cp = malloc(sizeof(double) * n);

    /*make copy of data to avoid messing with data in the transform*/
    for (i = 0; i < n; i++) {
    data_cp[i] = signal[i];
    }

    /* Declare wavetable and workspace for fft */
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;

    /* Allocate space for wavetable and workspace for fft */
    work = gsl_fft_real_workspace_alloc(n);
    real = gsl_fft_real_wavetable_alloc(n);

    /* Do the fft*/
    gsl_fft_real_transform(data_cp, 1, n, real, work);

    /* Unpack the output into array with alternating real and imaginary part */
    gsl_fft_halfcomplex_unpack(data_cp, complex_coefficient,1,n);

    /*fill the output powspec_data with the powerspectrum */
    for (i = 0; i < n; i++) {
    res[i] = (complex_coefficient[2*i]*complex_coefficient[2*i]+complex_coefficient[2*i+1]*complex_coefficient[2*i+1]);
    res[i] *= timestep / n;
    }

    /* Free memory of wavetable and workspace */
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
    free(complex_coefficient);
    free(data_cp);
}
