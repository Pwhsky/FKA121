#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#define PI 3.14159


double* generate_signal(double* signal, double*t, size_t len_t, double a, double f, double phi){
    for(int i; i< len_t; i++){
        signal[i] = a * cos(2*PI*f*t[i] + 0) + a*cos(2*PI*6*t[i]);
    }
    return signal;
}

void arange(double* array, double start, int len_t, double dt){
    for(int i = 0; i < len_t; i++){
        array[i] = start + i*dt;
    }
}

void write_to_file ( char* fname , double* time_array ,double* signal , int n_points ){
    FILE* fp = fopen ( fname , "w" );
    fprintf (fp ,"time,signal\n" );
    for ( int i = 0; i < n_points ; ++ i ){
        fprintf (fp,"%f,%f\n" , time_array [i] , signal [i]) ;
    }
    fclose(fp) ;
}

int main (){
    int    N = 250; 
    double dt = 0.1;
    double a = 1; 
    double f = 2; 
    double phi = 0;
    double time_array[N];
    arange(time_array,0,N,dt);
    double signal[N];
    generate_signal(signal,time_array,N,a,f,phi) ;
    write_to_file("signal.csv",time_array,signal,N) ;
    return 0;
}



