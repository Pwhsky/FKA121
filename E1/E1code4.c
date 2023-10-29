#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_const_mksa.h>
#include "fft.h"
#define  m_asu 9649

void calc_acc(double *a, double *u, double *m, double kappa, int size_of_u)
{
    int i;
    a[0] = kappa * (-2 * u[0] + u[1]) / m[0];
    a[size_of_u - 1] = kappa * (u[size_of_u - 2] - 2 * u[size_of_u - 1]) / m[size_of_u - 1];
    for (i = 1; i < size_of_u - 1; i++) {
        a[i] = kappa * (u[i - 1] - 2 * u[i] + u[i + 1]) / m[i];
    }
}

void velocity_verlet(int n_timesteps, int n_particles, double *v, double *q_1, double *q_2, double *q_3,
                                                                    double *kinetic, double *potential, double *total,
                                                                     double dt, double *m, double kappa){
    double q[n_particles];
    double a[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];
    calc_acc(a, q, m, kappa, n_particles);
    for (int i = 1; i < n_timesteps + 1; i++) {
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * v[j];
        }
        calc_acc(a, q, m, kappa, n_particles);
        for (int j = 0; j < n_particles; j++) {
            v[j] += dt * 0.5 * a[j];
        }
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];
        kinetic[i]    =   m[1]*0.5*( pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
        //potential[i]  = kappa*0.5*( pow(q_1[i],2) + pow(q_2[i],2) + pow(q_3[i],2));

        for (int j = 0; j < n_particles+1; j++) {
            if(j == 0) {
                potential[i] += pow(q[j], 2)*kappa/2.0;
            } else if(j == n_particles){
                potential[i] += pow(q[j-1], 2)*kappa/2.0;
            }else{
                potential[i] += pow(q[j]-q[j-1], 2)*kappa/2.0;
            }
        }
     
            total[i]  = kinetic[i] + potential[i];
        
        
    }
}

void write_to_file(char *fname,double* time, int n_timesteps,  double* q1, double* q2, double* q3,double* kinetic, double* potential, double* total) {
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time,q1,q2,q3\n");
    for (int i = 0; i < n_timesteps; ++i) {
        fprintf(fp, "%f,%f,%f,%f,%f,%f,%f\n", time[i], q1[i], q2[i],q3[i],kinetic[i],potential[i],total[i]);
    }
    fclose(fp);
}

void write_powerspectrum_to_file(double *q1, double *q2, double *q3, 
                              double *frequencies, int n_points){
    FILE *fp = fopen("powerspectrum.csv", "w");
    fprintf(fp, "q1,q2,q3,frequencies\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f,%f,%f\n", q1[i], q2[i], q3[i], frequencies[i]);
    }
    fclose(fp);
}

void write_qs_file(double *q1, double *q2, double *q3, 
                   double *U_kin, double *U_pot, 
                   double *timesteps, int n_points){
    FILE *fp = fopen("timetrail.csv", "w");
    fprintf(fp, "q1,q2,q3,U_kin,U_pot,time\n");
    for(int i = 0; i < n_points; ++i){
	    fprintf(fp, "%f,%f,%f,%f,%f,%f\n", q1[i], q2[i], q3[i], U_kin[i], U_pot[i], timesteps[i]);
    }
    fclose(fp);
}


int main(){
    int n_timesteps = 10000;
    double* time = calloc(sizeof(double),n_timesteps+1);
    int n_particles  = 3;
    double kappa = 1000/16.0218; double dt = 0.00001; double carbon_mass = 12.01/m_asu;
    for(int i = 0; i < n_timesteps+1;i++){
        time[i] = i*dt;
    }
    

    double* q1 = calloc(sizeof(double),n_timesteps+1);
    double* q2 = calloc(sizeof(double),n_timesteps+1);
    double* q3 = calloc(sizeof(double),n_timesteps+1);
    double* kinetic = calloc(sizeof(double),n_timesteps+1);
    double* potential = calloc(sizeof(double),n_timesteps+1);
    double* total = calloc(sizeof(double),n_timesteps+1);
    
    double* m = calloc(sizeof(double),n_particles);
     m[0] = carbon_mass; m[1] = carbon_mass;m[2] = carbon_mass;
  
    double* v = calloc(sizeof(double), n_particles);
    q1[0] = 0.01; q2[0] = 0.005; q3[0] = -0.005;

  
    velocity_verlet(n_timesteps, n_particles, v,q1,q2,q3, kinetic,potential,total, dt, m, kappa);

    write_qs_file(q1, q2, q3, kinetic, potential, time, n_timesteps);
    int N_POINTS = 12000;
    double fftd_q1[ N_POINTS];
    double fftd_q2[ N_POINTS];
    double fftd_q3[ N_POINTS];
    powerspectrum(q1, fftd_q1, N_POINTS); 
    powerspectrum_shift(fftd_q1,  N_POINTS);

    powerspectrum(q2, fftd_q2,  N_POINTS);
    powerspectrum_shift(fftd_q2,  N_POINTS);

    powerspectrum(q3, fftd_q3,  N_POINTS);
    powerspectrum_shift(fftd_q3,  N_POINTS);

    double frequencies[ N_POINTS];
    for(int i = 0; i <  N_POINTS; i++){
	    frequencies[i] = i / (dt * N_POINTS);
    }
    write_to_file("coupled_oscillator.csv",time, n_timesteps,q1,q2,q3,kinetic, potential,total);

    fft_freq_shift(frequencies, dt, N_POINTS);
    write_powerspectrum_to_file(fftd_q1, fftd_q2, fftd_q3, frequencies,  N_POINTS);
    free(q1);
    free(q2);
    free(q3);
    free(m);
    free(time);
    return 0;
}