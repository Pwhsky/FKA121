#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
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

void velocity_verlet(int n_timesteps, int n_particles, double *v, double *q_1, double *q_2, double *q_3, double dt, double *m, double kappa)
{
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
    }
}

void write_to_file(char *fname,double* time, int n_timesteps,  double* q1, double* q2, double* q3) {
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time,q1,q2,q3\n");
    for (int i = 0; i < n_timesteps; ++i) {
        fprintf(fp, "%f,%f,%f,%f\n", time[i], q1[i], q2[i],q3[i]);
    }
    fclose(fp);
}

int main(){
    int n_timesteps = 1500;
    double* time = calloc(sizeof(double),n_timesteps+1);
    int n_particles  = 3;
    double kappa = 1000/16.0218; double dt = 0.0001; double carbon_mass = 12.01/m_asu;
    for(int i = 0; i < n_timesteps+1;i++){
        time[i] = i*dt;
    }
    

    double* q1 = calloc(sizeof(double),n_timesteps+1);
    double* q2 = calloc(sizeof(double),n_timesteps+1);
    double* q3 = calloc(sizeof(double),n_timesteps+1);
    double* m = calloc(sizeof(double),n_particles);
     m[0] = carbon_mass; m[1] = carbon_mass;m[2] = carbon_mass;
  
    double* v = calloc(sizeof(double), n_particles);
    q1[0] = 0.01; q2[0] = 0.005; q3[0] = -0.005;

  
    velocity_verlet(n_timesteps, n_particles, v,q1,q2,q3, dt, m, kappa);


    write_to_file("coupled_oscillator.csv",time, n_timesteps,q1,q2,q3) ;

    free(q1);
    free(q2);
    free(q3);
    free(m);
    free(time);
    return 0;
}