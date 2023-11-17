#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_const_mksa.h>

#define  m_asu 9649
double calculate_potential_energy(double *positions, double kappa);
double calculate_kinetic_energy(double *velocities, double *masses);

void calc_acc(double *accelerations, double *displacements, double *m, double kappa)
{
    int N = 3;
    int i;
    accelerations[0] = kappa * ( displacements[1] - displacements[0]) / m[0];
    accelerations[1] = kappa * ( displacements[2] - 2*displacements[1]+displacements[0]) / m[1];
    accelerations[2] = kappa * (-displacements[2] + displacements[1])/m[2];
    
}

void velocity_verlet_multiple_steps(int n_timesteps, int n_particles, double *velocities,
                                     double *q_1, double *q_2, double *q_3, double dt, double *m, double kappa,double* Ekin,double* Epot){
    double q[n_particles];
    double accelerations[n_particles];
    q[0] = q_1[0];
    q[1] = q_2[0];
    q[2] = q_3[0];
    calc_acc(accelerations, q, m, kappa);
    for (int i = 1; i < n_timesteps + 1; i++) {
    
        for (int j = 0; j < n_particles; j++) {
            velocities[j] += dt * 0.5 * accelerations[j];
        }
        
        for (int j = 0; j < n_particles; j++) {
            q[j] += dt * velocities[j];
        }
        
        calc_acc(accelerations, q, m, kappa);
        for (int j = 0; j < n_particles; j++) {
            velocities[j] += dt * 0.5 * accelerations[j];
        }
        
        q_1[i] = q[0];
        q_2[i] = q[1];
        q_3[i] = q[2];
        Ekin[i] = calculate_kinetic_energy(velocities,m);
        Epot[i] = calculate_potential_energy(q,kappa);
    }
}

void write_to_file(char *fname,double* time, int n_timesteps,  double* q1, double* q2, double* q3,double* Ekin, double* Epot, double* Etot) {
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "time,q1,q2,q3,Ekin,Epot,Etot\n");
    for (int i = 0; i < n_timesteps; ++i) {
        fprintf(fp, "%f,%f,%f,%f,%f,%f,%f\n", time[i], q1[i], q2[i],q3[i],Ekin[i],Epot[i],Etot[i]);
    }
    fclose(fp);
}



int main(){
    int n_timesteps = 10000;
    double* time = calloc(sizeof(double),n_timesteps+1);
    int n_particles  = 3;
    double kappa = 1600/16.0218; double dt = 0.0001; double carbon_mass = 12.01/m_asu; double oxygen_mass = 15.999/m_asu;
    for(int i = 0; i < n_timesteps+1;i++){
        time[i] = i*dt;
    }
    

    double* q1 = calloc(sizeof(double),n_timesteps);
    double* q2 = calloc(sizeof(double),n_timesteps);
    double* q3 = calloc(sizeof(double),n_timesteps);
    double* Ekin = calloc(sizeof(double),n_timesteps);
    double* Epot = calloc(sizeof(double),n_timesteps);    
    double* Etot =  calloc(sizeof(double),n_timesteps);  
    double* m = calloc(sizeof(double),n_particles);
    double* v = calloc(sizeof(double), n_particles);
    m[0] = oxygen_mass; m[1] = carbon_mass;m[2] = oxygen_mass;
  

    q1[0] = 0.1; q2[0] = 0.005; q3[0] = -0.005;

  
    velocity_verlet_multiple_steps(n_timesteps, n_particles, v,q1,q2,q3, dt, m, kappa,Ekin,Epot);
    for(int i = 0; i<n_timesteps+1;i++){
        Etot[i] = Ekin[i] + Epot[i];
    }

    write_to_file("carbon_dioxide.csv",time, n_timesteps,q1,q2,q3,Ekin,Epot,Etot);

   
    free(q1);
    free(q2);
    free(q3);
    free(m);
    free(time);
    return 0;
}

double calculate_potential_energy(double *positions, double kappa)
{
    const unsigned int N = 3;  // There are three particles
    double Epot =0.0;
    Epot += kappa*pow((positions[1]-positions[0]),2)/2.0;
    Epot += kappa*pow((positions[2]-positions[1]),2)/2.0;
    return Epot;
    
}

double calculate_kinetic_energy(double *velocities, double *masses)
{
    const unsigned int N = 3;  // There are three particles
    double Ekin =0.0;
    Ekin += masses[0]*pow((velocities[0]),2)/2.0;
    Ekin += masses[1]*pow((velocities[1]),2)/2.0;
    Ekin += masses[2]*pow((velocities[2]),2)/2.0;
    return Ekin;
    
}
