#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "tools.h"
#include <math.h>
#define N 2000
 int burnoff = 1000;
 const double T = 300.0; 
 const double kb = 1.380649e-8; 
 const double m = 3.0134e-5;
 const double tau_high = 147.3e-3; //microseconds
 const double tau_low  = 48.5e-3; //microseconds
 const double w0 = 19.477; //rad per ms
 const double dt = 0.001;

typedef struct {
    double position;
    double velocity; 
} result_t;

result_t BD3(double initial_position, double initial_velocity, double dt, double eta, gsl_rng *k)
{
    result_t result;

    double c0 = exp(-eta*dt);
    double v_thermal  = sqrt(kb*T/m);
    double variance = 2.0*eta*kb*T/m;
    //Compute initial acceleration:
    double a = (-w0*w0 * initial_position);
    double new_velocity = 0.5*a*dt + 
                        sqrt(c0)*initial_velocity + 
                        v_thermal*sqrt(1-c0)*gsl_ran_gaussian(k,1.0);
    
    double new_position = initial_position + new_velocity*dt;
    double new_a = (-w0*w0 * new_position);

    new_velocity = 0.5*new_a*dt*sqrt(c0) + 
                        sqrt(c0)*new_velocity + 
                        v_thermal*sqrt(1.0-c0)*gsl_ran_gaussian(k,1.0);
                        
    result.position = new_position;
    result.velocity = new_velocity;
    return result;
}
void write_to_file(char *fname,double** array){
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "x,v,t\n");
    for (int i = 0; i < N; i++) {
       
        fprintf(fp, "%lf, %lf, %lf\n",array[i][0],array[i][1],array[i][2]);
    }
    fclose(fp);
}
int main(){
    const gsl_rng_type *T;
    gsl_rng *k;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    k = gsl_rng_alloc(T);
    gsl_rng_set(k, time(NULL));

    double** trajectory = create_2D_array(N,3);

    double eta = 1/tau_high;
    double init_pos = 0.0;
    double init_vel = 0.0;
    double a = -w0*w0*init_pos;
    result_t step = BD3(init_pos,init_vel,dt,eta,k);

    //brunoff 
    for (int i = 0; i< burnoff; i++){
        
        step = BD3(step.position,step.velocity,dt,eta,k);
        
    }

    //production run:
    for(int i = 0; i<N; i++){

        trajectory[i][0] = step.position;
        trajectory[i][1] = step.velocity;
        trajectory[i][2] = i*dt;
        step = BD3(step.position,step.velocity,dt,eta,k);
    }
    write_to_file("T3.csv",trajectory);


    return 0;
}