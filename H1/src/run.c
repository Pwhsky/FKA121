#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Isothermal compressability for aluminium
double kappa_T = 1.0/760000.0;
double kb = 8.6173e-5;

void verlet_step(double** positions, double** velocities, double** forces, int nAtoms, double dt,  double m, double L);
double get_Ekin(double** v, int nAtoms, double m);
double get_temp(double** velocities, int nAtoms, double m);
double get_pressure(double T, double W, double V,int nAtoms);
double get_alphaT(double T_eq, double T, double dt, double tau);
double get_alphaP(double P_eq,double P,double dt,double tau);
double get_MSD(double** positions,double** initial_positions,int nAtoms);

void task1(); //Lattice parameter
void task2(); //material relaxation
void task3() //Molecular dynamics
{
    int timeSteps = 200000;
    
    double dt = 0.001; double T_eq = 973.15; double P_eq =1; double tau_T = 300*dt; double tau_P = 300*dt; 
    int nAtoms = 256; double a0 = 4.03; double cell_length = 4.0*a0; double m = 26.0/9649.0; //aluminum mass


    double** positions = create_2D_array(256,3);
    init_fcc(positions,4,a0);
   
    // Create a random number generator
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default; 
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));
    double lower_bound = -0.065*a0;
    double upper_bound = 0.065*a0;
    double** forces    = create_2D_array(256,3);
    double** velocities = create_2D_array(256,3);
    //Randomly perturb the positions:
    for(int i = 0; i<nAtoms;i++){
        for(int j=0; j < 3;j++){
            double displacement = gsl_ran_flat(r, lower_bound, upper_bound);
            positions[i][j] += displacement;
            velocities[i][j] = 0.0;
        }
    }

    double** sample_trajectory = create_2D_array(timeSteps,3);
    double* temperature      = (double*)malloc(sizeof(double)*timeSteps);
    double* pressure         = (double*)malloc(sizeof(double)*timeSteps);
    double* lattice_constant = (double*)malloc(sizeof(double)*timeSteps);
    double alphaT = 1.0;
    double alphaP = 1.0;
    double V         = cell_length*cell_length*cell_length;
    get_forces_AL(forces,positions,cell_length,nAtoms);
    

    for(int t = 0; t< timeSteps;t++){


        if(t < 1000){
            T_eq = 1200.0;
        }else{
            T_eq = 973.15;}

        for(int i = 0; i < nAtoms; i++){
            for(int j = 0; j < 3; j++){
                velocities[i][j] += dt * 0.5 * forces[i][j]/m;

                positions[i][j] += dt*velocities[i][j];
            }
        }   
        get_forces_AL(forces,positions,cell_length,nAtoms);


        for(int i = 0; i < nAtoms; i++){
            for(int j = 0; j < 3; j++){
                velocities[i][j] += dt*0.5*forces[i][j]/m;
            }
        }   
        //Rescale:

     
        double W         = get_virial_AL(positions, cell_length,nAtoms)/4.0;
        double temp      = get_temp(velocities,nAtoms,m);
        double press     = get_pressure(temp,W,V,nAtoms);
        V                   = cell_length*cell_length*cell_length;
        alphaT              = sqrt(get_alphaT(T_eq,temp,dt,tau_T));
        alphaP              = cbrt(get_alphaP(P_eq,press,dt,tau_P));
        cell_length *= alphaP;

        for(int i = 0; i < nAtoms; i++){
            for(int j = 0; j < 3; j++){
                velocities[i][j] *= alphaT;
            }
        }   
        for(int i = 0; i < nAtoms; i++){
            for(int j = 0; j < 3; j++){
                positions[i][j] *= alphaP;
            }
        }   

        
        //Update alphas


        lattice_constant[t] = cell_length/4.0;
        temperature[t]      = temp;
        pressure[t]         = press;
        sample_trajectory[t][0] = positions[0][0];
        sample_trajectory[t][1] = positions[1][0];
        sample_trajectory[t][2] = positions[2][0];
    }
    
    FILE *fp = fopen("task3.csv", "w");
    fprintf(fp, "time,temperature,pressure,lattice,MSD\n");

    for (int i = 0; i < timeSteps; i++) {
        fprintf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",i*dt,temperature[i],pressure[i],lattice_constant[i],
                                                            sample_trajectory[i][0],sample_trajectory[i][1],sample_trajectory[i][2]);
    }
    fclose(fp);
}

int
run(
    int argc,
    char *argv[]
   )
{   
    //task1();
    //task2();
	task3();
    
    return 0;
}

void task1(){
    int nAtoms = 256;
    int n_trials = 15;
    
    double lattice_params[n_trials];
    double potential_energies[n_trials];
    for(int i = 0; i<n_trials; i++){
        lattice_params[i] = 4.05- (n_trials/2.0 - i)*0.01;
    }

    for(int i = 0; i<n_trials; i++){
    
        double** positions = create_2D_array(256,3);
        double a0 = lattice_params[i];
        double cell_length = 4.0*a0;
        init_fcc(positions,4,a0);
        double Epot = get_energy_AL(positions, cell_length,nAtoms);
        //Divide with total volume (4x4x4)
        Epot = Epot/(4.0*4.0*4.0);
        potential_energies[i] = Epot;

    }

    FILE *fp = fopen("task1.csv", "w");
    fprintf(fp, "Epot,a0\n");

    for (int i = 0; i < n_trials; i++) {
        fprintf(fp, "%lf,%lf\n",potential_energies[i],lattice_params[i]);
    }
    fclose(fp);


}
void task2(){
    int timeSteps = 100;
    
    double dt = 0.1; //0.001 0.01
    int nAtoms = 256; double a0 = 4.04; double cell_length = 4.0*a0; double m = 26.0/9649.0; //aluminum mass
    double** positions = create_2D_array(256,3);
    init_fcc(positions,4,a0);

    // Create a random number generator
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default; 
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));
    double lower_bound = -0.065*a0;
    double upper_bound = 0.065*a0;
    double** forces    = create_2D_array(256,3);
    double** velocities = create_2D_array(256,3);
    //Randomly perturb the positions:
    for(int i = 0; i<nAtoms;i++){
        for(int j=0; j < 3;j++){
            double displacement = gsl_ran_flat(r, lower_bound, upper_bound);
            positions[i][j] += displacement;
            velocities[i][j] = 0.0;
        }
    }
    double* kinetic_energy   = (double*)malloc(sizeof(double)*timeSteps);
    double* potential_energy = (double*)malloc(sizeof(double)*timeSteps);
    double* total_energy     = (double*)malloc(sizeof(double)*timeSteps);
    double* temperature      = (double*)malloc(sizeof(double)*timeSteps);

    double* time = (double*)malloc(sizeof(double)*timeSteps);

    get_forces_AL(forces,positions,cell_length,nAtoms);
    for(int t = 0; t< timeSteps;t++){
        time[t] = t*dt; 
        verlet_step(positions,velocities,forces, nAtoms, dt, m,  cell_length);
        kinetic_energy[t]   = get_Ekin(velocities,nAtoms,m);
        potential_energy[t] = get_energy_AL(positions,cell_length,nAtoms);
        total_energy[t]     = kinetic_energy[t]+potential_energy[t];
        temperature[t]      = get_temp(velocities,nAtoms,m);
    }
    
    FILE *fp = fopen("task2bigdt.csv", "w");
    fprintf(fp, "time,Ekin,Epot,Etot,temp\n");

    for (int i = 0; i < timeSteps; i++) {
        fprintf(fp, "%lf,%lf,%lf,%lf,%lf\n",time[i],kinetic_energy[i],potential_energy[i],total_energy[i],temperature[i]);
    }
    fclose(fp);

}

void verlet_step(double** positions, double** velocities, double** forces, int nAtoms, double dt,  double m, double L){

    for(int i = 0; i < nAtoms; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] += dt*forces[i][j]/(2.0*m);
            positions[i][j] += dt*velocities[i][j];
        }
    }

    get_forces_AL(forces,positions, L, nAtoms);

    for(int i = 0; i < nAtoms; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] += dt*forces[i][j]/(2.0*m);
        }
    }
}
void verlet_stepT3(double** positions, double** velocities, double** forces, int nAtoms, double dt,  double m, double L){

    for(int i = 0; i < nAtoms; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] += dt*forces[i][j]/(2.0*m);
            positions[i][j] += dt*velocities[i][j];
        }
    }

    get_forces_AL(forces,positions, L, nAtoms);

    for(int i = 0; i < nAtoms; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] += dt*forces[i][j]/(2.0*m);
        }
    }
}


double get_Ekin(double** velocities, int nAtoms, double m){
    double energy = 0.0;
    for(int i = 0; i < nAtoms; i++){
        for(int j = 0; j < 3; j++){
            energy += m*pow(velocities[i][j],2)/2.0;
        }
    }  
    return energy;
}
double get_temp(double** velocities, int nAtoms, double m){
    double energy = 0.0;
    for(int i = 0; i < nAtoms; i++){
        for(int j = 0; j < 3; j++){
            energy += m*pow(velocities[i][j],2)/2.0;
        }
    }  

    energy = energy * 2.0/((double)3.0*kb*nAtoms); //convert to kelvin
    return energy;
}
double get_alphaT(double T_eq, double T,double dt,double tau){

    return 1 + (2*dt/(tau))*(T_eq-T)/(T); 
}

double get_alphaP(double P_eq, double P, double dt, double tau){

    return 1 - (kappa_T*dt/(tau))*(P_eq-P); 
}
double get_pressure(double T, double W, double V,int nAtoms){

    return ((256*kb*T+W*1602000.0)/V);
}
