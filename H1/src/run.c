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
double avgTemp = 0.0;
int nAtoms = 256;

double* temperature;
double* pressure;
double* lattice_constant;
double** sample_trajectory; double** positions;
double** velocities; double** forces;
double* potential_energy; double* kinetic_energy; double* total_energy;
gsl_rng *r;
const gsl_rng_type *T;
int timeSteps;
void verlet_step(double** positions, double** velocities, double** forces,
                 int nAtoms, double dt,  double m, double L);
                 
double get_Ekin(double** v, int nAtoms, double m);
double get_temp(double** velocities, int nAtoms, double m);
double get_pressure(double T, double W, double V,int nAtoms);
double get_alphaT(double T_eq, double T, double dt, double tau);
double get_alphaP(double P_eq,double P,double dt,double tau);
void perturb(double** positions,double** velocities);
   
double  kinetic= 0, squaredAverage= 0,  averageSquared= 0,alphaT = 1.0, 
alphaP = 1.0; int timeSteps;
int burnoff = 5000;
double lower_bound = -0.065*4.03;
double upper_bound = 0.065*4.03;

void lattice_parameter_simulation(); //Lattice parameter
void energy_conservation_simulation(); //material relaxation
void solid_aluminum_simulation(); //Pressure and temperature for solid  AND heat capacity
void liquid_aluminum_simulation(); //Pressure and temperature for liquid AND heat capacity AND radial distribution function

int
run(
    int argc,
    char *argv[]
   )
{      
     // Create a random number generator
    
    gsl_rng_env_setup();
    T = gsl_rng_default; 
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, time(NULL));
    timeSteps = 20000;
    temperature      = (double*)malloc(sizeof(double)*timeSteps); 
    pressure         = (double*)malloc(sizeof(double)*timeSteps);
    lattice_constant = (double*)malloc(sizeof(double)*timeSteps);
    kinetic_energy   = (double*)malloc(sizeof(double)*timeSteps);
    potential_energy = (double*)malloc(sizeof(double)*timeSteps);
    total_energy     = (double*)malloc(sizeof(double)*timeSteps);

    positions  = create_2D_array(256,3);
    velocities = create_2D_array(256,3);
    forces     = create_2D_array(256,3);
    sample_trajectory = create_2D_array(timeSteps,3);

    //lattice_parameter_simulation();
    //energy_conservation_simulation();
	solid_aluminum_simulation();
    liquid_aluminum_simulation();
    
    return 0;  
}


void solid_aluminum_simulation() //Molecular dynamics
{
    double factor =3.0*nAtoms*kb/2.0;
    double dt = 0.001; double T_eq = 773.15; double P_eq = 1.0; 
    double tau_T = 100*dt; double tau_P = 300*dt; 
    double a0 = 4.04; double cell_length = 4.0*a0; double m = 26.0/9649.0; 

    init_fcc(positions,4,a0);
    perturb(positions,velocities);

    double V  = cell_length*cell_length*cell_length;


    double sum = 0.0;
    double pressureSum = 0.0;
    avgTemp = 0.0;
    get_forces_AL(forces,positions,cell_length,nAtoms);
    for(int t = 0; t< timeSteps;t++){

        verlet_step(positions, velocities, forces,nAtoms, dt, m, cell_length);
        V                = cell_length*cell_length*cell_length;
        double W         = get_virial_AL(positions, cell_length,nAtoms);
        double temp      = get_temp(velocities,nAtoms,m);
        double press     = get_pressure(temp,W,V,nAtoms);
        
        alphaT           = sqrt(get_alphaT(T_eq,temp,dt,tau_T));
        alphaP           = cbrt(get_alphaP(P_eq,press,dt,tau_P));
        cell_length *= alphaP;

        for(int i = 0; i < nAtoms; i++){
            for(int j = 0; j < 3; j++){
                velocities[i][j] *= alphaT;
                positions[i][j] *= alphaP;
            }
        }   
        //Update alphas
        lattice_constant[t] = cell_length/4.0;
        temperature[t]      = temp;
        pressure[t]         = press;
        for (int i = 0; i<3;i++)
            sample_trajectory[t][i]=positions[0][i];
        
      //Compute fluctuations and average pressure after equilibration:
      if (t > burnoff){
            double mean_kinetic = get_Ekin(velocities,nAtoms,m)/(double)nAtoms;
            
            for(int i = 0; i < nAtoms; i++){
                double atom_kinetic = 0.0;
                for(int j = 0; j < 3; j++)
                    atom_kinetic += m*pow(velocities[i][j],2)/2.0;

                double fluctuation = mean_kinetic-atom_kinetic;
                sum += fluctuation*fluctuation/(double)nAtoms;
            }
            pressureSum += press;
            avgTemp += temp;

        }
       
    }


    sum /= ((double)(timeSteps-burnoff));
    double heatCapacity = factor/
    (1- 2.0*sum/(3.0*nAtoms*kb*kb*avgTemp*avgTemp));

    avgTemp /= ((double)timeSteps-burnoff);
    pressureSum /= (double)(timeSteps-burnoff);
    printf("Average temperature during production run = %lf \n Average pressure during production run = %lf \n",avgTemp,pressureSum);
    printf("Heat Capacity for 773.15K = %lf \n",heatCapacity);



   //writing and avg temp computation
    FILE *fp = fopen("solid_aluminum_simulation.csv", "w");
    fprintf(fp, "time,temperature,pressure,lattice,pos1,pos2\n");
    
    for (int i = 0; i < timeSteps; i++){
        fprintf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",i*dt,temperature[i],
        pressure[i],lattice_constant[i],sample_trajectory[i][0],
        sample_trajectory[i][1],sample_trajectory[i][2]);
    } 
    fclose(fp);
}
void liquid_aluminum_simulation() //Molecular dynamics
{
    double factor =3.0*nAtoms*kb/2.0;
    double dt = 0.001; double T_eq = 973.15; double P_eq =1.0; 
    double tau_T = 100*dt; double tau_P = 300*dt; 
    double a0 = 4.04;  double cell_length = 4.0*a0; double m = 26.0/9649.0;
    init_fcc(positions,4,a0);

    perturb(positions,velocities);
    

    double V                = cell_length*cell_length*cell_length;
    double sum = 0.0;
    double pressureSum = 0.0;
    avgTemp = 0.0;
    get_forces_AL(forces,positions,cell_length,nAtoms);
    for(int t = 0; t< timeSteps;t++){

        //Melting routine
        if(t < 1000){ T_eq = 2000.0;}
        else if (t<2000){T_eq -= 0.82685;}
        else{ T_eq = 973.15;}

        verlet_step(positions, velocities, forces,nAtoms, dt, m, cell_length);
        V                = cell_length*cell_length*cell_length;
        double W         = get_virial_AL(positions, cell_length,nAtoms);
        double temp      = get_temp(velocities,nAtoms,m);
        double press     = get_pressure(temp,W,V,nAtoms);
   
        alphaT           = sqrt(get_alphaT(T_eq,temp,dt,tau_T));
        alphaP           = cbrt(get_alphaP(P_eq,press,dt,tau_P));
        cell_length *= alphaP;

        for(int i = 0; i < nAtoms; i++){
            for(int j = 0; j < 3; j++){
                velocities[i][j] *= alphaT;
                positions[i][j]  *= alphaP;
            }
        }   
        lattice_constant[t] = cell_length/4.0;
        temperature[t]      = temp;
        pressure[t]         = press;
        for (int i = 0; i<3;i++)
            sample_trajectory[t][i]=positions[0][i];
        //Compute fluctuations:
        if (t > burnoff){

            double mean_kinetic = get_Ekin(velocities,nAtoms,m)/(double)nAtoms;
            
            for(int i = 0; i < nAtoms; i++){
                double atom_kinetic = 0.0;
                for(int j = 0; j < 3; j++)
                    atom_kinetic += m*pow(velocities[i][j],2)/2.0;

                double fluctuation = mean_kinetic-atom_kinetic;
                sum += fluctuation*fluctuation/(double)nAtoms;
            }
            pressureSum += press;
            avgTemp += temp;

        }
        
    }
    
    FILE *fp = fopen("liquid_aluminum_simulation.csv", "w");
    fprintf(fp, "time,temperature,pressure,lattice,pos1,pos2\n");
    for (int i = 0; i < timeSteps; i++) {
        fprintf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",i*dt,
    temperature[i],pressure[i], lattice_constant[i],
      sample_trajectory[i][0],sample_trajectory[i][1],sample_trajectory[i][2]);

    }
    fclose(fp);
    sum /= ((double)(timeSteps-burnoff));
    avgTemp /= ((double)timeSteps-burnoff);
    double heatCapacity = factor/
    (1- 2.0*sum/(3.0*nAtoms*kb*kb*avgTemp*avgTemp));


     
    pressureSum /= (double)(timeSteps-burnoff);
    printf("Average temperature during production run = %lf \n Average pressure during production run = %lf \n",avgTemp,pressureSum);


    printf("Heat Capacity for 973.15K = %lf \n",heatCapacity);
}


void lattice_parameter_simulation(){

    int n_trials = 15;
    double lattice_params[n_trials], potential_energies[n_trials];
    
    for(int i = 0; i<n_trials; i++)
        lattice_params[i] = 4.05- (n_trials/2.0 - i)*0.01;

    for(int i = 0; i<n_trials; i++){
    
        double a0 = lattice_params[i];
        double cell_length = 4.0*a0;
        init_fcc(positions,4,a0);
        double Epot = get_energy_AL(positions, cell_length,nAtoms);
        //Divide with total volume (4x4x4)
        Epot = Epot/(4.0*4.0*4.0);
        potential_energies[i] = Epot;

    }

    FILE *fp = fopen("lattice_parameter_simulation.csv", "w");
    fprintf(fp, "Epot,a0\n");

    for (int i = 0; i < n_trials; i++) 
        fprintf(fp, "%lf,%lf\n",potential_energies[i],lattice_params[i]);
    fclose(fp);  


}
void energy_conservation_simulation(){
    
    double dt = 0.1; //0.001 0.01
    double a0 = 4.04; double cell_length = 4.0*a0; double m = 26.0/9649.0;
    init_fcc(positions,4,a0);
    perturb(positions,velocities);

    get_forces_AL(forces,positions,cell_length,nAtoms);
    for(int t = 0; t< timeSteps;t++){
 
        verlet_step(positions,velocities,forces, nAtoms, dt, m,  cell_length);
        kinetic_energy[t]   = get_Ekin(velocities,nAtoms,m);
        potential_energy[t] = get_energy_AL(positions,cell_length,nAtoms);
        total_energy[t]     = kinetic_energy[t]+potential_energy[t];
        temperature[t]      = get_temp(velocities,nAtoms,m);
    }
    
    //Writing
    FILE *fp = fopen("energy_conservation_simulationbigdt.csv", "w");
    fprintf(fp, "time,Ekin,Epot,Etot,temp\n");
    for (int i = 0; i < timeSteps; i++) 
        fprintf(fp, "%lf,%lf,%lf,%lf,%lf\n",i*dt,kinetic_energy[i],
                    potential_energy[i],total_energy[i],temperature[i]);
    fclose(fp);

}
void verlet_step(double** positions, double** velocities, double** forces,
                 int nAtoms, double dt,  double m, double L){
    for(int i = 0; i < nAtoms; i++){
        for(int j = 0; j < 3; j++){
            velocities[i][j] += dt*forces[i][j]/(2.0*m);
            positions[i][j] += dt*velocities[i][j];
        }
    }
    get_forces_AL(forces,positions, L, nAtoms);
    for(int i = 0; i < nAtoms; i++)
        for(int j = 0; j < 3; j++)
            velocities[i][j] += dt*forces[i][j]/(2.0*m);

}

double get_Ekin(double** velocities, int nAtoms, double m){
    double energy = 0.0;
    for(int i = 0; i < nAtoms; i++)
        for(int j = 0; j < 3; j++)
            energy += m*pow(velocities[i][j],2)/2.0;
    return energy;
}
double get_temp(double** velocities, int nAtoms, double m){
    double energy = get_Ekin(velocities,nAtoms,m);
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
    return ((nAtoms*kb*T+W*1602000)/V);
}
void perturb(double** positions,double**velocities){
    for(int i = 0; i<nAtoms;i++){
        for(int j=0; j < 3;j++){
            positions[i][j] += gsl_ran_flat(r, lower_bound, upper_bound);
            velocities[i][j] = 0.0;
        }
    }
}