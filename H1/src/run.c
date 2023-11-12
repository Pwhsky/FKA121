#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lattice.h"
#include "potential.h"
#include "tools.h"
void task1();
void task2(){

    int nAtoms = 256;
    double positions[256][3];
    double a0 = 4.05;

}

int
run(
    int argc,
    char *argv[]
   )
{   
    task1();

	
    return 0;
}

void task1(){
    int nAtoms = 256;
    int n_trials = 10;
    
    double lattice_params[n_trials];
    double potential_energies[n_trials];
    for(int i = 0; i<n_trials; i++){
        lattice_params[i] = 4.05- (n_trials/2.0 - i)*0.05;
    }
    for(int i = 0; i<n_trials; i++){
    
        double positions[256][3];
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