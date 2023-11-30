#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <tools.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#define kb 8.6173303e-5
#define PI 3.141592653589
double E_AA=-436e-3; 
double E_BB=-113e-3;
double E_AB=-294e-3;

void write_to_file(char *fname,double* array,double dT){
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "temp, val\n");
    for (int i = 0; i < 1e4-1; i++) {
        fprintf(fp, "%lf, %lf\n",i*dT,array[i]);
    }
    fclose(fp);
}


void task1();
int
run(
    int argc,
    char *argv[]
   )
{
    task1();
    //task2();
    return 0;
}

void task1(){
    double dT = 0.15;
    int M = 1e4;
    int i,j;
    int N = 9; // Number of unit cells in bcc lattice
    int n=(N+1)*(N+1)*(N+1); // number of atoms in sublattice

    double T= 0.0;


    
    double E_0 = 2*n*(E_AA+E_BB+2*E_AB);
    double dE = E_AA+E_BB-2*E_AB;
    
    double F,Fmin;
    double *Pmin = calloc(M,sizeof(double));
    double *U_T = calloc(M,sizeof(double));
    double *C_T = calloc(M-1,sizeof(double));
    double P=0;
    double T_c= 2*dE/kb;

    /* Calculating the minumum P */
    for(i=0;i<M;i++){
        T += dT; Fmin=0; P = 0; Pmin[i] = 0; F = 0;
        for(j=0;j<M-1;j++){
            P+=0.0001;
            F = E_0 - 2*n*P*P*dE - 2*n*kb*T*log(2) + 
                n*kb*T*((1+P)*log(1+P)+(1-P)*log(1-P));
            if (F<Fmin) {
                Fmin = F;
                Pmin[i] = P;
            }
        }
    }
    write_to_file("P.csv",Pmin,dT);

    printf("Tc = %e K\n",T_c);

    for (i=0;i<M;i++){
        U_T[i] = E_0 - 2*n*Pmin[i]*Pmin[i]*dE;
    }
    write_to_file("U.csv",U_T,dT);
    


T = 0.0;
    //averaging over 100 steps to calculate C:
    for (i=0; i<M-100; i++) {
        C_T[i] =(U_T[i+100] - U_T[i])/(100*dT);
        T +=dT;
    }
    write_to_file("C.csv",C_T,dT);

}

void InitLattice(int n, int lattice[], int type){
  int i;
  for( i = 0 ; i < n ; i++ ) {
    lattice[i] = type;
  }
}

void SwapAtoms(int n, int lattice_1[], int lattice_2[], int index1, int index2){
  int swap1, swap2;
  
  swap1 = index1<n ? lattice_1[index1] : lattice_2[index1-n];
  swap2 = index2<n ? lattice_1[index2] : lattice_2[index2-n];
  if( index1 < n ) {
    lattice_1[index1] = swap2;
  } else {
    lattice_2[index1-n] = swap2;
  }
  if( index2 < n ) {
    lattice_1[index2] = swap1;
  } else {
    lattice_2[index2-n] = swap1;
  }
}

