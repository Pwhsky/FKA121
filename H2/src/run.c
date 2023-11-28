#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <tools.h>
#define kb 8.6173303e-5
#define PI 3.141592653589
void write_to_file(char *fname,double* array){
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "temp, val\n");
    for (int i = 0; i < 1e4-1; i++) {
       
        fprintf(fp, "%lf, %lf\n",i*0.1,array[i]);
    }
    fclose(fp);
}

void task1(){
    int i,j;
    int N = 9; // Number of unit cells in bcc lattice
    int n=(N+1)*(N+1)*(N+1); // number of atoms in sublattice
    //int n_particles = 2*n;
    double T_start = 0, T = T_start;
    double dT = 0.1;
    double E_AA=-436e-3; // in eV
    double E_BB=-113e-3;
    double E_AB=-294e-3;
    
    double E_0 = 2*n*(E_AA+E_BB+2*E_AB);
    double dE = E_AA+E_BB-2*E_AB;
    int M = 1e4;
    double F,Fmin;
    double *Pmin = calloc(M-1,sizeof(double));
    double *U_T = calloc(M-1,sizeof(double));
    double *C_T = calloc(M-1,sizeof(double));
    double P=0; //from 0 to 1 but not 0 or 1
    double T_c= 2*dE/kb;

    /* Calculating the minumum P */
    for(i=0;i<M;i++){
        T += dT; Fmin=0; P = 0; Pmin[i] = 0; F = 0;
        for(j=0;j<M-1;j++){
            P+=0.0001;
            F = E_0 - 2*n*P*P*dE - 2*n*kb*T*log(2) + n*kb*T*((1+P)*log(1+P)+(1-P)*log(1-P));
            if (F<Fmin) {
                Fmin = F;
                Pmin[i] = P;
            }
        }
    }
    write_to_file("P.csv",Pmin);
    /* Calculate energy U */
    //T_c = 2*dE/kb;
    printf("Tc = %e K\n",T_c);
    for (i=0;i<M;i++){
        U_T[i] = E_0 - 2*n*Pmin[i]*Pmin[i]*dE;
    }
    write_to_file("U.csv",U_T);
    
//double *meanUsq = calloc(M-1,sizeof(double));
//double *meansqU = calloc(M-1,sizeof(double));
double meanUsq =0;
double meansqU = 0;
for (int i = 0; i < M-1; ++i)
{
   meanUsq += U_T[i]; // To calculate variance
   meansqU += U_T[i]*U_T[i];

}
meanUsq /= M;
meansqU /= M;
T = T_start;

    //Using finite difference method to calculate C:
    for (i=0; i<M-1; i++) {
        // 1/(kT^2)(<U^2> - <U>^2) 
       // C_T[i] =(U_T[i+1] - U_T[i])/dT;
        C_T[i]= (meansqU -meanUsq*meanUsq)/(kb*T*T);//
        T +=dT;
    }
    write_to_file("C.csv",C_T);


    
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
