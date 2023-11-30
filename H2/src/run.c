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
void task2(){


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
    double dT = 0.1;
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


void InitNeighbors(int n, int GA[][8], int GB[][8])
{
  int i;
  // NB! GA contains indices to LatticeB and GB indexes in latticeA
  // ie. the neighbours to all particles in
  // sublattice A are located in sublattice B.
  for(i=0; i<n; i++){
    // All inner points
    GA[i][0] = i;         GB[i][0] = i+111;
    GA[i][1] = i-1;       GB[i][1] = i+110;
    GA[i][2] = i-10;      GB[i][2] = i+101;
    GA[i][3] = i-11;      GB[i][3] = i+100;
    GA[i][4] = i-100;     GB[i][4] = i+11;
    GA[i][5] = i-101;     GB[i][5] = i+10;
    GA[i][6] = i-110;     GB[i][6] = i+1;
    GA[i][7] = i-111;     GB[i][7] = i;

    if(i<100){  // The bottom side
      GA[i][0] = i;          
      GA[i][1] = i-1;        
      GA[i][2] = i-10;        
      GA[i][3] = i-11;       
      GA[i][4] = i+n-100;     
      GA[i][5] = i+n-101;     
      GA[i][6] = i+n-110;    
      GA[i][7] = i+n-111;    
    }
        
    if(i%100<10){ // Front
      GA[i][0] = i;           //GB[i][0] = i+111;
      GA[i][1] = i-1;         //GB[i][1] = i+110;
      GA[i][2] = i+100-10;    //GB[i][2] = i+101;
      GA[i][3] = i+100-11;    //GB[i][3] = i+100;
      GA[i][4] = i-100;       //GB[i][4] = i+11;
      GA[i][5] = i-101;       //GB[i][5] = i+10;
      GA[i][6] = i-10;        //GB[i][6] = i+1;
      GA[i][7] = i-11;        //GB[i][7] = i;
    }

    // Next up, all the faces
    if((i%100) >= 90){ // Back
      /*GA[i][0] = i;*/         GB[i][0] = i-100+111;
      /*GA[i][1] = i-1;*/       GB[i][1] = i-100+110;
      /*GA[i][2] = i-10;*/      GB[i][2] = i+101;
      /*GA[i][3] = i-11;*/      GB[i][3] = i+100;
      /*GA[i][4] = i-100;*/     GB[i][4] = i-100+11;
      /*GA[i][5] = i-101;*/     GB[i][5] = i-100+10;
      /*GA[i][6] = i-110;*/     GB[i][6] = i+1;
      /*GA[i][7] = i-111;*/     GB[i][7] = i;
    }
    if((i+1)%10==0){ // Right
      /*GA[i][0] = i;*/         GB[i][0] = i+101;
      /*GA[i][1] = i-1;*/       GB[i][1] = i+110;
      /*GA[i][2] = i-10;*/      GB[i][2] = i+91;
      /*GA[i][3] = i-11;*/      GB[i][3] = i+100;
      /*GA[i][4] = i-100;*/     GB[i][4] = i+1;
      /*GA[i][5] = i-101;*/     GB[i][5] = i+10;
      /*GA[i][6] = i-110;*/     GB[i][6] = i-9;
      /*GA[i][7] = i-111;*/     GB[i][7] = i;
    }
    if(i%10==0){ // Left
      GA[i][0] = i;         //GB[i][0] = i+111;
      GA[i][1] = i+9;       //GB[i][1] = i+110;
      GA[i][2] = i-10;      //GB[i][2] = i+101;
      GA[i][3] = i-1;       //GB[i][3] = i+100;
      GA[i][4] = i-100;     //GB[i][4] = i+11;
      GA[i][5] = i-101+10;  //GB[i][5] = i+10;
      GA[i][6] = i-110;     //GB[i][6] = i+1;
      GA[i][7] = i-111+10;  //GB[i][7] = i;
    }
    if(i>(n-100)){ // Top
      /*GA[i][0] = i;*/         GB[i][0] = i+111-100;
      /*GA[i][1] = i-1;*/       GB[i][1] = i+110-100;
      /*GA[i][2] = i-10;*/      GB[i][2] = i+101-100;
      /*GA[i][3] = i-11;*/      GB[i][3] = i+100-100;
      /*GA[i][4] = i-100;*/     GB[i][4] = i+11;
      /*GA[i][5] = i-101;*/     GB[i][5] = i+10;
      /*GA[i][6] = i-110;*/     GB[i][6] = i+1;
      /*GA[i][7] = i-111;*/     GB[i][7] = i;
    }

    // Next up, edges between faces
    if(i<10){ //Främre undre raden 
      GA[i][0] = i;         //GB[i][0] = i+111;
      GA[i][1] = i-1;       //GB[i][1] = i+110;
      GA[i][2] = i+100-10;  //GB[i][2] = i+101;
      GA[i][3] = i+100-11;  //GB[i][3] = i+100;
      GA[i][4] = i+n-100;   //GB[i][4] = i+11;
      GA[i][5] = i+n-101;   //GB[i][5] = i+10;
      GA[i][6] = i+n-110+100;//GB[i][6] = i+1;
      GA[i][7] = i+n-111+100;//GB[i][7] = i;
    }
    if(i%100 == 0){
      GA[i][0] = i;         GB[i][0] = i+111;
      GA[i][1] = i+10-1;    GB[i][1] = i+110;
      GA[i][2] = i+100-10;  GB[i][2] = i+101;
      GA[i][3] = i-11+110;   GB[i][3] = i+100;
      GA[i][4] = i-100;     GB[i][4] = i+11;
      GA[i][5] = i+10-101;  GB[i][5] = i+10;
      GA[i][6] = i+100-110; GB[i][6] = i+1;
      GA[i][7] = i-111+110;  GB[i][7] = i; //?
    }
    if(i%100 == 99){ // Bakre högra raden
    GB[i][0] = i+111-110;
      /*GA[i][1] = i-1;*/       GB[i][1] = i+110-100;
      /*GA[i][2] = i-10;*/      GB[i][2] = i+101-10;
      /*GA[i][3] = i-11;*/      GB[i][3] = i+100;
      /*GA[i][4] = i-100;*/     GB[i][4] = i+11-110;
      /*GA[i][5] = i-101;*/     GB[i][5] = i+10-100;
      /*GA[i][6] = i-110;*/     GB[i][6] = i+1-10;
      /*GA[i][7] = i-111;*/     GB[i][7] = i;
    }
    if(i>=n-10){ //bakre övre raden
      /*GA[i][0] = i;*/         GB[i][0] = i-n-100+111;
      /*GA[i][1] = i-1;*/       GB[i][1] = i-n-100+110;
      /*GA[i][2] = i-10;*/      GB[i][2] = i-n+101;
      /*GA[i][3] = i-11;*/      GB[i][3] = i-n+100;
      /*GA[i][4] = i-100;*/     GB[i][4] = i+11-100;
      /*GA[i][5] = i-101;*/     GB[i][5] = i+10-100;
      /*GA[i][6] = i-110;*/     GB[i][6] = i+1;
      /*GA[i][7] = i-111;*/     GB[i][7] = i;
    }
    if(i>(n-100) && i%10==9){ // högra sidan övre raden
      /*GA[i][0] = i;*/         GB[i][0] = i-n-10+111;
      /*GA[i][1] = i-1;*/       GB[i][1] = i-n+110;
      /*GA[i][2] = i-10;*/      GB[i][2] = i-n-10+101;
      /*GA[i][3] = i-11;*/      GB[i][3] = i-n+100;
      /*GA[i][4] = i-100;*/     GB[i][4] = i-10+11;
      /*GA[i][5] = i-101;*/     GB[i][5] = i+10;
      /*GA[i][6] = i-110;*/     GB[i][6] = i-10+1;
      /*GA[i][7] = i-111;*/     GB[i][7] = i;
    }
    if(i<100 && i%10 == 0){ // vänstra sidan botten-raden
      GA[i][0] = i;         //GB[i][0] = i+111;
      GA[i][1] = i+10-1;    //GB[i][1] = i+110;
      GA[i][2] = i-10;      //GB[i][2] = i+101;
      GA[i][3] = i+10-11;   //GB[i][3] = i+100;
      GA[i][4] = i+n-100;   //GB[i][4] = i+11;
      GA[i][5] = i+n+10-101;//GB[i][5] = i+10;
      GA[i][6] = i+n-110;   //GB[i][6] = i+1;
      GA[i][7] = i+n+10-111;//GB[i][7] = i;
    }

    if(i == 99){ //Hörn - bakre, lägre högra
     GB[i][0] = 100;
     GB[i][1] = 109;
     GB[i][2] = 190;
     GB[i][3] = 199;
     GB[i][4] = 0;
     GB[i][5] = 9;
     GB[i][6] = 90;
     GB[i][7] = 99;
    }
                      
    if(i==0){      
     for (int j = 0; j<7;j++)
      GA[i][0] = 0;           
      GA[i][1] = -1+10;       
      GA[i][2] = -10+100;     
      GA[i][3] = -1+100;     
      GA[i][4] = -100+n;    
      GA[i][5] = -101+n+10;  
      GA[i][6] = n-10;    
      GA[i][7] = n-1;    
    }
    if(i == 999){
      GB[i][0] = 999-n-100-10+111; 
      GB[i][1] = 999-n-100+110;
      GB[i][2] = 999-n-10+101;
      GB[i][3] = 999-n+100;
      GB[i][4] = 999-100-10+11;
      GB[i][5] = 999-100+10;
      GB[i][6] = 999-10+1;
      GB[i][7] = 999;
    }
  }

}