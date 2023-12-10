
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define kb 8.617e-5
#define PI 3.141592653589
  int blocks = 75;
void write_blocks(char *fname,double* x, double* y ) {
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "Block_size,estimate\n");
    
    for (int i = 0; i < blocks; ++i) {
        fprintf(fp, "%f,%f\n",x[i] ,y[i] );
    }
    fclose(fp);
}

void write_to_file(char *filename, double *data, int len);
double autocorrelation(double *array, int nlines);
double** bcc(int N);
void init_bcc(int **nnList[][8], int n, int N);
void swing(int r1,int r2,int *species);
double get_energy(int r1, int r2,int *species, double **neighbour, double band_energies[3]);
double get_P(int *species, int N);
double get_r(int *species, int N, double **neighbour);
double E_initial(int species[], double **neighbour, int n, double band_energies[3]);
int accepted = 0;

double average(double *v1,unsigned int len);
double variance(double *v1, unsigned int len);


double block_average(double *data,
	                 int block_size,
	                 int data_len
	                )
{
    if (block_size == 0){
        block_size++;
    }
    int lenPerBlock = data_len/(block_size);
    double F[lenPerBlock];
    for (int i=0; i<lenPerBlock; i++)
        F[i] = 0;
    
    for (int i=0; i<lenPerBlock; i++){
        for (int k=0; k<block_size; k++){
            F[i] += data[k+i*block_size]/((double)block_size);
        }
    }
    double variance_F = 0;
    double mean_F = 0;
    for (int i=0; i<lenPerBlock; i++){
        variance_F += F[i]*F[i]/((double)lenPerBlock);
        mean_F += F[i]/((double)lenPerBlock);
    }
    variance_F = variance_F - mean_F*mean_F;
    double var = variance(data,data_len);

    double s = (double)block_size * variance_F / var;
    return s;
}



int main()
{
    int nsteps = 1e6; // Number of measure steps
    int Neq = 1e6; // Number of equlibration steps for the first temperature
    int Neql = 4e6; // Number of equilibration steps for all other temperatures
    int i,t,Eswitch;
    int N = 10; // Number of unit cells in one sublattice
    int n=(N)*(N)*(N); // number of atoms in sublattice
    double T_start = 0, T = T_start;
    double Tmax = 1500; // Max temp
    double dT = 20;

    int T_len = (int)(Tmax-T_start)/dT;
    
    // Bond energies
    double E_AA=-436e-3; // in eV
    double E_BB=-113e-3;
    double E_AB=-294e-3;
    double band_energies[3] = {E_AA,E_AB,E_BB};
    
    // Declare matrices
    double *P_avg = calloc(T_len,sizeof(double));
    double *U_avg = calloc(T_len,sizeof(double));
    double *r_avg = calloc(T_len,sizeof(double));
	    double s_avgP[T_len];
    double s_avgU[T_len];
    double s_avgr[T_len]; 
	
    double **P_n = malloc((nsteps)*sizeof(double));
    double **U_n = malloc((nsteps)*sizeof(double));
    double **r_n = malloc((nsteps)*sizeof(double));
    for (i=0; i<nsteps; i++) {
        P_n[i] = malloc(T_len*sizeof(double));
        U_n[i] = malloc(T_len*sizeof(double));
        r_n[i] = malloc(T_len*sizeof(double));
    }


    // Declare energy variables
    double E=0,dE=0, U0=0,Unew=0;

    double **neighbour;
    neighbour = bcc(N-1);

    const gsl_rng_type *q;
    gsl_rng *k;
    gsl_rng_env_setup();
    q = gsl_rng_default;
    k = gsl_rng_alloc(q);
    gsl_rng_set(k,time(NULL));
    int r1,r2;
    
    // Declarations for particletypes
    int *species = malloc((2*n)*sizeof(int));;
    
    // Assign particle type 0 for A and 1 for B
    for (i=0;i<n; i++) 
    {
        species[i] = 0;
        species[i+n] = 1;
    }


    //temperature for-loop
    for(t=0;t<T_len;t++){
        T = T_start + t*dT;
        
        // Calculate initial total energy
        E = E_initial(species, neighbour, n, band_energies);
        // equilibriate the first iteration
        if (t==1) {Neq = Neql;}
        
        // Iterate the system
        for (i=0; i<nsteps+Neq; i++) {
            // Assign two particles to switch places
           
            r1 = (int) gsl_rng_uniform_int(k, 2*n);
            r2 = (int) gsl_rng_uniform_int(k, 2*n);


            U0 = get_energy(r1, r2, species, neighbour, band_energies);//get_energy(neighbour, band_energies, species, n); 
            swing(r1,r2,species);
            Unew = get_energy(r1, r2, species, neighbour, band_energies);//get_energy(r1, r2, species, neighbour, band_energies);

            dE = Unew - U0;
            // Check trial configuration
            if (dE<=0 || exp(-dE/(kb*T))>gsl_rng_uniform(k)) { // 
                E += dE;
                accepted++;
            }
            else {
                swing(r1,r2,species);
            }
            
            if(i>=Neq){       
                P_n[i-Neq][t] = get_P(species,n);
                U_n[i-Neq][t] = (double)E;
                r_n[i-Neq][t] = (double)get_r(species, n, neighbour);
     
            }
        }
        for (int i = 0; i < nsteps; ++i)
        {
            U_avg[t] += U_n[i][t]/((double)nsteps);
            P_avg[t] += P_n[i][t]/((double)nsteps);
            r_avg[t] += r_n[i][t]/((double)nsteps);
        }
		
		
        
       
		s_avgP[t] = autocorrelation((P_n[t]), T_len);
		s_avgU[t] = autocorrelation((U_n[t]), T_len);
		s_avgr[t] = autocorrelation((r_n[t]), T_len);
         
		printf("Progress: %d/%d\n",t,T_len);

}   
   

/*
    double* block_sizes  = (double*)malloc(sizeof(double)*blocks);
    double means[3];
    means[0] = average(U_avg,T_len);
    means[1] = average(P_avg,T_len);
    means[2] = average(r_avg,T_len);
   //shift the mean:

    for (int i = 0; i < T_len; ++i){
        U_avg[i] -= means[0];
        P_avg[i] -= means[1];
        r_avg[i] -= means[2];
    }

    for(int t = 0; t<blocks;t++ ){
        int block_size = t;

        block_sizes[t] = block_size;
        s_avgU[t] =  block_average( U_avg, block_size, T_len);
        s_avgP[t] =  block_average( P_avg, block_size, T_len);
        s_avgr[t] =  block_average( r_avg, block_size, T_len);

    }
    write_blocks("sU.csv",block_sizes,s_avgU);
    write_blocks("sP.csv",block_sizes,s_avgP);
    write_blocks("sr.csv",block_sizes,s_avgr);
    */
    double C_T[T_len-1];// = malloc((T_len-1)*sizeof (double)); // Specific heat capacity of our system of atoms (alloy)
    double s_avgC[T_len -1];

    for (int i = 0; i<T_len-1; i++)
        {
        C_T[i] = (double)(U_avg[i+1]-U_avg[i])/(dT); // From the definition of the derivative

        }  
		
	

    FILE *file;
    file = fopen("T.csv","w");
    double current_temp = T_start;
    for (i=0; i<T_len; i++) {
        current_temp = T_start + i*dT;
        fprintf(file, "%e \n",current_temp);
    }
    fclose(file);
	
    write_to_file("P.csv", P_avg, T_len);
    write_to_file("U.csv", U_avg, T_len);
    write_to_file("r.csv", r_avg, T_len);
    write_to_file("C.csv", C_T, T_len-1);


    write_to_file("sP.csv", s_avgP, T_len);
    write_to_file("sU.csv", s_avgU, T_len);
    write_to_file("sr.csv", s_avgr, T_len);   
	//write_to_file("sC.csv", s_avgC, T_len-1);
    

    free(neighbour);
    free(P_n); free(P_avg);
    free(U_n); free(U_avg);
    free(r_n);  free(r_avg);
    free(species);

    printf("%lf\n",accepted/(double)((nsteps+Neq)*T_len));
    return 0;
}

void write_to_file(char *filename, double *data, int len){
    FILE *file;
    file = fopen(filename,"w");

    for (int i=0; i<len; i++) {
        fprintf(file, "%e\n",data[i]);
    }
    fclose(file);


}

double autocorrelation(double *array, int nlines){
    int k =400;

	double data[nlines];

    for (int i = 0; i < nlines; ++i)
    {
        data[i] = array[i];
    }

        //center around mean

      double main_mean=0;
      
      for (int i=0; i<nlines; i++){
        main_mean += data[i];
      }
      main_mean = main_mean / (double)nlines;
      for (int i=0; i<nlines; i++){
        data[i]=data[i]-main_mean;
      }
     

      double f_ik[k];
      double phi[k];
      int guess=0;
      double s = exp(-2.0);
      for (int i=0; i<k; i++){
        f_ik[i]=0;
        phi[i]=0;
      }
      double fsq = 0;
      double sqf =0;

      FILE * corrfunc;
      for (int i = 0; i<k; i++) {
        for (int j=0; j<nlines-i; j++){
          f_ik[i] += data[j]*data[i+j]/(double)(nlines-i);
        
        }
      }
      // averages
      for (int i=0; i<nlines; i++){
        fsq += data[i]*data[i]/nlines;
         sqf += data[i]/(nlines);
      }
      corrfunc = fopen("phi_ny.csv", "w");
      for (int i=0; i<k; i++){
        phi[i] = (f_ik[i] - sqf*sqf) / (fsq - sqf*sqf);
        fprintf(corrfunc, "%d \t %f\n", i, phi[guess]);
        

        if ( sqrt((s-phi[i])*(s-phi[i])) < sqrt((s-phi[guess])*(s-phi[guess] ))){
          guess = i;
        }
      }
      return guess;
  }
// Construct neighbour list
double** bcc(int N){

    int n=(N+1)*(N+1)*(N+1); // number of atoms in sublattice
    int i,j,k;
    int a;

    double **neighbour = malloc(2*n*sizeof(double));
    for (i=0; i<2*n; i++) {
        neighbour[i] = malloc((8)*sizeof(double));
    }
    
    double fi[8] = {0, 1, -(N+1)*(N+1), -(N+1)*(N+1)+1,  -(N+1)*(N+1)-(N+1), -(N+1)*(N+1)-(N+1)+1, -(N+1), -(N+1)+1}, f[8];
    
    // Create neighbour list 
    for (i=0; i<n; i++) {
        for (j=0; j<8; j++) {
            f[j] = fi[j];
        }
        
        // Check if on one wall 1
        int d=0;
        for (k = 0; k<N+1 ; k++) {
            
            if (d<=i && i<d+(N+1)) {
                for (j=4; j<8; j++) {
                    f[j] += (N+1)*(N+1);
                }
            }
            d += (N+1)*(N+1);
            
            }
        // end check for wall 1
        
        // check if on wall 2 
        if (i<(N+1)*(N+1)) {
            for (j=2; j<6; j++) {
                f[j] += (N+1)*(N+1)*(N+1);
            }
        }
        // end check for wall 2 
        
        for(j=0;j<8;j++){
            a=0;
            // Check if on wall 3 
            a = n + i + f[j];
            if((j%2)==0){a += ((a%(N+1))==0)*(N+1);}
            //7 end check for wall 3 
            neighbour[i][j] = a-1;
        }
    }
    
    for(i=0;i<n;i++){
        for (j=0; j<8; j++) {
            k = (int)neighbour[i][j];
           // printf("%d\n",k);
            neighbour[k][j]=i;
        }
    }
    return neighbour;
}

// Calculate initial energy
double E_initial(int species[], double **neighbour, int n, double band_energies[3]){
    int j,k,a,b;
    double Einitial=0;
    
    for (j = 0; j<n; j++) {
        for (k=0; k<8; k++) {
            a = neighbour[j][k];
            b = species[j] + species[a];
            Einitial += band_energies[b];
        }
    }
    
    return Einitial;
}
double get_energy(int r1, int r2,int *species, double **neighbour, double band_energies[3]){
    int a=0,b=0,k;
    double local_energy=0;
    
    // Calculate current bond energies for the switching particles
    for (k=0; k<8; k++) {
        a = neighbour[r1][k];
        b = species[r1] + species[a];
        local_energy += band_energies[b];
        
        a = neighbour[r2][k];
        b = species[r2] + species[a];
        local_energy += band_energies[b];
    }
    return local_energy;
}


void swing(int r1,int r2,int *species){
    int tmp;
    tmp = species[r1];
    species[r1] = species[r2];
    species[r2] = tmp;
    
}
double get_P(int *species, int N)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
        sum += species[i];

    return fabs(1 - 2 * sum / N);
}

double get_r(int *species, int N, double **neighbour)
{
    int index=0;
    double sum = 0.0;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            index = neighbour[i][j];
            if (species[i] + species[index] == 1)
                sum++;
        }
    }
    sum /= N;
    return fabs(1.0 / 4.0 * (sum - 4));
}

double
average(
        double *v1,
        unsigned int len
       )
{
    double sum = 0.0    ;
    for(int i = 0; i< len; i++){
        sum += v1[i];
    }
    sum = sum/(double)len;
    return sum;
}
double variance(double *v1, unsigned int len)
{
    double mean = average(v1,len);
    double sum = 0;
    for(int i = 0; i<len; i++){
        sum += pow((v1[i]-mean),2);
    }
    
    sum = sqrt(sum/(double)len);
    return sum;
}
