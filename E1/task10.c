#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <gsl/gsl_const_mksa.h>

int N = 32;
double PI = 3.141592;

double** transformation_matrix;
double* P;
double* Q;
double** E;
double* Etot;
double* Etot_dt;
double* omega;
double** create_2D_array(unsigned int column_size,unsigned int row_size);

void calc_acc(double *a, double *u, double m, double kappa,double alpha)
{
    a[0] = kappa*(u[1] - 2*u[0])/(m) + alpha*(u[1]*u[1]-2*u[0]*u[1])/m;
    a[N-1] = kappa*(- 2*u[N-1] +
                 u[N-2])/m + 
                 alpha*(2*u[N-1]*u[N-2]-
                 u[N-2]*u[N-2])/m;

    /* Calculating the acceleration of the inner points */
    for (int i = 1; i < N-1; i++){
        a[i] = kappa*(u[i + 1] - 2*u[i] + u[i - 1])/(m) + alpha*(u[i + 1]*u[i + 1] - u[i - 1]*u[i - 1] + 2*u[i]*(u[i - 1]-u[i + 1]))/m;
    }

    
}
void create_transformation_matrix(double** transformation_matrix){
    double normalize = 1.0/((double)N +1.0);
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            transformation_matrix[i][j] = sqrt(2*normalize)*sin((j+1)*(i+1)*PI*normalize);
        }
    }
}
void transform(double** transformation_matrix,double* positions, double* Q){
    
    for(int i = 0; i<N;i++){
        double sum = 0.0;
        for(int j = 0; j<N;j++){
            sum += positions[j]*transformation_matrix[i][j];
        }
        Q[i] = sum;
    }
}


void velocity_verlet_multiple_steps(int n_timesteps, double *velocities,
                                     double *q, double dt, double m, double kappa, double alpha){

    double accelerations[N];

    FILE *fp = fopen("T10.csv", "w");
    fprintf(fp, "time,1,2,3,4,5\n");

    double T;
    double temp[N];
    calc_acc(accelerations, q, m, kappa,alpha);
    int counter = 0;
    for (int i = 0; i < n_timesteps; i++) {
    
        for (int j = 0; j < N; j++) {
            velocities[j] += dt * 0.5 * accelerations[j];
        }
        
        for (int j = 0; j < N; j++) {
            q[j] += dt * velocities[j];
        }
        
        calc_acc(accelerations, q, m, kappa,alpha);
        for (int j = 0; j < N; j++) {
            velocities[j] += dt * 0.5 * accelerations[j];
        }
        

        
        if (i%1000 == 0 || i ==1){
            transform(transformation_matrix,q,Q);
            transform(transformation_matrix,velocities,P);
            //Compute time average:
            for(int ix = 0; ix < N; ix++){
                E[ix][i] = 0.5 * (pow(P[ix], 2) + pow(omega[ix]*Q[ix], 2) );
                if(counter != 0){
                    T = counter*dt*1000;
                    temp[ix] += E[ix][i]*dt;
                    Etot[ix] = temp[ix]/T;
                }
                
            }
            
                    fprintf(fp, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                i * dt, Etot[0], Etot[1], Etot[2], Etot[3], Etot[4], Etot[5], Etot[6], Etot[7], Etot[8],
                Etot[9], Etot[10], Etot[11], Etot[12], Etot[13], Etot[14], Etot[15], Etot[16], Etot[17], Etot[18],
                Etot[19], Etot[20], Etot[21], Etot[22], Etot[23], Etot[24], Etot[25], Etot[26], Etot[27], Etot[28],
                Etot[29], Etot[30], Etot[31]);
            counter++;
        }
        
    }
}




int main(){
    int n_timesteps = 2e7;

    double kappa = 1; double dt = 0.15; double alpha = 0.1;double m = 1.0;


    
    double* q = calloc(sizeof(double),N);
    
    double* v = calloc(sizeof(double), N);
    P = calloc(sizeof(double),N);
    P[0] = sqrt(2.0*N);
    Q = calloc(N,sizeof(double));
    omega = calloc(sizeof(double),N);
    Etot = calloc(sizeof(double),N);

    E = create_2D_array(N,n_timesteps);
    

    double factor = 1/((double)(N+1));
    //Set initial conditions:
    
    for(int i = 0; i < N; i++){
        q[i] = 0.0;
        v[i] = sqrt(2*factor) * (P[0]/sqrt(m)) * sin((i+1)*PI*factor);
    }

    
    

            
   
    for(int ix = 0; ix < N; ix++)
    {   
        omega[ix] = 2.0 * sin( (ix + 1) * PI *0.5 * factor ) ;
    }


    transformation_matrix = create_2D_array(N,N);

    create_transformation_matrix(transformation_matrix);

    transform(transformation_matrix,q,Q);
    transform(transformation_matrix,v,P);

    velocity_verlet_multiple_steps(n_timesteps, v,q, dt, m, kappa,alpha);

    return 0;
}


double **
create_2D_array(
                unsigned int column_size,
                unsigned int row_size
               )
{
    double* arrayEntries = (double*)malloc(sizeof(double)*row_size*column_size);
    double** array       = (double**)malloc(column_size*sizeof(double*));
    for (int i = 0,j=0; i < column_size; ++i,j+=row_size) {
        array[i] = arrayEntries +j;
    }
    return array;
}