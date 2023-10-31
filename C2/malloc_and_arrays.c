#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "linalg.h"

int size;
/*

double scalarProduct(double* x1, double* x2){
    double sum = 0;
	for(int i = 0; i<size; i++){
        sum += x1[i]*x2[i];
	}	
    return sum;
}

double computeDistance(double* x1, double* x2){
    double distance = sqrt( pow(x1[0]-x2[0],2) + 
                            pow(x1[1]-x2[1],2) +
                            pow(x1[2]-x2[2],2));
    return distance;

}
*/

void write_to_file(char *fname,double** matrix){
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "x,y,z\n");
    for (int i = 0; i < size; i++) {
       
        fprintf(fp, "%lf,%lf,%lf\n",matrix[i][0],matrix[i][1],matrix[i][2]);
        
    }
    fclose(fp);
}
int main(int argc, char** argv){
    size = atoi(argv[1]);
    double* array1 = (double*)malloc(sizeof(double)*size);
    double* array2 = (double*)malloc(sizeof(double)*size);

    double* coordinate = (double*)malloc(sizeof(double)*3*size);
    double** coordinates = (double**)malloc(sizeof(double*)*size);
    for (int i = 0, j = 0; i< size; i++, j+=3){
        coordinates[i] =  coordinate + j;
    }
    for(int i = 0; i < size; i++){
        for(int j = 0; j < 3;j++){
            coordinates[i][j] = 3.0+i;
        }
    }
    double distance = computeDistance(coordinates[0],coordinates[size-1]);
    printf("%lf\n",distance);
    write_to_file("matrix.csv",coordinates);
      free(coordinates);
    free(coordinate);
  

    return 0;
}
