#include <stdio.h>
#include <math.h>
#include <time.h>


double scalarProduct(double* x1, double* x2,int size){
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
