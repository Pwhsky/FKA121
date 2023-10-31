#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
void write_to_file(char *fname,double* array){
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "x\n");
    for (int i = 0; i < 300; i++) {
       
        fprintf(fp, "%lf\n",array[i]);
    }
    fclose(fp);
}

int main(){
    const gsl_rng_type *T;
    gsl_rng *B;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    B = gsl_rng_alloc(T);
    gsl_rng_set(B, time(NULL));

    double* array = (double*)malloc(sizeof(double)*300);
    
    for(int i = 0; i< 300; i++){
        double rand = gsl_rng_uniform(B);
        array[i] = rand;
    }
    write_to_file("distribution.csv",array);


    return 0;
}
