#include <stdlib.h>
#include <stdio.h>

int
main()
{
    float *heap_array = (float *)malloc(sizeof(float) * 10000);
    if(heap_array == NULL){
	perror("malloc failed");
	exit(1);
    }
    heap_array[99] = 99;
    printf("%f\n", heap_array[99]);
    free(heap_array);
    return 0;
}
