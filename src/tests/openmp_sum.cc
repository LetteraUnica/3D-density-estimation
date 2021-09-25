#include <stdio.h>
// #include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>


void fill_array(float *array, float value, size_t n_elements)
{
    for(size_t i = 0; i < n_elements; i++){
        array[i] = value;
    }
}

int main(int argc, char **argv) {

    size_t N = 10000000000;

    printf("%ld\n", sizeof(unsigned long long));
    
    double register sum = 0.;
    #pragma omp parallel for reduction (+:sum) schedule(dynamic, N/omp_get_num_threads()/1000)
    for(size_t i = 0; i < N; i++){
        sum += i;
    }

    double Nd = (double) N;
    printf("%lf %lf\n", sum, (Nd*(Nd-1) / 2));

    return 0;
}