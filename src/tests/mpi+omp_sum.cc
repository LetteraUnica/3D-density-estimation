#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>


int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    #pragma omp parallel
    {
        printf("%d, %d\n", my_rank, omp_get_max_threads());
    }

    if (my_rank == 0) {
        size_t N = 10000000000;
        
        double register sum = 0.;
        #pragma omp parallel for reduction (+:sum)
        for(size_t i = 0; i < N; i+=1){
            sum += (double)i;
        }

        double Nd = (double) N;
        printf("%lf %lf\n", sum, (Nd*(Nd-1) / 2));
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    
    return 0;
}