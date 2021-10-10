#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>

#include <random>
#include <iostream>

#include "utils/io_utils.h"
#include "utils/array.h"
#include "utils/math_utils.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

typedef u_int8_t byte;

using namespace array;
using namespace std;

void compute_density_chunk(u_int32_t *__restrict density, float *__restrict points,
                           const size_t n_points, const register size_t N,
                           const size_t N2, const float N_inv,
                           const float R, const float R2,
                           const int my_task, const int n_tasks)
{
    // Compute start and size of each chunk
    const size_t start = (N / n_tasks) * my_task + MIN(N % n_tasks, my_task);
    const size_t Nx = N / n_tasks + (size_t)(N % n_tasks > my_task);
    size_t Nx_range[2]{start, start + Nx};
    size_t default_range[2]{0lu, N};

    // Initialize density to 0
    memset(&density[start * N2], 0, Nx * N2 * sizeof(uint32_t));

    // Possible improvement is avoid checking all points for each task,
    // maybe with a presorting of the points by the x coordinate and/or
    // allocate/migrate the points to the memory bank closest to the task
    // like in the MPI version of the code

    // printf("my_id: %d, start: %ld, end: %ld, step: %ld, n_points: %ld, R: %f\n",
    // omp_get_thread_num(), Nx_range[0], Nx_range[1], Nx, n_points, R);
    for (register size_t i = 0; i < n_points * 3; i += 3)
    {
        update_density_matrix_OMP(density, &points[i], N, N2, N_inv, R, R2, Nx_range, default_range);
    }
}

int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng;
    std::uniform_real_distribution<float> udist(0, 1);

    if (argc < 3)
    {
        printf("Error: you must provide a value for the grid number N and the radius R\n");
        return 1;
    }

    const u_int32_t N = atol(argv[1]);
    const size_t N2 = N * N;
    const size_t N3 = N2 * N;
    const register float R = (float)atof(argv[2]);
    const unsigned int n_dims = 3;

    float *points;
    uint32_t n_points;

    // Read input file
    if (argc == 4)
    {
        // Open input file
        FILE *file = fopen(argv[3], "rb");

        // Read number of points
        n_points = get_number_of_points(file);

        // Allocate memory for the points
        const size_t n_bytes = n_points * n_dims * sizeof(float);
        points = (float *)malloc(n_bytes);

        // Read the 3d points
        const size_t n_bytes_read = read_bytes(file, (byte *)points, n_bytes);

        assert(n_bytes_read == n_bytes &&
               "Error: The file has less points than expected\n");

        fclose(file);
    }

    // Generate random points
    else
    {
        n_points = N3 / 10;
        points = (float *)malloc(n_points * sizeof(float) * n_dims);
        set_random_seed(rng, rd);
        array::rand_array(points, n_points * n_dims, rng, udist);
    }

    // Density matrix computation
    uint32_t *density = (uint32_t *)malloc(N3 * sizeof(uint32_t));
    {
        const register float N_inv = 1.f / (float)N;

        #pragma omp parallel shared(density, points)
        {
            #pragma omp single nowait
            {
                const int n_threads = omp_get_num_threads();
                const int n_tasks = n_threads * 4;
                const float R2 = R * R;
                for (int i = 0; i < n_tasks; i++)
                {
                    // Since density is shared we need to avoid false sharing, so
                    // we schedule the tasks with different priorities
                    #pragma omp task priority(i % n_threads)
                    compute_density_chunk(density, points, n_points,
                                          N, N2, N_inv, R, R2, i, n_tasks);
                }
            }
        }
    }

    // Write density matrix to file
    {
        FILE *file = fopen("density.bin", "wb");
        const size_t dtype_size = sizeof(uint32_t);

        // Write the grid number N
        write_bytes(file, (byte *)&N, dtype_size);

        // Write the density
        write_bytes(file, (byte *)density, N3 * dtype_size);

        fclose(file);
    }

    free(points);
    free(density);

    return 0;
}