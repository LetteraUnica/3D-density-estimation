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

    u_int32_t N = atol(argv[1]);
    size_t N2 = N * N;
    size_t N3 = N2 * N;
    float R = (float)atof(argv[2]);
    unsigned int n_dims = 3;

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
        size_t n_bytes = n_points * n_dims * sizeof(float);
        points = (float *)malloc(n_bytes);

        // Read the 3d points
        size_t n_bytes_read = read_bytes(file, (byte *)points, n_bytes);

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
    size_t Nx_range[2] = {0lu, N};
    uint32_t *density = (uint32_t *)calloc(N3, sizeof(uint32_t));
    float N_inv = 1.f / (float)N;
    for (size_t i = 0; i < n_points; i += 1)
    {
        fast_update_density_matrix(density, &points[3 * i], N, N_inv, R, Nx_range);
    }

    // Write density matrix to file
    {
        FILE *file = fopen("density.bin", "wb");
        size_t dtype_size = sizeof(uint32_t);

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