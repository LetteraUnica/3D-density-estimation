#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>

#include <random>
#include <iostream>

#include "io_utils.h"
#include "array.h"
#include "math_utils.h"

#if SIZE_MAX == UCHAR_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define MY_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define my_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "what is happening here?"
#endif

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

    size_t N = atol(argv[1]);
    size_t N2 = N * N;
    size_t N3 = N * N * N;
    float R = (float)atof(argv[2]);
    unsigned int n_dims = 3;

    uint32_t *density = (uint32_t *)calloc(N3, sizeof(uint32_t));
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
        points = (float *)malloc(n_points * sizeof(float) * n_dims);

        // Read the 3d points
        size_t points_per_read = 1024;
        size_t conversion_factor = n_dims * sizeof(float);
        size_t read_points = 0;

        byte *buffer = create_empty_buffer(points_per_read * conversion_factor);

        while (read_points < n_points)
        {
            size_t points_to_read = MIN(points_per_read, n_points - read_points);
            size_t n_bytes_read = read_bytes(file, buffer,
                                             points_to_read * conversion_factor);

            assert(points_to_read == n_bytes_read / conversion_factor &&
                   "Error: The file has less points than expected\n");

            read_points += points_to_read;

            // Insert points in the DS
            memcpy(&points[read_points * n_dims], buffer,
                   points_to_read * conversion_factor);
        }

        free(buffer);
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
    for (size_t i = 0; i < n_points; i += 1)
    {
        update_density_matrix(density, &points[3 * i], N, R, Nx_range);
    }

    // Write density matrix to file
    {
        FILE *file = fopen("density.bin", "wb");
        size_t conversion_factor = sizeof(uint32_t);
        fseek(file, Nx_range[0] * conversion_factor * N2, SEEK_SET);

        size_t cells_per_write = 1024;
        size_t total_cells = N3;
        for (size_t cells_written = 0; cells_written < total_cells; cells_written += cells_per_write)
        {
            size_t cells_to_write = MIN(cells_per_write, total_cells - cells_written);
            write_bytes(file, (byte *)&density[cells_written],
                        cells_to_write * conversion_factor);
        }

        fclose(file);
    }

    free(points);
    free(density);

    return 0;
}