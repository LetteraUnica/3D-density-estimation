#include <stdio.h>
#include <mpi.h>
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

#define BATCH_SIZE 1024

typedef u_int8_t byte;

using namespace array;
using namespace std;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    std::random_device rd;
    std::mt19937 rng;
    std::uniform_real_distribution<float> udist(0, 1);

    int my_rank, n_processors;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processors);

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

    float *local_points;
    uint32_t *local_density;
    uint32_t local_n_points;

    size_t Nx = N / n_processors + (size_t)(N % n_processors > my_rank);

    if (my_rank == 0)
    {
        uint32_t n_points;
        // Create the data structure
        size_t max_points = 1024;
        ResizableArray *DS = create_empty_data_structure(n_processors, max_points);

        // Read input file
        if (argc == 4)
        {
            // Open input file
            FILE* file = fopen(argv[3], "rb");

            // Read number of points
            n_points = get_number_of_points(file);

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
                insert_points((float *)buffer, points_to_read, DS, R, n_processors);
            }

            free(buffer);
            fclose(file);
        }

        // Generate random points
        else
        {
            set_random_seed(rng, rd);
            n_points = N3 / 10;

            size_t batch_size = 1024;
            size_t conversion_factor = n_dims * sizeof(float);
            size_t generated_points = 0;

            float *points = (float *)malloc(batch_size * conversion_factor);
            while (generated_points < n_points)
            {
                size_t current_batch_size = MIN(batch_size, n_points - generated_points);
                array::rand_array(points, current_batch_size * n_dims, rng, udist);

                // Insert points in the DS
                insert_points(points, current_batch_size, DS, R, n_processors);
                generated_points += batch_size;
            }

            free(points);
        }

        // Send local_n_points to each processor
        for (int i = 1; i < n_processors; i++)
        {
            MPI_Send(&DS[i].cur_points, 1, MY_MPI_SIZE_T, i, 0, MPI_COMM_WORLD);
        }

        // Send the points to each processor
        for (int i = 1; i < n_processors; i++)
        {
            MPI_Send(DS[i].data, DS[i].cur_points * 3, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
            free(DS[i].data);
        }

        // Allocate the density matrix on the master processor
        local_density = (uint32_t *)calloc(N2 * Nx, sizeof(uint32_t));
        local_n_points = DS[0].cur_points;
        local_points = DS[0].data;

        free(DS);
    }

    if (my_rank != 0)
    {
        // Allocate the density matrix on the slaves, in the meantime the master will
        // be reading the input file
        local_density = (uint32_t *)calloc(N2 * Nx, sizeof(uint32_t));

        // Receive points
        MPI_Status status;
        printf("Line 166 sets the number of processors to 0\n");
        MPI_Recv(&local_n_points, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        MPI_Comm_size(MPI_COMM_WORLD, &n_processors);

        local_points = (float *)malloc(local_n_points * sizeof(float) * 3);

        MPI_Recv(local_points, local_n_points * 3, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    // Density matrix computation
    size_t start = (N / n_processors) * my_rank + MIN(N % n_processors, my_rank);
    size_t Nx_range[2] = {start, start + Nx};

    printf("%ld, %ld, %ld\n\n", Nx_range[0], Nx_range[1], Nx);
    for (size_t i = 0; i < local_n_points; i += 1)
    {
        update_density_matrix(local_density, &local_points[3 * i], N, R, Nx_range);
    }

    // Write density matrix to file
    FILE *file = fopen("density.bin", "wb");
    size_t conversion_factor = sizeof(uint32_t);
    fseek(file, start * N2 * conversion_factor, SEEK_SET);

    size_t cells_per_write = 1024;
    size_t total_cells = Nx * N2;
    for (size_t cells_written = 0; cells_written < total_cells; cells_written += cells_per_write)
    {
        size_t cells_to_write = MIN(cells_per_write, total_cells - cells_written);
        write_bytes(file, (byte *)&local_density[cells_written],
                    cells_to_write * conversion_factor);
    }

    free(local_points);
    free(local_density);
    fclose(file);

    MPI_Finalize();
    return 0;
}