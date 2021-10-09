#include <stdio.h>
#include <mpi.h>
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

    u_int32_t N = atol(argv[1]);
    size_t N2 = N * N;
    size_t N3 = N2 * N;
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
        uint32_t max_points = 1024;
        ResizableArray *DS = create_empty_data_structure(n_processors, max_points);

        // Read input file
        if (argc == 4)
        {
            // Open input file
            FILE *file = fopen(argv[3], "rb");

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

                generated_points += batch_size;

                // Insert points in the DS
                insert_points(points, current_batch_size, DS, R, n_processors);
            }

            free(points);
        }

        // Send local_n_points to each processor
        for (int i = 1; i < n_processors; i++)
        {
            MPI_Send(&DS[i].cur_points, 1, MPI_UINT32_T, i, 0, MPI_COMM_WORLD);
        }

        // Send the points to each processor
        for (int i = 1; i < n_processors; i++)
        {
            MPI_Send(DS[i].data, DS[i].cur_points * n_dims, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
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
        MPI_Recv(&local_n_points, 1, MPI_UINT32_T, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        local_points = (float *)malloc(local_n_points * sizeof(float) * n_dims);

        MPI_Recv(local_points, local_n_points * n_dims, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    // Density matrix computation
    size_t start = (N / n_processors) * my_rank + MIN(N % n_processors, my_rank);
    size_t Nx_range[2] = {start, start + Nx};
    float N_inv = 1.f / (float)N;
    for (size_t i = 0; i < local_n_points; i += 1)
    {
        fast_update_density_matrix(local_density, &local_points[3 * i], N, N_inv, R, Nx_range);
    }

    // Write density matrix to file
    MPI_File file;
    int access_mode = MPI_MODE_CREATE    /* Create the file if it does not exist */
                      | MPI_MODE_WRONLY; /* With write access */
    MPI_File_open(MPI_COMM_WORLD, "density.bin", access_mode, MPI_INFO_NULL, &file);

    // Write the grid number N
    MPI_Status status;
    if (my_rank == 0)
    {
        MPI_File_write(file, &N, 1, MPI_UINT32_T, &status);
    }

    // Write the density matrix
    size_t offset = (start * N2 + 1) * sizeof(uint32_t);
    size_t count = Nx * N2;
    MPI_File_write_at_all(file, offset, local_density, count, MPI_FLOAT, &status);

    MPI_File_close(&file);

    free(local_density);
    free(local_points);

    MPI_Finalize();

    return 0;
}