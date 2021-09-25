#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <io_utils.h>
#include <array.h>

#include <random>

#include <stdint.h>
#include <limits.h>


typedef u_int8_t byte;

using namespace array;

std::random_device rd;
std::mt19937 rng;
std::uniform_real_distribution<float> udist(0, 1);

void set_random_seed()
{
    rng.seed(rd());
}

void create_input_file(const char *filename, int n_points)
{
    set_random_seed();

    FILE *ptr = fopen("test.bin", "rb");
    write_bytes(ptr, (byte *)&n_points, 4);

    for (int i = 0; i < n_points; i++)
    {
        float *random_points = array::rand_array(3);
        byte *buffer = float_to_byte(random_points, 3);
        write_bytes(ptr, buffer, 3 * sizeof(float));
    }
}

void insert_points(float *points, size_t n_points, ResizableArray *DS, float R, int n_processors)
{
    for (size_t i = 0; i < n_points; i++)
    {
        insert_point(DS, points[3 * i], points[3 * i + 1], points[3 * i + 2], R, n_processors);
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int my_rank, n_processors;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_processors);

    if (my_rank == 0 && argc < 4) {
        printf("Error: you must provide a value for the grid number N and the radius R");
        return 1;
    }

    unsigned long N = atol(argv[1]);
    float R = (float)atof(argv[2]);
    unsigned int n_dims = 3;

    if (my_rank == 0)
    {
        // Get number of points
        int n_points = N/10;
        FILE *file;
        if (argc == 4)
        {
            // Open input file
            file = fopen(argv[3], "rb");

            // Read number of points
            byte *buffer = read_input_file(file, 4);
            memcpy(&n_points, buffer, 4);
        }

        size_t max_points = n_points * (1./(float)n_processors + 2*(float)R) + 1024;
        ResizableArray *DS = create_empty_data_structure(n_processors, max_points);

        // Read input file
        if (argc == 4)
        {
            // Read the 3d points
            size_t points_per_read = 1024;
            size_t conversion_factor = sizeof(float) * 3;
            size_t block_size = points_per_read * conversion_factor;
            size_t actual_n_points = 0;

            for (size_t i = 0; i < n_points; i += points_per_read)
            {
                byte *buffer = create_empty_buffer(block_size);

                size_t n_bytes_read = read_bytes(file, buffer, block_size);
                size_t n_points_read = n_bytes_read / conversion_factor;

                assert(n_bytes_read == n_points_read * conversion_factor &&
                       "Warning: The coordinates of the last point are missing,"
                       "discarding them\n");

                actual_n_points += n_points_read;

                // Convert to float
                float *points = byte_to_float(buffer, n_points_read);

                // Insert points in the DS
                insert_points(points, n_points_read, DS, R, n_processors);
            }

            if (actual_n_points != n_points)
            {
                printf("Warning: The number of points found in the file are"
                       "different than what is specified in the file header,"
                       "specified: %ld, found: %ld, the program will continue with"
                       "the number of found points\n",
                       n_points, actual_n_points);

                n_points = actual_n_points;
            }

            fclose(file);
        }

        // Generate random points
        else
        {
            set_random_seed();

            size_t batch_size = 1024;
            size_t conversion_factor = 3;

            for (size_t i = 0; i < n_points; i += batch_size)
            {
                float *points = array::rand_array(batch_size * conversion_factor);

                // Insert points in the DS
                insert_points(points, batch_size, DS, R, n_processors);
            }
        }


        for (int i=1; i<n_processors; i++) {
            MPI_Send(&DS[i].cur_points, 1, MPI_UNSIGNED_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD);
        }
        for (int i=1; i<n_processors; i++) {
            MPI_Send(DS[i].data, DS[i].cur_points*3, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
        }

        unsigned int *local_density = (unsigned int*)calloc(N*N*N / n_processors, sizeof(unsigned int));

        unsigned long int n_local_points = DS[0].cur_points;

    }

    if (my_rank != 0) {
        unsigned int *local_density = (unsigned int*)calloc(N*N*N / n_processors, sizeof(unsigned int));
        
        unsigned long int n_local_points;
        MPI_Status status;
        MPI_Recv(&n_local_points, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        float* local_points = (float*)malloc(n_local_points * sizeof(float) * 3);

        MPI_Recv(local_points, n_local_points*3, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }
}