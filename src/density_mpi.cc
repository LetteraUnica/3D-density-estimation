#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>

#include <random>

#include <io_utils.h>
#include <array.h>

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

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

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
    write_bytes(ptr, (byte *)&n_points, sizeof(int));

    float *random_points = (float *)malloc(3 * sizeof(float));
    for (int i = 0; i < n_points; i++)
    {
        array::rand_array(random_points, 3);
        byte *buffer = float_to_byte(random_points, 3);
        write_bytes(ptr, buffer, 3 * sizeof(float));
    }

    free(random_points);
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

    if (my_rank == 0 && argc < 4)
    {
        printf("Error: you must provide a value for the grid number N and the radius R");
        return 1;
    }

    size_t N = atol(argv[1]);
    float R = (float)atof(argv[2]);
    unsigned int n_dims = 3;
    size_t local_n_points;
    float *local_points;
    unsigned int *local_density;

    size_t Nx = N / n_processors + N % n_processors < my_rank;

    if (my_rank == 0)
    {
        // Get number of points
        int n_points = N / 10;
        FILE *file;
        if (argc == 4)
        {
            // Open input file
            file = fopen(argv[3], "rb");

            // Read number of points
            byte *buffer = read_input_file(file, 4);
            memcpy(&n_points, buffer, 4);
        }

        size_t max_points = n_points * (1. / (float)n_processors + 2 * (float)R) + 1024;
        ResizableArray *DS = create_empty_data_structure(n_processors, max_points);

        // Read input file
        if (argc == 4)
        {
            // Read the 3d points
            size_t points_per_read = 1024;
            size_t conversion_factor = sizeof(float) * n_dims;
            size_t block_size = points_per_read * conversion_factor;
            size_t actual_n_points = 0;

            byte *buffer = create_empty_buffer(block_size);

            for (size_t i = 0; i < n_points; i += points_per_read)
            {
                size_t n_bytes_read = read_bytes(file, buffer, block_size);
                size_t n_points_read = n_bytes_read / conversion_factor;

                assert(n_bytes_read == n_points_read * conversion_factor &&
                       "Warning: The coordinates of the last point are missing,"
                       "discarding them\n");

                actual_n_points += n_points_read;

                // Convert to float
                float *points = byte_to_float(buffer, n_bytes_read);

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

            free(buffer);
            fclose(file);
        }

        // Generate random points
        else
        {
            set_random_seed();

            size_t batch_size = 1024;
            size_t conversion_factor = n_dims * sizeof(float);

            float *points = (float *)malloc(batch_size * conversion_factor);
            for (size_t i = 0; i < n_points; i += batch_size)
            {
                array::rand_array(points, batch_size * n_dims);

                // Insert points in the DS
                insert_points(points, batch_size, DS, R, n_processors);
            }

            free(points);
        }

        for (int i = 1; i < n_processors; i++)
        {
            MPI_Send(&DS[i].cur_points, 1, MPI_UNSIGNED_LONG, i, MPI_ANY_TAG, MPI_COMM_WORLD);
        }

        // Send the points to each processor
        for (int i = 1; i < n_processors; i++)
        {
            MPI_Send(DS[i].data, DS[i].cur_points * 3, MPI_FLOAT, i, MPI_ANY_TAG, MPI_COMM_WORLD);
            free(DS[i].data);
        }

        local_density = (unsigned int *)calloc(Nx * N * N, sizeof(unsigned int));

        local_n_points = DS[0].cur_points;

        local_points = DS[0].data;
        free(DS);
    }

    if (my_rank != 0)
    {
        local_density = (unsigned int *)calloc(Nx * N * N / n_processors, sizeof(unsigned int));

        MPI_Status status;
        MPI_Recv(&local_n_points, 1, MPI_UNSIGNED_LONG, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        local_points = (float *)malloc(local_n_points * sizeof(float) * 3);

        MPI_Recv(local_points, local_n_points * 3, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    // Density matrix computation
    size_t start = N / n_processors * my_rank + MIN(N % n_processors, my_rank);
    size_t Nx_range[2] = {start, start + Nx};

    for (size_t i = 0; i < local_n_points; i+=3) {
        float point[3]{local_points[3*i], local_points[3*i+1], local_points[3*i+2]};
        update_density_matrix(local_density, point, N, R, Nx_range);
    }

    // Write density matrix to file
    FILE* file = fopen("density.bin", "wb");

    // Read the 3d points
    size_t cells_per_write = 1024;
    size_t conversion_factor = sizeof(float);
    size_t block_size = cells_per_write * conversion_factor;

    fseek(file, N*N*Nx * conversion_factor, SEEK_SET);
    byte *buffer = create_empty_buffer(block_size);

    for (size_t i = 0; i < N * N * Nx; i += cells_per_write)
    {
        generic_convert(&local_density[i], cells_per_write);
        size_t n_bytes_read = write_bytes(file, buffer, block_size);
        size_t n_points_read = n_bytes_read / conversion_factor;

        assert(n_bytes_read == n_points_read * conversion_factor &&
                "Warning: The coordinates of the last point are missing,"
                "discarding them\n");

        actual_n_points += n_points_read;

        // Convert to float
        float *points = byte_to_float(buffer, n_bytes_read);

        // Insert points in the DS
        insert_points(points, n_points_read, DS, R, n_processors);
    }
}


float* get_cell_center(size_t Nx, size_t Ny, size_t Nz, size_t N) {
    float center[3];
    float conversion_factor = 1.0 / (float)N;
    center[0] = Nx * conversion_factor + 0.5 * N;
    center[1] = Ny * conversion_factor + 0.5 * N;
    center[2] = Nz * conversion_factor + 0.5 * N;

    return center;
}

float squared_point_distance(const float* point1, const float* point2) {
    float x = point1[0] - point2[0];
    float y = point1[1] - point2[1];
    float z = point1[2] - point2[2];

    return x*x + y*y + z*z;
}

void update_density_matrix(unsigned int* local_density, float* point, size_t N, float R, size_t* Nx_range) {
    for (size_t Nx = Nx_range[0]; Nx < Nx_range[1]; Nx++) {
        for(size_t Ny = 0; Ny < N; Ny++) {
            for(size_t Nz = 0; Nz < N; Nz++) {
                float* center = get_cell_center(Nx, Ny, Nz, N);
                float d_2 = squared_point_distance(center, point);

                if (d_2 < R*R) {
                    local_density[(Nx-Nx_range[0]) * N * N + Ny * N + Nz] += 1;
                }
            }
        }
    }
}
