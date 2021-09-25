#include <stdio.h>
// #include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory> // unique_ptr
#include <utility>
#include <string>
#include <cstring>
#include <cassert>

typedef u_int8_t byte;

byte *float_to_byte(float *float_array, int n_elements)
{
    byte *byte_array = (byte *)malloc(n_elements * sizeof(float) * sizeof(byte));

    memcpy(byte_array, float_array, n_elements * sizeof(float));

    return byte_array;
}

float *byte_to_float(byte *byte_array, int n_elements)
{
    assert(n_elements / 4 * 4 == n_elements);

    float *float_array = (float *)malloc(n_elements / 4 * sizeof(float));

    memcpy(float_array, byte_array, n_elements * sizeof(byte_array[0]));

    return float_array;
}

byte *create_empty_buffer(size_t n_bytes)
{
    return (byte *)malloc(n_bytes * sizeof(byte));
}

size_t read_bytes(FILE *file, byte *buffer, size_t n_bytes)
{
    return fread(buffer, sizeof(byte), n_bytes, file);
}

size_t write_bytes(FILE *file, byte *buffer, size_t n_bytes)
{
    return fwrite(buffer, sizeof(byte), n_bytes, file);
}

void fill_array(float *array, float value, size_t n_elements)
{
    for(int i = 0; i < n_elements; i++){
        array[i] = value;
    }
}

std::random_device rd;
std::mt19937 rng;
std::uniform_real_distribution<float> udist(0, 1);

float* rand_array(size_t n_elements)
{
    float* array = (float *)malloc(n_elements * sizeof(float));
    for (size_t i = 0; i < n_elements; i++) {
        array[i] = udist(rng);
    }

    return array;
}

void set_random_seed() {
    rng.seed(rd());
}

void create_input_file(const char *filename, int n_points) {
    set_random_seed();

    FILE* ptr = fopen("test.bin", "rb");
    write_bytes(ptr, (byte*)&n_points, 4);

    for (int i = 0; i < n_points; i++) {
        float* random_points = rand_array(3);
        byte* buffer = float_to_byte(random_points, 3);
        write_bytes(ptr, buffer, 3 * sizeof(float));
    }
}

float* read_input_file(const char *filename);

float* generate_random_3d_points(size_t n_points);


struct ResizableArray {
    size_t N;
    float* data;
};

int main(int argc, char **argv)
{
    if (argc < 4) {
        printf("Error: you must provide a value for the grid number N and the radius R");
        return 1;
    }

    unsigned long N = atol(argv[1]);
    float R = (float)atof(argv[2]);

    int n_processors = 4;
    ResizableArray* data_structure = (ResizableArray*) malloc(n_processors*sizeof(ResizableArray));

    if (argc==4) {
        // Read input file
        FILE* file = fopen("test.bin", "rb");

        byte buffer[4];
        int n_points;
        fread(buffer, 1, 4, file);
        memcpy(&n_points, buffer, 4);

        size_t points_per_read = 1024;
        size_t conversion_factor = sizeof(float)*3;
        size_t block_size = points_per_read*conversion_factor;
        size_t actual_n_points = 0;

        for (size_t i=0; i<n_points; i+=block_size) {
            byte* buffer = create_empty_buffer(block_size);
            
            size_t n_bytes_read = read_bytes(file, buffer, block_size);
            size_t n_points_read = n_bytes_read / conversion_factor;

            if (n_bytes_read != n_points_read*conversion_factor) {
                printf("Warning: The coordinates of the last point"
                "are missing, discarding them\n");
            }
            actual_n_points += n_points_read;

            // Process the buffer
            float* points = byte_to_float(buffer, n_points_read);

        }   

        if (actual_n_points !=n_points) {
            printf("Warning: The number of points found in the file are"
            "different than what is specified in the file header,"
            "specified: %ld, found: %ld, the program will continue with"
            "the number of found points\n", n_points, actual_n_points);

            n_points = actual_n_points;
        }

        n_po

    }

    else {
        // Generate n points
        unsigned long n = 1000;
    }
}