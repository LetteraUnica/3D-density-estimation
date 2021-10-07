#ifndef IO_UTILS_H_
#define IO_UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "array.h"

#include <random>

typedef u_int8_t byte;
using namespace std;


template <typename O, typename I>
O *generic_convert(const I *input, size_t n_bytes)
{
    O *output = (O *)malloc(n_bytes);
    memcpy(output, input, n_bytes);

    return output;
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

byte *read_input_file(FILE *file, size_t n_bytes)
{
    byte *buffer = create_empty_buffer(n_bytes);
    read_bytes(file, buffer, n_bytes);

    return buffer;
}

u_int32_t get_number_of_points(FILE *file)
{
    byte *buffer = read_input_file(file, 4);
    u_int32_t n_points;
    memcpy(&n_points, buffer, 4);
    free(buffer);

    return n_points;
}

void create_input_file(const char *filename, int n_points, mt19937 &rng,
                       uniform_real_distribution<float> &udist)
{
    FILE *ptr = fopen(filename, "rb");
    write_bytes(ptr, (byte *)&n_points, sizeof(int));

    float *random_points = (float *)malloc(3 * sizeof(float));
    for (int i = 0; i < n_points; i++)
    {
        array::rand_array(random_points, 3, rng, udist);
        write_bytes(ptr, (byte *)random_points, 3 * sizeof(float));
    }

    fclose(ptr);
    free(random_points);
}


#endif