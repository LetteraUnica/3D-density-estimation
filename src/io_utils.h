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

size_t get_file_size(FILE *file)
{
    size_t previous_position = ftell(file);
    fseek(file, 0L, SEEK_END);
    size_t size = ftell(file);
    fseek(file, previous_position, SEEK_SET);

    return size;
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