#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <random>

typedef u_int8_t byte;

byte *float_to_byte(float *float_array, int n_elements)
{
    byte *byte_array = (byte *)malloc(n_elements * sizeof(float) * sizeof(byte));

    memcpy(byte_array, float_array, n_elements * sizeof(float));

    return byte_array;
}

float *byte_to_float(byte *byte_array, int n_bytes)
{
    size_t float_size = sizeof(float);
    assert(n_bytes / float_size * float_size == n_bytes);

    float *float_array = (float *)malloc(n_bytes);

    memcpy(float_array, byte_array, n_bytes);

    return float_array;
}

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