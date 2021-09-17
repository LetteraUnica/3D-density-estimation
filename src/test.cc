#include <stdio.h>
// #include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory> // unique_ptr
#include <utility>
#include <string>
#include <cstring>
#include <cassert>

using byte = unsigned char;
constexpr size_t byte_size = sizeof(byte);

byte *float_to_byte(float *float_array, int n_elements)
{
    byte *byte_array = (byte *)malloc(n_elements * 4 * sizeof(byte));

    memcpy(byte_array, float_array, n_elements * sizeof(float_array[0]));

    return byte_array;
}

float *byte_to_float(byte *byte_array, int n_elements)
{
    assert(n_elements / 4 * 4 == n_elements);

    float *float_array = (float *)malloc(n_elements / 4 * sizeof(float));

    memcpy(float_array, byte_array, n_elements * sizeof(byte_array[0]));

    return float_array;
}

byte *create_empty_buffer(int n_bytes)
{
    return (byte *)malloc(n_bytes * sizeof(byte));
}

int read_bytes(FILE *file, byte *buffer, int n_bytes)
{
    return fread(buffer, sizeof(byte), n_bytes, file);
}

int write_bytes(FILE *file, byte *buffer, int n_bytes)
{
    return fwrite(buffer, sizeof(byte), n_bytes, file);
}

int main(int argc, char **argv)
{
    const int N = 1000000;

    FILE *ptr;
    // Opens a file in read binary mode
    ptr = fopen("test.bin", "rb");
}