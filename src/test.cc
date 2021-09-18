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

void rand_array(float *array, size_t n_elements, size_t seed)
{
    for (size_t i = 0; i < n_elements; i++) {
        array[i] = udist(rng);
    }
}

void set_random_seed() {
    rng.seed(rd());
}


int main(int argc, char **argv)
{
    const int N = 10;

    float* float_array = (float *)malloc(N * sizeof(float));
    rand_array(float_array, N, 1697);
    FILE *ptr;
    // Opens a file in write binary mode
    ptr = fopen("test.bin", "wb+");
    fseek(ptr, 0, SEEK_SET);

    byte* bytes = float_to_byte(float_array, N);
    fwrite(bytes, 1, 4*N, ptr);

    byte* buffer = create_empty_buffer(4*N);
    fseek(ptr, 0, SEEK_SET);
    fread(buffer, 1, 4*N, ptr);

    float* new_float_array = byte_to_float(buffer, N*4);


    // float *float_array = (float *)malloc(N * sizeof(float));
    // memcpy(float_array, int_array, N * sizeof(int_array[0]));
    
    for (int i = 0; i < N; i++){
        printf("%f  %f\n", float_array[i], new_float_array[i]);
    }
}