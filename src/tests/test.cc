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

void print_array(float* array, int n_elements) {
    for (size_t i = 0; i < n_elements; i++){
        printf("%f, %p\t", array[i], &array[i]);
    }
    printf("\n\n");
}

void print_array(byte* array, int n_elements) {
    for (size_t i = 0; i < n_elements; i++){
        printf("%c, %p\t", array[i], &array[i]);
    }
    printf("\n\n");
}


int main(int argc, char **argv)
{
    const int N = 10;

    float* float_array = (float *)malloc(N * sizeof(float));
    rand_array(float_array, N, 1697);
    
    print_array(float_array, N);

    byte* byte_array = (byte*) float_array;

    print_array(byte_array, N*4);

    print_array((float*)byte_array, N);
}