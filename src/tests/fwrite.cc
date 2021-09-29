#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <random>

typedef u_int8_t byte;

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

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))


int main(int argc, char **argv) {
    size_t N = 1000000;
    FILE* file = fopen("prova.bin", "wb");
    float* float_array = (float*) malloc(N * sizeof(float));
    // byte* byte_array = create_empty_buffer(N*sizeof(float));
    // write_bytes(file, byte_array, N*sizeof(float));

    size_t cells_per_write = 1024;
    size_t conversion_factor = sizeof(float);
    size_t block_size = cells_per_write * conversion_factor;

    for (size_t i = 0; i < N; i += cells_per_write)
    {
        size_t current_block_size = MIN(block_size, (N - i) * conversion_factor);
        size_t n_bytes = write_bytes(file, (byte *)&float_array[i], current_block_size);
    }

    fclose(file);

    file = fopen("prova.bin", "rb");

    byte* byte_array = create_empty_buffer(N*sizeof(float));
    read_bytes(file, byte_array, 4*N);

    size_t sum = 0;
    for (size_t i = 0; i < 4*N; i++) {
        sum += byte_array[i];
    }
    printf("%ld\n", sum);

    return 0;
}