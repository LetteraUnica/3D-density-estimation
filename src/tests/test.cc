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

template <typename T>
T clamp(T a, T min, T max)
{
    T res = a < min ? min : a;
    return res < max ? res : max;
}

template <typename T>
T clamp(T a, T* range)
{
    return clamp(a, range[0], range[1]);
}

size_t get_cell(float x, size_t N, size_t *Nx_range)
{
    return clamp((long)(N * x), (long)Nx_range[0], (long)Nx_range[1]);
}

float get_point(size_t cell, float N_inv)
{
    return (cell + 0.5f) * N_inv;
}

size_t low_bound(float x, size_t N, float N_inv, size_t *Nx_range)
{
    size_t Nx = get_cell(x, N, Nx_range);
    float point = get_point(Nx, N_inv);

    return clamp(Nx + (size_t)(x > point), Nx_range);
}

size_t high_bound(float x, size_t N, float N_inv, size_t *Nx_range)
{
    size_t Nx = get_cell(x, N, Nx_range);
    float point = get_point(Nx, N_inv);

    return clamp(Nx + (size_t)(x > point), Nx_range);
}

int main(int argc, char **argv)
{
    size_t N = 5;
    float N_inv = 1.f/N;
    size_t Nx_range[2]{0lu, N};
    float x = atof(argv[1]);
    float R = atof(argv[2]);

    printf("x: %f, N: %ld, R: %f, low: %ld, high: %ld", x, N, R, low_bound(x-R, N, N_inv, Nx_range), high_bound(x+R, N, N_inv, Nx_range));

    return 0;
}