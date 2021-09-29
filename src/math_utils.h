#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

template <typename T>
T clamp(T a, T min, T max)
{
    T res = a < min ? min : a;
    return res < max ? res : max;
}

float *cell_to_point(size_t Nx, size_t Ny, size_t Nz, size_t N)
{
    float *center = (float *)malloc(3 * sizeof(float));
    float conversion_factor = 1.0f / (float)N;

    center[0] = (Nx + 0.5f) * conversion_factor;
    center[1] = (Ny + 0.5f) * conversion_factor;
    center[2] = (Nz + 0.5f) * conversion_factor;

    return center;
}

size_t *point_to_cell(float x, float y, float z, size_t N, size_t *Nx_range)
{
    size_t *cell = (size_t *)malloc(3 * sizeof(size_t));

    cell[0] = clamp((size_t)(N * x), Nx_range[0], Nx_range[1]);
    cell[1] = clamp((size_t)(N * y), 0lu, N);
    cell[2] = clamp((size_t)(N * z), 0lu, N);

    return cell;
}

float squared_distance(const float *point1, const float *point2)
{
    float x = point1[0] - point2[0];
    float y = point1[1] - point2[1];
    float z = point1[2] - point2[2];

    return x * x + y * y + z * z;
}

#endif