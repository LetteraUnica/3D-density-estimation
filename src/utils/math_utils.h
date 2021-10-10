#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <omp.h>

template <typename T>
T clamp(T a, T min, T max)
{
    T res = a < min ? min : a;
    return res < max ? res : max;
}

template <typename T>
T clamp(T a, T *range)
{
    return clamp(a, range[0], range[1]);
}

void set_random_seed(mt19937 &rng, random_device &rd)
{
    rng.seed(rd());
}

size_t get_cell(const float x, const size_t N, size_t *__restrict Nx_range)
{
    return clamp((long)(N * x), (long)Nx_range[0], (long)Nx_range[1]);
}

float get_point(const size_t cell, const float N_inv)
{
    return (cell + 0.5f) * N_inv;
}

size_t low_bound(const float x, size_t N, const float N_inv, size_t *__restrict Nx_range)
{
    size_t Nx = get_cell(x, N, Nx_range);
    float point = get_point(Nx, N_inv);

    return clamp(Nx + (size_t)(x > point), Nx_range);
}

size_t high_bound(const float x, size_t N, const float N_inv, size_t *__restrict Nx_range)
{
    size_t Nx = get_cell(x, N, Nx_range);
    float point = get_point(Nx, N_inv);

    return clamp(Nx + (size_t)(x > point), Nx_range);
}

void update_density_matrix(u_int32_t *__restrict local_density, float *__restrict point,
                           register const size_t N, const size_t N2,
                           register const float N_inv, register const float R, const float R2,
                           size_t *Nx_range, size_t *default_range)
{
    size_t register low_x = low_bound(point[0] - R, N, N_inv, Nx_range);
    size_t high_x = high_bound(point[0] + R, N, N_inv, Nx_range);
    for (; low_x < high_x; low_x++)
    {
        size_t register a = (low_x - Nx_range[0]) * N2;

        float x = get_point(low_x, N_inv) - point[0];
        float x2 = x * x;
        float radius_x2 = R2 - x2;
        float register radius_x = sqrtf(radius_x2);

        size_t register low_y = low_bound(point[1] - radius_x, N, N_inv, default_range);
        size_t high_y = high_bound(point[1] + radius_x, N, N_inv, default_range);
        for (; low_y < high_y; low_y++)
        {
            size_t register b = a + low_y * N;

            float y = get_point(low_y, N_inv) - point[1];
            float y2 = y * y;
            float register radius_xy = sqrtf(radius_x2 - y2);

            size_t register low_z = low_bound(point[2] - radius_xy, N, N_inv, default_range);
            size_t high_z = high_bound(point[2] + radius_xy, N, N_inv, default_range);
            for (; low_z < high_z; low_z++)
            {
                local_density[b + low_z] += 1lu;
            }
        }
    }
}

void update_density_matrix_OMP(u_int32_t *__restrict local_density, float *__restrict point,
                               register const size_t N, const size_t N2,
                               register const float N_inv, register const float R, const float R2,
                               size_t *Nx_range, size_t *default_range)
{
    size_t register low_x = low_bound(point[0] - R, N, N_inv, Nx_range);
    size_t high_x = high_bound(point[0] + R, N, N_inv, Nx_range);
    for (; low_x < high_x; low_x++)
    {
        size_t register a = low_x * N2;

        float x = get_point(low_x, N_inv) - point[0];
        float x2 = x * x;
        float radius_x2 = R2 - x2;
        float register radius_x = sqrtf(radius_x2);

        size_t register low_y = low_bound(point[1] - radius_x, N, N_inv, default_range);
        size_t high_y = high_bound(point[1] + radius_x, N, N_inv, default_range);
        for (; low_y < high_y; low_y++)
        {
            size_t register b = a + low_y * N;

            float y = get_point(low_y, N_inv) - point[1];
            float y2 = y * y;
            float register radius_xy = sqrtf(radius_x2 - y2);

            size_t register low_z = low_bound(point[2] - radius_xy, N, N_inv, default_range);
            size_t high_z = high_bound(point[2] + radius_xy, N, N_inv, default_range);
            for (; low_z < high_z; low_z++)
            {
                local_density[b + low_z] += 1lu;
            }
        }
    }
}

#endif