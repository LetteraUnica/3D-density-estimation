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

    long Nl = (long)N;
    cell[0] = clamp((long)(N * x), (long)Nx_range[0], (long)Nx_range[1]);
    cell[1] = clamp((long)(N * y), 0l, Nl);
    cell[2] = clamp((long)(N * z), 0l, Nl);

    return cell;
}

float squared_distance(const float *point1, const float *point2)
{
    float x = point1[0] - point2[0];
    float y = point1[1] - point2[1];
    float z = point1[2] - point2[2];

    return x * x + y * y + z * z;
}

void set_random_seed(mt19937 &rng, random_device &rd)
{
    rng.seed(rd());
}

void update_density_matrix(u_int32_t *local_density, float *point, size_t N, float R, size_t *Nx_range)
{
    float R2 = R * R;
    R = R + 1;

    size_t *lows = point_to_cell(point[0] - R, point[1] - R, point[2] - R, N, Nx_range);
    size_t *highs = point_to_cell(point[0] + R, point[1] + R, point[2] + R, N, Nx_range);

    for (size_t Nx = lows[0]; Nx < highs[0]; Nx++)
    {
        size_t a = (Nx - Nx_range[0]) * N * N;
        for (size_t Ny = lows[1]; Ny < highs[1]; Ny++)
        {
            size_t b = a + Ny * N;

            for (size_t Nz = lows[2]; Nz < highs[2]; Nz++)
            {
                float *center = cell_to_point(Nx, Ny, Nz, N);
                float d2 = squared_distance(center, point);

                local_density[b + Nz] += (u_int32_t)(d2 < R2);

                free(center);
            }
        }
    }

    free(lows);
    free(highs);
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

void fast_update_density_matrix(u_int32_t *local_density, float *point, size_t N, float N_inv, float R, size_t *Nx_range)
{
    float R2 = R * R;

    size_t default_range[2]{0lu, N};

    size_t low_x = low_bound(point[0] - R, N, N_inv, Nx_range);
    size_t high_x = high_bound(point[0] + R, N, N_inv, Nx_range);
    for (size_t Nx = low_x; Nx < high_x; Nx++)
    {
        size_t a = (Nx - Nx_range[0]) * N * N;

        float x = get_point(Nx, N_inv);
        float x2 = (x - point[0]) * (x - point[0]);
        float radius_x2 = R2 - x2;
        float radius_x = sqrtf(radius_x2);

        size_t low_y = low_bound(point[1] - radius_x, N, N_inv, default_range);
        size_t high_y = high_bound(point[1] + radius_x, N, N_inv, default_range);
        for (size_t Ny = low_y; Ny < high_y; Ny++)
        {
            size_t b = a + Ny * N;

            float y = get_point(Ny, N_inv);
            float y2 = (y - point[1]) * (y - point[1]);
            float radius_xy = sqrtf(radius_x2 - y2);

            size_t low_z = low_bound(point[2] - radius_xy, N, N_inv, default_range);
            size_t high_z = high_bound(point[2] + radius_xy, N, N_inv, default_range);
            for (size_t Nz = low_z; Nz < high_z; Nz++)
            {
                local_density[b + Nz] += 1;
            }
        }
    }
}

void fast_update_density_matrix_omp(u_int32_t *local_density, float *point, size_t N, float N_inv, float R, size_t *Nx_range)
{
    float R2 = R * R;

    size_t default_range[2]{0lu, N};

    size_t low_x = low_bound(point[0] - R, N, N_inv, Nx_range);
    size_t high_x = high_bound(point[0] + R, N, N_inv, Nx_range);
    for (size_t Nx = low_x; Nx < high_x; Nx++)
    {
        size_t a = Nx * N * N;

        float x = get_point(Nx, N_inv);
        float x2 = (x - point[0]) * (x - point[0]);
        float radius_x2 = R2 - x2;
        float radius_x = sqrtf(radius_x2);

        size_t low_y = low_bound(point[1] - radius_x, N, N_inv, default_range);
        size_t high_y = high_bound(point[1] + radius_x, N, N_inv, default_range);
        for (size_t Ny = low_y; Ny < high_y; Ny++)
        {
            size_t b = a + Ny * N;

            float y = get_point(Ny, N_inv);
            float y2 = (y - point[1]) * (y - point[1]);
            float radius_xy = sqrtf(radius_x2 - y2);

            size_t low_z = low_bound(point[2] - radius_xy, N, N_inv, default_range);
            size_t high_z = high_bound(point[2] + radius_xy, N, N_inv, default_range);
            for (size_t Nz = low_z; Nz < high_z; Nz++)
            {
                local_density[b + Nz] += 1;
            }
        }
    }
}

#endif