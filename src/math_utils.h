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

    cell[0] = clamp((size_t)abs(N * x), Nx_range[0], Nx_range[1]);
    cell[1] = clamp((size_t)abs(N * y), 0lu, N);
    cell[2] = clamp((size_t)abs(N * z), 0lu, N);

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

void update_density_matrix(unsigned int *local_density, float *point, size_t N, float R, size_t *Nx_range)
{
    float R_2 = R * R;

    size_t *lows = point_to_cell(point[0] - R, point[1] - R, point[2] - R, N, Nx_range);
    size_t *highs = point_to_cell(point[0] + R, point[1] + R, point[2] + R, N, Nx_range);
    printf("%f, %f, %f,,, ", point[0], point[1], point[2]);
    printf("%ld, %ld, %ld,,, %ld, %ld, %ld\n", lows[0], lows[1], lows[2], highs[0], highs[1], highs[2]);

    for (size_t Nx = lows[0]; Nx < highs[0]; Nx++)
    {
        size_t a = (Nx - Nx_range[0]) * N * N;
        for (size_t Ny = lows[1]; Ny < highs[1]; Ny++)
        {
            a += Ny * N;
            for (size_t Nz = lows[2]; Nz < highs[2]; Nz++)
            {
                float *center = cell_to_point(Nx, Ny, Nz, N);
                float d_2 = squared_distance(center, point);

                if (d_2 < R_2)
                {
                    local_density[a + Nz] += 1;
                }

                free(center);
            }
        }
    }

    free(lows);
    free(highs);
}

#endif