#ifndef ARRAY_H_
#define ARRAY_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <random>

using namespace std;

namespace array
{
    void rand_array(float *array, size_t n_elements, mt19937 &rng,
                    uniform_real_distribution<float> &udist)
    {
        for (size_t i = 0; i < n_elements; i++)
        {
            array[i] = udist(rng);
        }
    }

    template <typename T>
    void fill_array(T *array, T value, size_t n_elements)
    {
        for (size_t i = 0; i < n_elements; i++)
        {
            array[i] = value;
        }
    }

    template <typename T>
    void print_array(T *array, size_t start, size_t end) {
        for (; start < end; start++)
        {
            printf("%f ", array[start]);
        }
        printf("\n");
    }

    struct ResizableArray
    {
        size_t max_points;
        size_t cur_points;
        float *data;
    };

    ResizableArray create_resizable_array(size_t max_points)
    {
        ResizableArray array{max_points, 0, (float *)malloc(3 * max_points * sizeof(float))};
        return array;
    }

    void resize_array(ResizableArray *array, size_t new_points)
    {
        array->max_points = new_points;
        array->data = (float *)realloc(array->data, 3 * new_points * sizeof(float));
    }

    void fill_resizable_array(ResizableArray *array, size_t start, size_t end,
                              float fill_value)
    {
        for (; start < end; start++)
        {
            array->data[start] = fill_value;
        }
    }

    ResizableArray *create_empty_data_structure(int n_processors, size_t max_points)
    {
        ResizableArray *data_structure = (ResizableArray *)malloc(n_processors * sizeof(ResizableArray));
        for (int i = 0; i < n_processors; i++)
        {
            data_structure[i] = create_resizable_array(max_points);
        }

        return data_structure;
    }

    void insert_point(ResizableArray *DS, float x, float y, float z, float R, int n_processors)
    {
        float step = 1.f / (float)n_processors;
        float low = 0.f;
        float high = low + step;

        for (int j = 0; j < n_processors; j++)
        {
            if (x >= low - R && x < high + R)
            {
                if (DS[j].max_points < DS[j].cur_points + 1)
                {
                    resize_array(&DS[j], DS[j].max_points * 2);
                }

                DS[j].data[3 * DS[j].cur_points] = x;
                DS[j].data[3 * DS[j].cur_points + 1] = y;
                DS[j].data[3 * DS[j].cur_points + 2] = z;
                DS[j].cur_points += 1;
            }

            low = high;
            high += step;
        }
    }

    void insert_points(float *points, size_t n_points, ResizableArray *DS, float R, int n_processors)
    {
    for (size_t i = 0; i < n_points; i++)
    {
        insert_point(DS, points[3 * i], points[3 * i + 1], points[3 * i + 2], R, n_processors);
    }
}
}

#endif