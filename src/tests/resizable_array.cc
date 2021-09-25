#include <stdio.h>
// #include <mpi.h>
#include <stdlib.h>
#include <string.h>

struct ResizableArray
{
    size_t max_size;
    size_t cur_size;
    float *data = (float *)malloc(max_size * sizeof(float));
};

void resize_array(ResizableArray *array, size_t new_size)
{
    array->max_size = new_size;
    array->data = (float *)realloc(array->data, new_size * sizeof(float));
}

void fill_resizable_array(ResizableArray *array, size_t start, size_t end,
                          float fill_value)
{
    for (; start < end; start++)
    {
        array->data[start] = fill_value;
    }
}

int main(int argc, char **argv)
{

    int n_processors = 4;
    ResizableArray *data_structure = (ResizableArray *)malloc(n_processors * sizeof(ResizableArray));

    printf("%ld\n", sizeof(ResizableArray));

    size_t n_elements = 10;
    for (int i = 0; i < n_processors; i++)
        data_structure[i] = ResizableArray{n_elements};

    for (int j = 0; j < n_processors; j++)
    {
        for (int i = 0; i < n_elements; i++)
        {
            printf("%f ", data_structure[j].data[i]);
        }
        printf("\n");
    }
    return 0;
}