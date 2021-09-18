#include <stdio.h>
// #include <mpi.h>
#include <stdlib.h>
#include <string.h>

struct ResizableArray
{
    size_t N;
    float *data = (float *)malloc(N * sizeof(float));
};

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