#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>

#include <random>
#include <iostream>

#include "utils/io_utils.h"

typedef u_int8_t byte;

using namespace array;
using namespace std;

int main(int argc, char **argv)
{
    std::random_device rd;
    std::mt19937 rng;
    std::uniform_real_distribution<float> udist(0, 1);

    if (argc < 2)
    {
        printf("Error: You must provide the number of points to be generated\n");
        return 1;
    }

    const u_int32_t n_points = atol(argv[1]);
    
    random_input_file("points.bin", n_points, rng, udist);
    return 0;
}