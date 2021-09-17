#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <memory> // unique_ptr
#include <utility>
#include <string>
#include <cstring>

int get_number_of_points(std::iostream file) {
    file.read()
}

float *read_input(char *filename)
{
}

template <typename T>
void print_vec(std::vector<T> *vector)
{
    int length = sizeof(vector) / sizeof(vector[0]);
    for (int i = 0; i < length; i++)
    {
        printf("%f ", vector[i]);
    }
    printf("\n");
}

int main(int argc, char **argv)
{
    int N = 1000;
    std::unique_ptr<std::vector<double>> vector{new std::vector<double>(N)};
}