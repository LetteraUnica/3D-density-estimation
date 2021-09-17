#include <stdio.h>
// #include <mpi.h>
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
#include <cassert>


template <typename T>
float *write_vec(std::vector<T> *vector, std::string filename)
{   
    ofstream file;
    file.open(filename, ios::out | ios::binary);
    memblock = new char
}

template <typename T>
float *read_vec(std::vector<T> *vector, std::string filename)
{   
    ofstream file;
    file.open(filename, ios::out | ios::binary);
}

template <typename T>
void print_vec(std::vector<T> *vector, std::iostream &stream = std::cout)
{
    for (int i = 0; i < vector.size(); i++)
    {
        stream << vector[i];
    }
    stream << std::endl;
}

template <typename T>
void fill_vec(std::vector<T> *vector, T value)
{
    for (int i = 0; i < vector.size(); i++)
    {
        vector[i] = value;
    }
}

using byte = unsigned char;
using namespace std;

byte* float_to_byte(float* float_array, int n_elements) {
    byte* byte_array = (byte*)malloc(n_elements*4*sizeof(byte));

    std::memcpy(byte_array, float_array, sizeof(float_array));

    return byte_array;
}

float* byte_to_float(byte* byte_array, int n_elements) {
    assert(n_elements / 4 * 4 == n_elements);

    float* float_array = (float*)malloc(n_elements/4*sizeof(float));

    std::memcpy(float_array, byte_array, sizeof(byte_array));

    return float_array;
}


int main(int argc, char **argv)
{
}