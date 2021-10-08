#!/bin/bash

flags="-Wall -Wextra -O3 -march=native"

g++ $flags -o generate_points.x generate_points.cc
g++ $flags -o density_serial.x density_serial.cc

mpicxx $flags -o density_mpi.x density_mpi.cc

#g++ $flags -fopenmp -o density_omp.x density_omp.cc