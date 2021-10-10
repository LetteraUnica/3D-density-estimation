#!/bin/bash

R=0.261

echo mpi:
rm -rf density.bin
mpirun -np 4 ./density_mpi.x 20 $R points.bin 
sha256sum density.bin 

echo serial:
rm -rf density.bin
./density_serial.x 20 $R points.bin 
sha256sum density.bin

echo omp:
rm -rf density.bin
./density_omp.x 20 $R points.bin 
sha256sum density.bin 