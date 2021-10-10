# How to compile
The code can be compiled by running `./compile` or `sh compile`.  
The script will produce 3 executables:  
1. `density_omp.x` Is the OpenMP version of the code  
2. `density_mpi.x` Is the MPI version of the code  
3. `density` The MPI version of the code, just named differently  
These executables when run correctly will produce as output a file named density.bin which is the 3D density matrix of the point distribution.  

Note: The script will probably trigger some warnings, however these could be safely ignored, they are caused by compiling pragmas without openmp or comparing different integer types.


# Performance evaluation
## Introduction
The OpenMP code was compiled with GCC version 10.3.0 while the MPI code was compiled with MPICH version 3.3.2

I tested the code using 4 MPI processes for the MPI version and 4 omp threads for the OpenMP version, this number corresponds to the number of physical threads of my CPU, an Intel(R) Core(TM) i5-4690 CPU @ 3.50GHz

## Elapsed time per n points
I fixed the grid number to $N=512$, the radius to $R=1/n^{1/3}$ and I varied the number of points from $10^6$ to $10^7$

![](Images/n_points_linear.png)

We can see that the elapsed time scales linearly with the number of points, this is expected because we didn't use a brute force approach but we only check the grid points that are close to a certain point. The overall complexity of this algorithm, without counting communication and I/O, is $\Theta(N^3 R^3 n)$, since $R=1/n^{1/3}$ we should have $\Theta(N^3)$ which is a constant. However this doesn't account for the time to read the input file and the time to communicate the points around, which are both operations that scale with $\Theta(n)$. So we improve the speed of the algorithm but it is still linear in the number of points.

## Elapsed time per grid number
I fixed the radius to $R=0.05$, the number of points to $n=2^{20} \approx 10^{6}$ and I varied the grid number from $16$ to $256$

![](Images/grid_number_linear.png)

We can see that the elapsed time scales with the cube of the grid number, as expected.
The green line represents a fitted polynomial of degree 3  

## Comments on MPI vs OpenMP
MPI is faster than OpenMP in these tests, this could be caused by the fact that the OpenMP program doesn't parallelize the I/O operations while the MPI program does, this could also be seen in the elapsed time per grid number graph, where the OpenMP program is faster than the MPI version up until $N\sim 80$, while for $N > 80$ the MPI version is faster