## How to compile
The code can be compiled by running `./compile` or `sh compile`.  
The script will produce 4 executables: 
* `density_serial.x` Is the serial version of the code, mainly useful for testing purpuses
* `density_omp.x` Is the OpenMP version of the code
* `density_mpi.x` Is the MPI version of the code
* `density` The MPI version of the code, just named differently

The script will probably trigger some warnings, however these could be safely ignored, they are caused by compiling pragmas without openmp or comparing different integer types.


## Performance evaluation
### Introduction
The OpenMP code was compiled with GCC version 10.3.0 while the MPI code was compiled with MPICH version 3.3.2

I tested the code using 4 MPI processes for the MPI version and 4 omp threads for the OpenMP version, this number corresponds to the number of physical threads of my CPU, an Intel(R) Core(TM) i5-4690 CPU @ 3.50GHz

### Elapsed time per n points
I fixed the grid number to $N=512$, the radius to $R=1/n^{1/3}$ and I varied the number of points from $10^6$ to $10^7$

![](n_points_linear.png)

We can see that the elapsed time scales linearly with the number of points, this is expected because we didn't use a brute force approach but we only check the grid points that are close to a certain point. The overall complexity of this algorithm, without counting communication and I/O, is $\Theta(N^3 R^3 n)$, since $R=1/n^{1/3}$ we should have $\Theta(N^3)$ which is a constant. However this doesn't account for the time to read the input file and the time to communicate the points around, which are both operations that scale with $\Theta(n)$. So we improve the speed of the algorithm but it is still linear in the number of points.

### Elapsed time per grid number
I fixed the radius to $R=0.05$, the number of points to $n=2**20 = 1048576$ and I varied the grid number from $16$ to $256$

![](grid_number_linear.png)

We can see that the elapsed time scales with the cube of the grid number, as expected.
The green line represents the function