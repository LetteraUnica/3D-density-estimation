# FHPC project report

## Contents
* Problem definition
* Single core solution
* Parallel implementation
* Summary
* Further development

## Problem definition
The problem is a specific version of the [range-searching problem](https://en.wikipedia.org/wiki/Range_searching), which in the most general case is defined as:  
> The range searching problem most generally consists of preprocessing a set $S$ of objects, in order to determine which objects from $S$ intersect with a query object, called a range

In our case we have an input distribution of points in 3D and we are asked to compute the local density at the node points of a regular 3D grid with grid number $N$, the density is defined as the number of points that lie within a sphere of radius $R$ centred at every grid node point; the grid number and the radius are parameters acquired at the command line.  

We can see that our case considers $S$ to be a distribution of $n$ points in $\R^3$ in which each coordinate is constrained to be in $[0, 1]$. We can rigously define the problem as follows:
* $S = \{x_i | x_i\in \R^3 \;\; i=1,...,n\}$
* $x_{ij} \in [0,1] \;\; \forall \; i,j$  

While the query object is a sphere centred at every node point

## Single core solution
### Brute force
Count for every node point the number of points that are within a radius R from it, this requires $\Theta(N^3n)$ operations.  

Pros:
* Simple to implement
* No cache misses

Cons:
* $\Theta(N^3n)$ complexity

### K-d tree
Use a k-d tree to preprocess the data, this approach has a complexity of $\Theta(\log n)$ per each grid point, for a total complexity of $\Theta(N^3\log n)$, however it has a lot of cache misses.  

Pros:  
* $\Theta(N^3\log n)$ complexity  

Cons:  
* Lots of cache misses ($\Theta(\log n)$ per node)  
* Difficult to implement

### Restrict search space  
This method does the opposite: instead of counting the number of $x_i$ "near" a grid point simply count the number of grid points near a $x_i$, that is $D_{ijk} = D_{ijk} + ||N_{ijk} - x_l||_2 < R$.  
At first glance this method is the same as the brute force one, however we can only update the part of the density $D$ which needs to be updated, resulting in a complexity of $\Theta((NR)^3 n)$, because the ratio of the volume occupied by the sphere goes like $\Theta(R^3)$ and the total number of grid points is $N^3$.  
Note that $R<1$, so this is better than the brute force method and the lower $R$ the more efficient it becomes.

Pros:  
* $\Theta((NR)^3 n)$ complexity, lower than brute force and for low $R$ better than kd-trees
* No cache misses if programmed correctly

Cons:
* Harder to implement than brute force but simpler than kd-trees
* Needs to keep the whole $D$ matrix in memory

## Parallel implementation
### Introduction
The above analysis assumes that all the $x_i$ are in the memory and sometimes even that the density matrix $D$ is in memory as well. In a parallel implementation we would like to avoid this, because we want each core to be as independent as possible and we want to minimize the amount of communication between them. Below we present some possible approaches to perform the domain decomposition for this problem

### Domain decomposition
Here we have two possible approaches:  
#### Split the $x_i$  
This method has a clear advantage which is the reduced memory usage, furthermore it improves the performance of all algorithms because $n \rightarrow n/p$ where $p$ is the number of processors used. A disadvantage is that the points need to first be sorted by performing a bucket sort, which is a single core operation  

1. Read the $x_i$ from the file ($t_{seek} + n \cdot t_{read}$)
2. Split the $x_i$ (can be done at the same time as 1.)
3. Communicate the splits to all cores ($n \cdot t_{comm}$)
4. Let each core compute its density matrix $D$, complexity:
    - Brute-force $\Theta(N^3n / p^2)$
    - K-d tree $\Theta(N^3\log(n/p)/p)$
    - Restrict search space $\Theta((NR)^3 n/p)$
5. Print $D$ to the file ($p \cdot t_{seek} + N^3 \cdot t_{write}$)  
* Space complexity: 
    - Brute-force $\Theta(batch \cdot n/p)$
    - K-d tree $\Theta(batch \cdot n/p)$
    - Restrict search space $\Theta(n/p + N^3/p)$
  
#### Don't split the $x_i$  
This method is simpler than the first, however it needs a copy of the $x_i$ and the full density matrix $D$ in each core, the main advantage of the method is that the cores can start computing as soon as they have finished reading from the file, however more compute is later needed to merge the local matrices $D$ into a global one

1. Read the $x_i$ from the file ($t_{seek} + n \cdot t_{read}$)
2. Let each core compute its density matrix $D$, complexity:
    - Brute-force $\Theta(N^3n / p)$
    - K-d tree $\Theta(N^3\log(n/p))$
    - Restrict search space $\Theta((NR)^3 n/p)$
3. Sum all the density matrices ($N^3 \log p \cdot t_{comm} + N^3 \log p \cdot t_{sum}$)
4. Print $D$ to the file ($p \cdot t_{seek} + N^3 \cdot t_{write}$)

* Space complexity: $\Theta(n/p + pN^3)$ 


## Summary
The first method that I will implement is restrict search space with $x_i$ splitting in MPI, I will make some assumptions about the machine the method is implemented on:  

* The memory can contain the matrix $D$ and all the points $n$
* A single core can saturate all the bandwidth of the disk read and write

## Further development
Below I list some features to add to this algorithm in order of difficulty  

* Allow for D matrices bigger than the available memory 
* Allow for $n$ points bigger than the available memory 
* Implement the k-d tree algorithm