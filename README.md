# Phase field fracture (PFF) 


## About


This is the source code of phase field fracture (or crack) simulation program written in C++ and makes use of [deal.II](https://dealii.org/) finite element library.

## Program usage


##### Running the program


1. `mkdir debug`
2. `mkdir results`
3. `rm -rf CMakeCache.txt CMakeFiles`
4. `cmake .`
5. `make release`
6. `make`
7. `mpirun -np N ./pfc`

##### Rerunning program with only parameter changes


1. `make runclean`
2. `make release`
3. `make`
4. `mpirun -np N ./pfc`

## References

## Acknowledgement


