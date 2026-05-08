# Phase field fracture (PFF) 


## About


This is the source code of phase field fracture (or crack) simulation program written in C++ and makes use of [deal.II](https://dealii.org/) finite element library.

## Program usage

## File structure

## Installing Dealii


##### Running the program

```
mkdir debug
mkdir results
rm -rf CMakeCache.txt CMakeFiles
cmake .
make release
make
mpirun -np N ./pfc
```

##### Rerunning program with only parameter changes

```
make runclean
make release
make
mpirun -np N ./pfc
```

## References

## Acknowledgement


