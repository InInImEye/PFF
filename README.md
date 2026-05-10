# Phase field fracture/ crack program (PFF) 


## About


This is the source code of phase field fracture (or crack) simulation program written in C++ and makes use of [deal.II](https://dealii.org/) finite element library.

This program is written based on the phase field crack models presented in the following papers:
- [A review on phase-field models of brittle fracture and a new fast hybrid formulation](https://doi.org/10.1007/s00466-014-1109-y)
- [Phase-Field Modeling of Ductile Fracture](https://doi.org/10.1007/s00466-015-1151-4)

Dealii 9.5.1 [ref](https://doi.org/10.1515/jnma-2023-0089) [ref2](https://github.com/dealii/dealii/tree/dealii-9.5)

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

## Acknowledgement

Dr. Sabyasachi Chatterjee [1](https://web.iitd.ac.in/~sabyasachi/ "Webpage") [2](https://sites.google.com/view/sabyasachichatterjee/home "Personal Webpage") [3](https://www.researchgate.net/profile/Sabyasachi-Chatterjee-2 "Researchgate") 

Dr. Anup Basak [1](https://old.iittp.ac.in/dr-anup-basak "Webpage") [2](https://scholar.google.com/citations?user=m_TDGD8AAAAJ&hl=en "Google scholar")
