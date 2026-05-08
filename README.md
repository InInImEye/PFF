# Phase field fracture/ crack program (PFF) 


## About


This is the source code of phase field fracture (or crack) simulation program written in C++ and makes use of [deal.II](https://dealii.org/) finite element library.

This program is written based on the phase field crack models presented in the following papers:
- [A review on phase-field models of brittle fracture and a new fast hybrid formulation](https://doi.org/10.1007/s00466-014-1109-y)
- [Phase-Field Modeling of Ductile Fracture](https://doi.org/10.1007/s00466-015-1151-4)

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

Dr. Sabyasachi Chatterjee [Webpage](https://web.iitd.ac.in/~sabyasachi/) [Personal Webpage](https://sites.google.com/view/sabyasachichatterjee/home) [Researchgate](https://www.researchgate.net/profile/Sabyasachi-Chatterjee-2) 

Dr. Anup Basak [Webpage](https://old.iittp.ac.in/dr-anup-basak) [Google scholar](https://scholar.google.com/citations?user=m_TDGD8AAAAJ&hl=en)
