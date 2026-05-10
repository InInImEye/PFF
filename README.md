# Phase field fracture/ crack program (PFF) 


## About


This is the source code of phase field fracture (or crack) simulation program written in C++ and makes use of [deal.II](https://dealii.org/) finite element library.

This program is written based on the phase field crack models presented in the following papers:
- [A review on phase-field models of brittle fracture and a new fast hybrid formulation](https://doi.org/10.1007/s00466-014-1109-y)
- [Phase-Field Modeling of Ductile Fracture](https://doi.org/10.1007/s00466-015-1151-4)

Dealii 9.5.1 [ref](https://doi.org/10.1515/jnma-2023-0089) [ref2](https://github.com/dealii/dealii/tree/dealii-9.5)

## Program usage

## File structure

## Installation of Libraries

Several libraries are required to run this code. These include P4est, PetSc, Parmetis, and Deal.ii. Deal.ii makes it possible to install all these with the already made available [candi](https://github.com/dealii/candi) script, which we can use to install the libraries and automatically link to each other. As of February 2025, this code is compatible with deal.ii version 9.5.1 **(Compatibility with the latest version has not been checked, the code might require forward porting with minor changes or possibly no changes at all.)**. All other libraries versions are automatically selected by candi during installation. We here describe the installation procedure with the Ubuntu OS.

First check the repositories are accessible with the `apt` and they are not giving any error. Try:
```
sudo apt update
```
If no error shows up then you are good to proceed. Next make sure the requirements following are met:
1. Use gcc version > 5. Check the version using
    ```
    gcc -v
    ```
2. Have MPICH installed with the gcc > 5. Check using
    ```
    which mpicc
    ```
    or
    ```
    mpirun --version
    ```
### Install MPICH
1. Download the latest version of MPICH from its [website](https://www.mpich.org/downloads/).
2. Install by following the [Guide](https://www.mpich.org/documentation/guides/). At the time of this code development Version [4.1.2](https://www.mpich.org/static/downloads/4.1.2/mpich-4.1.2-installguide.pdf) was available. Installation commands for `bash` terminal are summarized below:

    ```
    mkdir ~/myloc                         # Make a directory where you want to extract files
    
    tar xfz mpich-4.1.2.tar.gz -C ~/myloc # Unpack the tar file to the location.

    mkdir ~/mpich-install                 # Choose an installation directory 

    mkdir /tmp/you/mpich-4.1.2            # Choose a build directory

    # Configure MPICH, specifying the installation directory, and running
    the configure script in the source directory

    cd /tmp/you/mpich-4.1.2
    ~/myloc/mpich-4.1.2/configure -prefix=~/mpich-install 2>&1 | tee c.txt

    make 2>&1 | tee m.txt                 # Build MPICH

    make install 2>&1 | tee mi.txt        # Install the MPICH 

    export PATH="$HOME/myloc/mpich-4.1.2/bin:$PATH"  # Set PATH by adding this to end of ~/.bashrc file

    # Restart the terminal and check that the MPICH is accessible by running:
    which mpicc
    which mpiexec
    ```

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
