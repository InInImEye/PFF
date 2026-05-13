# Phase field fracture/ crack program (PFF) 


## About


This is the source code of phase field fracture (or crack) simulation program written in C++ and makes use of [deal.II](https://dealii.org/) finite element library.

This program is written based on the phase field crack models presented in the following papers:
- [A review on phase-field models of brittle fracture and a new fast hybrid formulation](https://doi.org/10.1007/s00466-014-1109-y)
- [Phase-Field Modeling of Ductile Fracture](https://doi.org/10.1007/s00466-015-1151-4)

deal.II 9.5.1 [ref](https://doi.org/10.1515/jnma-2023-0089) [ref2](https://github.com/dealii/dealii/tree/dealii-9.5)

## Program usage

## File structure

## Installation of Libraries

Several libraries are required to run this code. These include P4est, PetSc, Parmetis, and Deal.ii. Deal.ii makes it possible to install all these with the already made available [candi](https://github.com/dealii/candi) script, which we can use to install the libraries and automatically link to each other. As of February 2025, this code is compatible with deal.II version 9.5.1 **(Compatibility with the latest version has not been checked, the code might require forward porting with minor changes or possibly no changes at all.)**. All other libraries versions are automatically selected by candi during installation. We here describe the installation procedure with the Ubuntu OS.

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
    Follow the step below if not installed.

### Install MPICH

1. Before installing MPICH also ensure that the other dependency for deal.II is fulfilled. For the Ubuntu OS you can check [here](https://github.com/dealii/candi/blob/master/deal.II-toolchain/platforms/supported/ubuntu.platform). Later, the same prompt will show up at the start of deal.II installation procedure which you can safely proceed without interrupting the installation. For other OS the requirement is mentioned in the [same](https://github.com/dealii/candi/blob/master/deal.II-toolchain/platforms/supported/) folder.
2. Download the latest version of MPICH from its [website](https://www.mpich.org/downloads/).
3. Install by following the [Guide](https://www.mpich.org/documentation/guides/). At the time of this code development Version [4.1.2](https://www.mpich.org/static/downloads/4.1.2/mpich-4.1.2-installguide.pdf) was available. Installation commands for `bash` terminal are summarized below:

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

    # Set PATH by adding this to end of ~/.bashrc file
    export PATH="$HOME/mpich-install/bin:$PATH"  

    # Restart the terminal and check that the MPICH is accessible by running:
    which mpicc
    which mpiexec
    mpirun --version
    ```

### Install deal.II and other libraries

Now, these pre-requisites are fulfilled according to the previous sections:

1. Installation of the dependencies, according to your supported OS as mentioned in: `candi/deal.II-toolchain/platforms/supported/`
2. Installation of MPICH using the instructions.

Next, we install the rest of the libraries using candi.

1. Enter the commands for creating a download and installation directory, and download the source code as given below. Do this in the directory `your_location` where you want your deal.II and related libraries to be installed.
    ```
    mkdir CL
    cd CL
    mkdir FEM
    git clone https://github.com/dealii/candi
    cd candi
    ```

2. Open the file `candi.cfg` to configure installed libraries using text editor of your choice, we use `vim` here.
    ```
    vim candi.cfg
    ```

3. Edit the file according to the packages and configuration required and save it.
   We are mainly going to be using the following packages, so uncomment them:
    ```
    PACKAGES="${PACKAGES} once:parmetis
    PACKAGES="${PACKAGES} once:hdf5
    PACKAGES="${PACKAGES} once:p4est
    PACKAGES="${PACKAGES} once:petsc
    PACKAGES="${PACKAGES} dealii
    ```
4. In the same file modify the version to install the correct deal.II version.
    ```
    DEAL_II_VERSION=v9.5.1
    ```

5. Install:
    
    ```
    ./candi.sh --prefix=/your_location/CL/FEM --platform=./deal.II-toolchain/platforms/supported/YOUR_OS.platform
    ```
    **Make sure `YOUR_OS` is in the supported platforms list.** For our system we used `ubuntu.platform`

6.  Once the setup is complete you can proceed to run simulations in the `your_location/CL/FEM/deal.II-v9.5.1` directory. If you want to use it outside the installation directory you need to add a `PATH` to it, similar to the MPICH installation step.
    ```
    export DEAL_II_DIR=/your_location/CL/FEM/deal.II-v9.5.1/
    ```

## Running the program

We have mentioned the steps for compiling and running the program under same section, since it is pretty straightforward.

##### First run of the program

Some directories need to be made before running the program since program looks for it, to output various information.

```
mkdir debug
mkdir results
rm -rf CMakeCache.txt CMakeFiles
cmake .
make release
make
mpirun -np N ./pfc
```

##### Rerunning the program

`./update.sh`

copy the results and debug and params files before rerunning

if you make changes only to the the params.in then make not required

if you make changes to the program then make 


```
make runclean
make release
make
mpirun -np N ./pfc
```

## Acknowledgement

This research was funded by [SERB](https://serb.gov.in/) (currently ANRF) via project No. RP04378G. Gratitude to [IRD](https://ird.iitd.ac.in/) [IITD](https://home.iitd.ac.in/) for the opportunity to work on the project. Dr. Sabyasachi Chatterjee ([1](https://web.iitd.ac.in/~sabyasachi/ "Webpage") [2](https://sites.google.com/view/sabyasachichatterjee/home "Personal Webpage") [3](https://www.researchgate.net/profile/Sabyasachi-Chatterjee-2 "Researchgate")) provided continuous guidance throughout the project through insightful discussions, relevant literature, and all necessary support. The discussions with Dr. Anup Basak ([1](https://old.iittp.ac.in/dr-anup-basak "Webpage") [2](https://scholar.google.com/citations?user=m_TDGD8AAAAJ&hl=en "Google scholar")), Mechanical Engineering, IIT Tirupati, on implementation of the model is greatly aknowledged.
