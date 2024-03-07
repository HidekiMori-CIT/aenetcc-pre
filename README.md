# Installation

0. This manual explains only clean install case using gnu g++, gfortran, and openmpi under linux OS.

2. Clone this package to your appropriate directory.
``` 
git clone https://github.com/HidekiMori-CIT/aenetcc-pre.git
```

2. Download and unpack LAMMPS package (currently stable_2Aug2023.tar.gz) from Github to same directory of lammps-extend:  
``` 
wget https://github.com/lammps/lammps/archive/refs/tags/stable_2Aug2023.tar.gz
tar -xvzf stable_2Aug2023.tar.gz
``` 

3. Copy USER-ML-AENET-CC/ in lammps-extend/ to lammps-stable_2Aug2023/src/ and complie LAMMPS with aenet-cc module.
```
cp -r ./aenetcc-pre/USER-ML-AENET-CC/ ./lammps-stable_2Aug2023/src/  
cd ./lammps-stable_2Aug2023/src/  
make yes-user-ml-aenet-cc
make mpi
```
_note_: If you use intel compiler, replace mpi to icc_openmpi. Other compile option of LAMMPS, please see LAMMPS manual.
