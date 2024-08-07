# Update
2024/07/23: Add Fe05_zbl.   
2024/07/23: Rename pair_aenetcc.* pair_aenet_cc.* and add weight_functions.h.  
2024/03/16: Update pair_aenetcc.* and aenet_sfb.h.  

# Installation
0. This manual explains only clean install case using gnu g++, gfortran, and openmpi under linux OS.

2. Clone this package to your appropriate directory.
``` 
git clone https://github.com/HidekiMori-CIT/aenetcc-pre.git
```

2. Download and unpack LAMMPS package (currently stable_2Aug2023.tar.gz) from Github to same directory of aenetcc-pre:  
``` 
wget https://github.com/lammps/lammps/archive/refs/tags/stable_2Aug2023.tar.gz
tar -xvzf stable_2Aug2023.tar.gz
``` 

3. Copy USER-ML-AENET-CC/ in aenetcc-pre/ to lammps-stable_2Aug2023/src/ and complie LAMMPS with aenet-cc module.
```
cp -r ./aenetcc-pre/USER-ML-AENET-CC/ ./lammps-stable_2Aug2023/src/  
cd ./lammps-stable_2Aug2023/src/  
make yes-user-ml-aenet-cc
make mpi
```
_note_: If you use intel compiler, replace mpi to icc_openmpi. Other compile option of LAMMPS, please see LAMMPS manual.

# Citing of this package and ANN potential
[1] H. Mori and T. Ozaki, Phys. Rev. Mater., **4**, 040601 (2020).  
[2] H. Mori, _et. al._, Phys. Rev. Matter., **7**, 063605 (2023).
