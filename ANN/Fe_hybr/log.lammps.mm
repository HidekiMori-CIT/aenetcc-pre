LAMMPS (2 Aug 2023)
atom_style      atomic
units           metal
boundary        p p p

lattice         bcc 2.9
Lattice spacing in x,y,z = 2.9 2.9 2.9
region          region1 block 0 1 0 1 0 1 units lattice
create_box      1 region1
Created orthogonal box = (0 0 0) to (2.9 2.9 2.9)
  1 by 1 by 1 MPI processor grid

create_atoms    1 box
Created 2 atoms
  using lattice units in orthogonal box = (0 0 0) to (2.9 2.9 2.9)
  create_atoms CPU = 0.000 seconds

pair_style      hybrid/scaled 0.9 aenet/cc 0.1 aenet/cc
pair_coeff      * * aenet/cc 1 v01 Fe 10tw-10tw.ann Fe
pair_coeff      * * aenet/cc 2 v03 Fe 10sw-10sw.nn Fe
mass            1 55.845

thermo_style    custom step pe lx vol press
thermo          10

fix             f1 all box/relax iso 0.0

minimize        1.0e-20 0.0 1000 100000
Neighbor list info ...
  update: every = 1 steps, delay = 0 steps, check = yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 8.5
  ghost atom cutoff = 8.5
  binsize = 4.25, bins = 1 1 1
  2 neighbor lists, perpetual/occasional/extra = 2 0 0
  (1) pair aenet/cc, perpetual
      attributes: full, newton on
      pair build: full/bin/atomonly
      stencil: full/bin/3d
      bin: standard
  (2) pair aenet/cc, perpetual, trim from (1)
      attributes: full, newton on, cut 7.5
      pair build: trim
      stencil: none
      bin: none
WARNING: Energy due to 1 extra global DOFs will be included in minimizer energies
 (../min.cpp:225)
Per MPI rank memory allocation (min/avg/max) = 4.747 | 4.747 | 4.747 Mbytes
   Step         PotEng           Lx           Volume         Press     
         0  -8959.71        2.9            24.389        -105699.76    
        10  -8959.7148      2.8971         24.315906     -103052.88    
        20  -8959.7194      2.8942         24.242958     -100122.82    
        30  -8959.7239      2.8913         24.170157     -96931.291    
        40  -8959.7282      2.8884         24.097501     -93497.551    
        50  -8959.7323      2.8855         24.024991     -89838.412    
        60  -8959.7363      2.8826         23.952627     -85968.29     
        70  -8959.7401      2.8797         23.880408     -81899.341    
        80  -8959.7437      2.8768         23.808334     -77641.668    
        90  -8959.7471      2.8739         23.736406     -73203.578    
       100  -8959.7503      2.871          23.664622     -68591.877    
       110  -8959.7532      2.8681         23.592984     -63812.19     
       120  -8959.7559      2.8652         23.52149      -58869.275    
       130  -8959.7585      2.8623         23.450141     -53767.322    
       140  -8959.7607      2.8594         23.378936     -48510.217    
       150  -8959.7628      2.8565         23.307875     -43101.752    
       160  -8959.7645      2.8536         23.236959     -37545.783    
       170  -8959.7661      2.8507         23.166186     -31846.319    
       180  -8959.7674      2.8478         23.095558     -26007.575    
       190  -8959.7684      2.8449         23.025073     -20033.982    
       200  -8959.7691      2.842          22.954732     -13930.209    
       210  -8959.7696      2.8391         22.884534     -7701.2172    
       220  -8959.7698      2.8362         22.814479     -1352.3925    
       224  -8959.7698      2.835589       22.799737     -4.8988385e-07
Loop time of 0.956395 on 1 procs for 224 steps with 2 atoms

99.5% CPU use with 1 MPI tasks x no OpenMP threads

Minimization stats:
  Stopping criterion = energy tolerance
  Energy initial, next-to-last, final = 
     -8959.71000752341  -8959.76979676124  -8959.76979676124
  Force two-norm initial, final = 4.8270177 2.1388908e-11
  Force max component initial, final = 4.8270177 2.1388906e-11
  Final line search alpha, max atom move = 1 2.1388906e-11
  Iterations, force evaluations = 224 226

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.9517     | 0.9517     | 0.9517     |   0.0 | 99.51
Neigh   | 0          | 0          | 0          |   0.0 |  0.00
Comm    | 0.0010101  | 0.0010101  | 0.0010101  |   0.0 |  0.11
Output  | 0.00061606 | 0.00061606 | 0.00061606 |   0.0 |  0.06
Modify  | 0          | 0          | 0          |   0.0 |  0.00
Other   |            | 0.00307    |            |       |  0.32

Nlocal:              2 ave           2 max           2 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:            557 ave         557 max         557 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:              0 ave           0 max           0 min
Histogram: 1 0 0 0 0 0 0 0 0 0
FullNghs:          360 ave         360 max         360 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 360
Ave neighs/atom = 180
Neighbor list builds = 0
Dangerous builds = 0

Total wall time: 0:00:00
