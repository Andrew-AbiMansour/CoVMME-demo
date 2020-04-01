In this file you can describe your simulation and include instructions on how to execute the simulation [script](Script) files. 
E.g.

The simuation consisted of 3 stages: energy minimization / structure optimization, equilibration (NVT), and production (NPT). The [charmm36](Files/Forcefields/par_all36_prot.prm) forcefield was used.
The input files consist of:
- Simulation-specific (NAMD) [files](Files)
- A bash script [file](Scripts/run.sh) for execution

## Structure optimization
Starting from the [prepared structure](../Structural_Data/Final/6lu7.autopsf.pdb) of protease, the system was subjected to energy minmization using the conjugate gradient algorithm for a total of 10,000 steps. Links:
- Optimized final [structure](../Simulation_Output_Files/6lu7.minimized.pdb)
- Trajectory [file](../Simulation_Output_Files/traj_coords.minimized.dcd)

## NVT equilibration


## NPT production run
