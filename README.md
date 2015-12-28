# lj_potential
Lennard Jones potential and forces for each atom of the system  

Program calculates Lennard-Jones potential and forces for X,Y and Z dimensions acting on each atom of the system. 
The system is described as .xyz file.  
Input:  
- E - depth of the potential well; 
- S (sigma) - finite distance at which the inter-particle potential is zero; 
- name of an input .xyz file; 
- name of an output forces file. 

The format of output file is: 

Fx Fy Fz En  

where Fx, Fy, Fz are forces and En is energy.
