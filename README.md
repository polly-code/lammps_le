# Forked LAMMPS
## Introduction
We developed a module for the LAMMPS [0] to simulate loop extrusion.
Loop extrusion is when active ATP-dependent motor (multiprotein complex) cohesin extrudes a loop from DNA [1]. We developed this module to study the effect of cohesin and CTCF on chromosome dynamics [2].
## Mechanics
We simulate extruders (cohesins) as additional sliding links. With the predefined probability (P1) loop can be initiated between {i, i+2} beads and then every predefined number of steps (N1) we will try to shift it for 1 bead each base: [i, j] -> [i-1, j+1] -> [i-2, j+2] and so on. There is a predefined probability of releasing the extruder from the chain (P2). Extruders cannot pass through each other, stalling until one of them is released.
We also have three barriers: left, right, and roadblocks. Left barrier blocks loop coming from the right, left barrier blocks loop coming from the left, and roadblock blocks both.
## Fix description
> fix loop all extrusion 17500 1 2 3 1.0 2 4

Here is standard LAMMPS syntax:
fix - keyword
loop - the name of the fix
all - group of particles to which we apply this fix
extrusion - specific name of the fix
17500 - N1, amount of steps we try to shift links
1 - the type of neutral type monomers from the chain
2 - left barrier
3 - right barrier
1.0 - the probability of going through the block
2 - the bond type of the extruder
4 - the type of roadblocks
To run the code, one needs proper loading and unloading of extruders. There are two more fixes:
> fix loading all ex_load 7000 1 1 1.12 2 prob 0.001 684474 iparam 1 1 jparam 1 1

> fix unloading all ex_unload 7000 2 0.5 prob 0.001 456456

Their syntax is similar to creating and breaking bonds from the MC package from the LAMMPS website.

https://docs.lammps.org/fix_bond_create.html

https://docs.lammps.org/fix_bond_break.html

## References
[0]. Plimpton, Steve. "Fast parallel algorithms for short-range molecular dynamics." Journal of computational physics 117.1 (1995): 1-19.

[1]. Fudenberg, Geoffrey, et al. "Formation of chromosomal domains by loop extrusion." Cell reports 15.9 (2016): 2038-2049.

[2]. Mach, Pia, et al. "Live-cell imaging and physical modeling reveal control of chromosome folding dynamics by cohesin and CTCF." BioRxiv (2022).
