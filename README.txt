# MPI Version of ANGULA

This repository provides an MPI-parallelized version of the ANGULA package.

ANGULA was originally developed by Dr. Luis Carlos Pardo  
(Group of Characterization of Materials, Polytechnic University of Catalonia).

Original project:
https://gcm.upc.edu/en/members/luis-carlos/angula

---

## Motivation

The original serial implementation of ANGULA can become computationally expensive for large systems or long trajectories.  

This MPI version enables parallel processing of configurations, significantly reducing execution time depending on:
- the number of available processors
- the system size

The implementation was tested on an AMD Ryzen 7 processor using 16 MPI processes. OS: Linux Ubuntu 24.04
Note: It was not tested on other OS windows or mac or other linux distos.
---

## Integration with ANGULA

To use this version:

1. Place the MPI source file inside the ANGULA source directory:

2. Ensure it uses the same modules as the original `angula.for`.

---

## Compilation

Compile using an MPI-enabled Fortran compiler:

```bash
mpif90 -ffixed-line-length-none -std=legacy -fallow-argument-mismatch angula_mpi_optionB.for -o angula_mpi.exe


Note: Compilation may produce warnings due to legacy Fortran structures, but the executable is expected to work correctly.




## Configuration Limits

The maximum system size is controlled by the following parameters, which must be updated consistently across:

ANGULA source files
Traductor source files

Parameters:

dimat → maximum number of atoms per molecule
dimmol → maximum number of molecules per molecule type
diEat → maximum number of atomic energy terms associated with dimat
## Performance

Benchmark results (cellulose system, 900 configurations):

- Serial execution: ~6 days 15 hours  
- MPI (8 processes, configuration-level parallelism): ~20 hours  
- MPI (8 processes, full parallel per configuration): ~16 hours  which this code upgrade uses.

Significant speedup achieved through parallel execution.

## Notes:

1- This implementation focuses on improving performance through MPI parallelization.
2- The code is based on legacy Fortran and may generate compiler warnings.
3- Contributions and improvements are welcome.

## Author

Rashed M. R. Aldughmi

## Contributions

- MPI parallelization of ANGULA
- Performance optimization for large systems
