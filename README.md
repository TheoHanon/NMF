

---

# Simulating Convection

## Description

This project contains all the code used to write the report of the Homework assigment. The simulation are conducted in `C` code, while all the post processing tasks are done with `Python`.  

## Directory Structure

The project is organized into several directories:

- `src/`: Contains the C source files for the simulation.
- `header/`: Contains the header files with function declarations and structure definitions.
- `post_processing/`: Contains Python scripts for processing and visualizing the simulation results. The `data` subdirectory within it is used to store output data files.
- `makefile`: The Makefile used to compile the source code, run the simulation with specified parameters, and initiate data processing and visualization.

## Prerequisites

- GCC compiler for compiling C code.
- Python with necessary libraries for running post-processing scripts.

## Compilation

To compile the source code into an executable, use the following command:

```bash
make
```

This will compile the `src/hw1.c` file (along with any other C source files in the `src/` directory) into an executable named `hw1`, using the GCC compiler with the `-Wall` and `-Wextra` flags for additional warnings.

## Running Simulations

To run a simulation with specific parameters, use the `make run` command with the following options:

- `N`: The grid size.
- `SCHEME_TYPE`: The finite difference scheme type (`E2`, `E4`, `I4`, `ED`).
- `GRID_TYPE`: The type of grid to use (`UNIF` for uniform, `NONUNIF` for non-uniform).
- `INITIAL_TYPE`: The initial condition (`GAUSSIAN`, `WAVEPACKET`).

For example:

```bash
make run N=128 SCHEME_TYPE=ED GRID_TYPE=NONUNIF INITIAL_TYPE=GAUSSIAN
```

This will run the `hw1` executable with the specified parameters and generate output data files in the `post_processing/data` directory. In addition it will display an animation. Note that if you want to tune any other parameters (e.g. L, U, etc.) of the simulation you should go to `src\hw1.c` and change them directly there. 

## Data Processing and Visualization

You can visualize the plots of the howemork by simply using the following command:

```bash
make plot
```

This will run the `PARTI.py`, `PARTII.py`, and `PARTIII.py` scripts. Note that these scripts should be launched only through this command otherwise you will have issues. This command is self-sufficient it will run the necessary simulation to create the plots. 

## Cleaning Up

To clean up the generated executable and data files, use the following command:

```bash
make clean
```

This will remove the `hw1` executable and all `.txt` data files in the `post_processing/data` directory.

## Help

For usage instructions, you can use the `make help` command:

```bash
make help
```

This will display information about how to run simulations with different parameters and options available.

---
