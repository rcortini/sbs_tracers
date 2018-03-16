# sbs_tracers
---

This package contains the Python scripts to run and analyze the simulations that
were presented in the paper

R. Cortini and G. Filion ["Theoretical principles of transcription factor
traffic on folded chromatin"](https://doi.org/10.1101/164541).

# Table of contents
1. [Introduction](#intro)
2. [Environment setup](#envsetup)
    1. [Requirements](#requirements)
    2. [Loading the Python module](#pythonmodule)
3. [Running a simulation](#running)
4. [Analyzing a simulation](#analyzing)

## Introduction

The code presented here is organized in two main subdirectories: `sims` and
`analysis`.

The **sims** folder contains several subdirectories:

1. `monovalent`: code for simulating monovalent tracers (Supplementary Note 4)
2. `multivalent`: code for simulating multivalent tracers (Figures 2, 3, 4, 5 of
   Main Text, Supplementary Notes 1, 2, 3, 7)
3. `crowding`: code for simulating tracers in the presence of a crowded
   environment (Supplementary Note 6)
4. `r_sigma`: code for simulating tracers with varying diameter, and comes in
   two flavors
    1. `monovalent`: the tracers are monovalent
    2. `multivalent`: the tracers are multivalent
    Both are in Supplementary Note 5.

The **analysis** folder contains two subdirectories:

1. `sbs_tracers_analysis`: a Python module that contains the core elements to
   run and analyze a simulation
2. `scripts`: a few scripts that perform a complete analysis of a simulation
   trajectory

## Environment setup

First let's set up the environment to run and analyze the simulations. First of
all, clone the repository to your local machine. From a shell:
```
cd /path/to/download
git clone https://github.com/rcortini/sbs_tracers
```

### Requirements

The software prerequisites for the package to work are the following

- [HOOMD-blue v2.1+](http://hoomd-blue.readthedocs.io/en/stable/)
- [MDAnalysis v0.17+](https://www.mdanalysis.org/)
- [GSD](https://bitbucket.org/glotzer/gsd)
- [NumPy](http://www.numpy.org/)
- [SciPy](https://www.scipy.org/)

Once the software prerequisites are correctly installed, you should be able to
run correctly an `import` from Python of all the following: `hoomd`,
`MDAnalysis`, `gsd`, `numpy`, `scipy`. If all the `import` commands succeed, you
are ready to go.

### Loading the Python module

The next step is to manually add the `sbs_tracers_analysis` subfolder to your
Python environment. In UNIX-like environment, you set your Python path by
setting the environment variable `PYTHONPATH`. If you create a symbolic link to
that directory, you should be able to import the module. For example, in my case

```
ln -s /home/rcortini/soft/sbs_tracers/analysis/sbs_tracers_analysis /home/rcortini/soft/python
```
That is because I cloned the git repository in `/home/rcortini/soft/sbs_tracers`
and my `PYTHONPATH` contains the folder `/home/rcortini/soft/python`. If you do
that, you should be able to do `import sbs_tracers_analysis` from a Python
shell. You are then ready to run the simulations and analyze them.

## Running a simulation
The `sims` directory contains several subdirectories, each of which contains an
`sbs_tracers.py` file. From a shell, for example you could run
```
cd /path/to/download/sims/multivalent
python sbs_tracers.py <options>
```
which will run a simulation. The output will be a `.gsd` file, which contains
both the topology and the trajectory of the simulated particles. The output file
name will contain also the parameters passed to the script through the
`<options>`.  See the individual README files for the options to each of the
simulation files.

All the simulation files contain the options `init_seed` and `integrate_seed`,
which are the seeds for the random number generator for initialization of the
system and integration of the equations of motion, respectively. Having two
separate seeds is useful if you want to check different trajectories starting
from the same initial conditions. The seeds are not written to the trajectory
file names.

## Analyzing a simulation
In the `analysis/scripts` subdirectory the scripts needed to analyze
the simulation trajectories are located. In general, you invoke those scripts
from a shell like follows:
```
python script.py trajectory.gsd <options>
```
See the README in that folder for the description of each of the analysis
scripts.
