# Analysis scripts

This subdirectory contains scripts that perform the main analysis routines that
were used to obtain the results of the paper. All the outputs of the scripts are
in `.npy` format, which is the standard form for storing `numpy` arrays.

The main argument of the analysis scripts is the trajectory file, which
typically has `.gsd` extension.

All the analysis scripts accept the following two parameters:
- `teq`: the "equilibration time", that is the number of frames of the
  trajectory at the beginning of the simulation that should be excluded from the
  analysis
- `tsample`: the "sampling time", that is the interval between the frames that
  is used to perform the calculations

This will result in a number of frames `nframes` which is equal to
`(tot_nframes-teq)/tsample`, where `tot_nframes` is the total number of
simulation frames in the stored trajectory.

## tracers_analysis.py
This scripts performs a full analysis of the tracers' binding patterns to the
polymer, and therefore constitutes the core of the paper's results.

Outputs several files

- `DKL_t.npy`: an array of `nframes` containing the value of the
  Kullback-Leibler divergence between the tracer traffic and the polymer
  contacts at each time step of the simulation.
- `polymer_contacts.npy`: an `N`x`N` matrix containing the total number of
  contacts between each monomer and each other monomer, summed over all the time
  steps. Is akin to the Hi-C matrix.
- `traffic.npy`: a `N`-sized vector which contains the total number of contacts
  that the tracers made with each monomer, summed over tracers and summed over
  time frames.
- `coverage.npy`: `ntracers`-sized vector which contains the percentage of
  monomer sites that each of the tracers of the simulation visited during the
  simulation.
- `ps.npy`: is the average probability of contact between monomers as a function
  of the linear separation of the monomers.
- `r.npy`: Pearson correlation coefficient between the polymer contacts and the
  traffic of the tracers.

### Invocation

```
python tracers_analysis.py trajectory [options...]
```

### Options

See above for the `teq` and `tsample` options.

- `p_threshold`: the threshold to determine whether two polymer monomers are in
  contact or not.
- `t_threshold`: the threshold to determine whether a tracer is in contact with
  the polymer or not.
- `threshold`: supply this option if you want the thresholds above to be equal.
