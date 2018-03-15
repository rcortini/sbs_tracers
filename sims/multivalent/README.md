# Multivalent tracers

The main results of the paper were obtained by running this simulation script.

## Invocation
```
python sbs_tracers.py phi e n init_seed integrate_seed
```

## Parameters
- `phi`: the fraction of binding sites over the total number of monomers in the
  polymer ($\phi$ in the paper)
- `e`: the depth of the Lennard-Jones attractive potential between the tracers
  and the polymer, in $k_B T$ ($\varepsilon$ in the paper)
- `n`: the simulation replicate number
- `init_seed`: random number generator seed for initializing the system
- `integrate_seed`: random number generator seed for integration of the
  equations of motion in the system

## Example invocation

```
python sbs_tracers.py 0.02 2.7 8 76435 671438 
```
will generate $10^8$ time steps and 10 000 simulation frames stored in the file
`sbs_tracers-phi-0.02-e-2.7-8.gsd`.

## Number of tracers
The number of tracers in the simulation is set in the script through the
`ntracers` variable. The main results of the paper were obtained using 10
tracers, but results for 200 tracers are also reported.
