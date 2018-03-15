# Crowding

This simulation prototype allows to simulate a system in the presence of a
certain amount of crowding particles. The diameter of the crowders is set to be
double the diameter of the monomers, binders and tracers. The total volume
fraction of the system is $\rho$. In order for the system to be easily
initializable, and the integration to be feasible, $\rho$ should be less than
30%.

## Invocation
```
python sbs_tracers.py phi e rho init_seed integrate_seed
```

## Parameters
- `phi`: the fraction of binding sites over the total number of monomers in the
  polymer ($\phi$ in the paper)
- `e`: the depth of the Lennard-Jones attractive potential between the tracers
  and the polymer, in $k_B T$ ($\varepsilon$ in the paper)
- `rho`: the volume fraction of the simulation ($\rho$)
- `init_seed`: random number generator seed for initializing the system
- `integrate_seed`: random number generator seed for integration of the
  equations of motion in the system

## Example invocation

```
python sbs_tracers.py 0.02 2.7 0.25 785478 98787463
```
will generate $10^8$ time steps and 10 000 simulation frames stored in the file
`sbs_tracers-phi-0.02-e-2.7-rho-0.25.gsd`.
