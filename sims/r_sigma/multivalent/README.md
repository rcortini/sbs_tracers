# r_sigma multivalent

This simulation script allows to simulate tracers with varying diameter. The
tracers are multivalent, so only a special precaution needs to be taken when
assigning the cutoff length of the Lennard-Jones interactions.

## Invocation
```
python sbs_tracers.py phi e sigma init_seed integrate_seed
```

## Parameters
- `phi`: the fraction of binding sites over the total number of monomers in the
  polymer ($\phi$ in the paper)
- `e`: the depth of the Lennard-Jones attractive potential between the tracers
  and the polymer, in $k_B T$ ($\varepsilon$ in the paper)
- `sigma`: the desired diameter of the tracers ($\sigma_t$ in the paper)
- `init_seed`: random number generator seed for initializing the system
- `integrate_seed`: random number generator seed for integration of the
  equations of motion in the system

## Example invocation

```
python sbs_tracers.py 0.02 2.7 0.8 7665 9890423
```
will generate $10^8$ time steps and 10 000 simulation frames stored in the file
`sbs_tracers-phi-0.02-e-2.7-sigma-0.8.gsd`.
