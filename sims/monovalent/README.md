# Monovalent tracers

The method for simulating monovalent tracers is akin to the one presented by
Brackley et al. ["Facilitated Diffusion on Mobile DNA: Configurational Traps and
Sequence Heterogeneity"](https://doi.org/10.1103/PhysRevLett.109.168103). That
is to say, the tracers are now *composite* particles, made of a central particle
and another, smaller particle which is rigidly attached to the main one. The
integration of the system is now based on the rigid body integrator of
HOOMD-blue. Everything else is equal to the case of the multivalent tracers.

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
