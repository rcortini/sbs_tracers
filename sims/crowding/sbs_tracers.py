from __future__ import print_function
import os
import numpy as np
import hoomd, hoomd.deprecated, hoomd.md
import sbs_tracers_analysis as sbs
import sys

# check for proper invocation
if len(sys.argv) < 6 :
    print("python sbs_tracers.py <phi> <e> <n> <init_seed> <integrate_seed>")
    sys.exit(1)

# initialize hoomd-blue
hoomd.context.initialize()

###############################
# parameters to vary
###############################
# concentration of binding sites
phi = float(sys.argv[1])
# tracer-to-polymer affinity (kT)
e = float(sys.argv[2])
# simulation replicate number
rho = float(sys.argv[3])
# random number generator seeds
init_seed = int(sys.argv[4])
integrate_seed = int(sys.argv[5])

###############################
# physical parameters
###############################
kT = 1.0                 # temperature in units of energy
N = 1024                 # number of beads in the polymer
nbs = int(N*phi)         # number of binding sites
sep = 0.8                # minimum separation of beads at init
b_polymer = 1.7          # more than twice the separation
polymer_r0 = 1.2         # polymer bond length, local units
ntracers = 200           # number of tracer molecules
k = 330.0                # polymer stiffness constant
sigma = 1.0              # particle diameter, local units
sigmac = 2.0             # crowder diameter, local units
e_repulsion = 1.0        # repulsion LJ epsilon
e_attraction = 10.0      # binder-to-binding sites affinity (kT)
r_cut = 3.0*sigmac       # LJ interaction cutoff distance

###############################
# general parameters
###############################
name = "sbs_tracers-phi-%.2f-e-%.1f-rho-%.1f"%(phi,e,rho)
gsd_freq = 10000
run_gsd = "%s.gsd"%name

###############################
# simulation parameters
###############################
dt = 0.005
run_length = 1e8

###############################
# simulation box init
###############################
L = 50.0                        # local dimensions of the simulation box
nfactors = nbs                  # there are always an equal number of binders and binding sites
box = hoomd.data.boxdim(L=L)

# number of crowders to obtain assigned rho value
ncrowders = 6*L**3*rho/(np.pi*sigmac**3) - (1+2*phi)*sigma**3*N/sigmac**3

###############################
# polymer init
###############################
polymer_type = ['p']*N
np.random.seed(seed=init_seed)
bs = np.random.choice (N,size=nbs)
for i in bs :
    polymer_type[int(i)] = 'a'
polymer1 = dict(bond_len=b_polymer, type=polymer_type,
                bond="linear", count=1)

###############################
# init diffusing beads
###############################
factors = dict(bond_len=b_polymer, type=['b'],
                bond="linear", count=nfactors)
tracers = dict(bond_len=b_polymer, type=['t'],
                bond="linear", count=ntracers)
crowders = dict(bond_len=b_polymer, type=['c'],
                bond="linear", count=ncrowders)

###############################
# system init
###############################
# check if initial file exists
if os.path.exists (run_gsd) :
    # this is the case in which we need to load the system starting from the
    # trajectory file
    import gsd.hoomd
    t = gsd.hoomd.open(name=run_gsd,mode='rb')
    fn = len(t)-1
    system = hoomd.init.read_gsd(filename=run_gsd,frame=fn)
else :
    # init the system
    system = hoomd.deprecated.init.create_random_polymers(box=box,
                                polymers=[polymer1,factors,tracers,crowders],
                                separation=dict(p=sep,
                                                a=sep,
                                                b=sep,
                                                t=sep,
                                                c=sep), seed=init_seed)
    for particle in system.particles :
       if particle.type == 'c' :
           particle.diameter = sigmac

# polymer bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=k, r0=polymer_r0)

# LJ interactions
nlist = hoomd.md.nlist.cell ()
lj = hoomd.md.pair.lj(r_cut,nlist)

# interactions of the p particles
lj.pair_coeff.set('p', 'p', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 'a', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 'b', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)

# interactions of the a particles
lj.pair_coeff.set('a', 'a', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('a', 'b', epsilon=e_attraction, sigma=sigma, alpha=1.0)
lj.pair_coeff.set('a', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)

# interactions of the b particles
lj.pair_coeff.set('b', 'b', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('b', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)

# interactions of the t particles
lj.pair_coeff.set('t', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)

# crowders interact with all particles with the same repulsive potential
sigmaLJ_crowders = 0.5 * (sigma + sigmac)
for t in ['p','a','b','t','c'] :
    lj.pair_coeff.set('c', t, epsilon=e_repulsion, sigma=sigmac, alpha=0.0)

# trajectory output
all = hoomd.group.all()
hoomd.dump.gsd(filename=run_gsd, period=gsd_freq, 
               group=hoomd.group.all(),phase=0)

###############################
# RUN!
###############################

# integrate NVT for a bunch of time steps
hoomd.md.integrate.mode_standard(dt=dt)
hoomd.md.integrate.langevin(group=all, kT=kT, seed=integrate_seed)
hoomd.run_upto(run_length)
