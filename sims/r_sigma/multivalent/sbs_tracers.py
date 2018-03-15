from __future__ import print_function
import os
import numpy as np
import hoomd, hoomd.deprecated, hoomd.md
import sbs_tracers_analysis as sbs
import sys

# check for proper invocation
if len(sys.argv) < 6 :
    print("python sbs_tracers.py <phi> <e> <sigma> <init_seed> <integrate_seed>")
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
sigmat = float(sys.argv[3])
# random number generator seeds
init_seed = int(sys.argv[4])
integrate_seed = int(sys.argv[5])

###############################
# physical parameters
###############################
kT = 1.0                 # temperature in units of energy
N = 1024                 # number of beads in the polymer
nbs = int(N*phi)         # number of binding sites
sep = sigmat/2. - 0.1    # minimum separation of beads at init
b_polymer = 1.7          # polymer bond length, local units
polymer_r0 = 1.2         # polymer bond length, local units
ntracers = 10            # number of tracer molecules
k = 330.0                # polymer stiffness constant
sigma = 1.0              # particle diameter, local units
e_repulsion = 1.0        # repulsion LJ epsilon
e_attraction = 10.0      # binder-to-binding sites affinity (kT)
r_cut = 2.5*sigmat       # LJ interaction cutoff distance

###############################
# general parameters
###############################
name = "sbs_tracers-phi-%.2f-e-%.1f-sigma-%.1f"%(phi,e,sigmat)
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
tracers = dict(bond_len=sigmat, type=['t'],
                bond="linear", count=ntracers)

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
                                polymers=[polymer1,factors,tracers],
                                separation=dict(p=sep,
                                                a=sep,
                                                b=sep,
                                                t=sep), seed=init_seed)

# polymer bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=k, r0=b_polymer)

# neighbor list init
nlist = hoomd.md.nlist.cell ()

def LJsigma(s1,s2):
    return 0.5*(s1+s2)

# LJ interactions init
lj = hoomd.md.pair.lj(r_cut,nlist)

# interactions of the p particles
lj.pair_coeff.set('p', 'p', epsilon=e_repulsion, sigma=LJsigma(sigma,sigma), alpha=0.0)
lj.pair_coeff.set('p', 'a', epsilon=e_repulsion, sigma=LJsigma(sigma,sigma), alpha=0.0)
lj.pair_coeff.set('p', 'b', epsilon=e_repulsion, sigma=LJsigma(sigma,sigma), alpha=0.0)
lj.pair_coeff.set('p', 't', epsilon=e_repulsion, sigma=LJsigma(sigma,sigmat), alpha=0.0)

# interactions of the a particles
lj.pair_coeff.set('a', 'a', epsilon=e_repulsion, sigma=LJsigma(sigma,sigma), alpha=0.0)
lj.pair_coeff.set('a', 'b', epsilon=e_attraction, sigma=LJsigma(sigma,sigma), alpha=1.0)
lj.pair_coeff.set('a', 't', epsilon=e_repulsion, sigma=LJsigma(sigma,sigmat), alpha=0.0)

# interactions of the b particles
lj.pair_coeff.set('b', 'b', epsilon=e_repulsion, sigma=LJsigma(sigma,sigma), alpha=0.0)
lj.pair_coeff.set('b', 't', epsilon=e_repulsion, sigma=LJsigma(sigma,sigmat), alpha=0.0)

# interactions of the t particles
lj.pair_coeff.set('t', 't', epsilon=e_repulsion, sigma=LJsigma(sigmat,sigmat), alpha=0.0)

# trajectory output
all = hoomd.group.all()
hoomd.dump.gsd(filename=run_gsd, period=gsd_freq, group=all,
               phase=0)

# integrator setup and run
hoomd.md.integrate.mode_standard(dt=dt)
hoomd.md.integrate.langevin(group=all, kT=1.0, seed=integrate_seed)
hoomd.run_upto(run_length)
