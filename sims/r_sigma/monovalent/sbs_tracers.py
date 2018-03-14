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
r_cut = 2.5*sigma        # LJ interaction cutoff distance
mt = 1.0                 # mass of t particles (needed for moment of inertia)

# geometry of the composite particles
delta = 0.15             # how much the small D1 particles emerge from D
sigmat1 = 0.5            # diameter of D1 particles
dt1 = delta - sigmat1/2. + sigmat/2.

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
    # moment of inertia of the t particles
    Ix = Iy = Iz = 2.0/5.0 * mt*(sigmat/2.0)**2
    # init random number generator seed
    np.random.seed(init_seed)
    # fix t particle moment of inertia and orientations
    for particle in system.particles :
       if particle.type == 't' :
           particle.mass = mt
           particle.diameter = sigmat
           particle.moment_inertia = [Ix,Iy,Iz]
           particle.orientation = sbs.random_quaternion()
    # Add consituent particles of type t1 and create the particles
    system.particles.types.add('t1')

# create the rigid bodies
rigid = hoomd.md.constrain.rigid()
rigid.set_param('t',
               types=['t1'],
               positions=[(dt1,0,0)],
               diameters=[sigmat1])
rigid.create_bodies()

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
lj.pair_coeff.set('p', 't1', epsilon=e, sigma=LJsigma(sigma,sigmat1), alpha=1.0)

# interactions of the a particles
lj.pair_coeff.set('a', 'a', epsilon=e_repulsion, sigma=LJsigma(sigma,sigma), alpha=0.0)
lj.pair_coeff.set('a', 'b', epsilon=e_attraction, sigma=LJsigma(sigma,sigma), alpha=1.0)
lj.pair_coeff.set('a', 't', epsilon=e_repulsion, sigma=LJsigma(sigma,sigmat), alpha=0.0)
lj.pair_coeff.set('a', 't1', epsilon=e, sigma=LJsigma(sigma,sigmat1), alpha=1.0)

# interactions of the b particles
lj.pair_coeff.set('b', 'b', epsilon=e_repulsion, sigma=LJsigma(sigma,sigma), alpha=0.0)
lj.pair_coeff.set('b', 't', epsilon=e_repulsion, sigma=LJsigma(sigma,sigmat), alpha=0.0)
lj.pair_coeff.set('b', 't1', epsilon=e_repulsion, sigma=LJsigma(sigma,sigmat1), alpha=0.0)

# interactions of the t particles
lj.pair_coeff.set('t', 't', epsilon=e_repulsion, sigma=LJsigma(sigmat,sigmat), alpha=0.0)
lj.pair_coeff.set('t', 't1', epsilon=e_repulsion, sigma=LJsigma(sigmat,sigmat1), alpha=0.0)

# interactions of the t1 particles
lj.pair_coeff.set('t1', 't1', epsilon=e_repulsion, sigma=LJsigma(sigmat1,sigmat1), alpha=0.0)

# trajectory output
hoomd.dump.gsd(filename=run_gsd, period=gsd_freq, group=hoomd.group.all(),
               phase=0)

# define which particles to integrate: all of them except the P particles
p_group = hoomd.group.type('p')
a_group = hoomd.group.type('a')
b_group = hoomd.group.type('b')
t_group = hoomd.group.type('t')
polymer_group = hoomd.group.union(name='polymer',a=p_group,b=a_group)
factors_group = hoomd.group.union(name='factors',a=b_group,b=t_group)
integrate_group = hoomd.group.union(name='integrate_group',a=polymer_group,
                                    b=factors_group)

# integrator setup and run
hoomd.md.integrate.mode_standard(dt=dt)
hoomd.md.integrate.langevin(group=integrate_group, kT=1.0, seed=integrate_seed)
hoomd.run_upto(run_length)
