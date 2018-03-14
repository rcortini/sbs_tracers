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

def random_quaternion () :
    phi = 2*np.pi*np.random.random()
    z = -1.0 + 2.0*np.random.random()
    s = np.sqrt (1.-z*z)
    a = np.array([s*np.cos(phi),s*np.sin(phi),z])
    angle = np.pi*np.random.rand()
    l = np.sin(angle)
    return np.array([np.cos(angle),a[0]*l,a[1]*l,a[2]*l])

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
sim = int(sys.argv[3])
# random number generator seeds
init_seed = int(sys.argv[4])
integrate_seed = int(sys.argv[5])

###############################
# physical parameters
###############################
kT = 1.0                 # temperature in units of energy
N = 1024                 # number of beads in the polymer
nbs = int(N*phi)         # number of binding sites
sep = 0.35               # minimum separation of beads at init
b_polymer = 1.2          # more than twice the separation
polymer_r0 = 1.2         # polymer bond length, local units
ntracers = 200           # number of tracer molecules
k = 330.0                # polymer stiffness constant
sigma = 1.0              # particle diameter, local units
e_repulsion = 1.0        # repulsion LJ epsilon
e_attraction = 10.0      # binder-to-binding sites affinity (kT)
r_cut = 3.0              # LJ interaction cutoff distance
mt = 1.0                 # mass of t particles (needed for moment of inertia)
dt1 = 0.4                # distance between center of t and t1 particles
sigmat = 1.0             # diameter of t particles
sigmat1 = 0.5            # diameter of t1 particles

###############################
# general parameters
###############################
name = "sbs_tracers-phi-%.2f-e-%.1f-%d"%(phi,e,sim)
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
           particle.orientation = random_quaternion()
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

# set interaction radius of t1-all-the-rest pair
sigmat1_interaction = 0.5 * (sigmat1+sigma)

# neighbor list init
nlist = hoomd.md.nlist.cell ()

# LJ interactions init
lj = hoomd.md.pair.lj(r_cut,nlist)

# interactions of the p particles
lj.pair_coeff.set('p', 'p', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 'a', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 'b', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 't1', epsilon=e, sigma=sigmat1_interaction, alpha=1.0)

# interactions of the a particles
lj.pair_coeff.set('a', 'a', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('a', 'b', epsilon=e_attraction, sigma=sigma, alpha=1.0)
lj.pair_coeff.set('a', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('a', 't1', epsilon=e, sigma=sigmat1_interaction, alpha=1.0)

# interactions of the b particles
lj.pair_coeff.set('b', 'b', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('b', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('b', 't1', epsilon=e_repulsion, sigma=sigmat1_interaction, alpha=0.0)

# interactions of the t particles
lj.pair_coeff.set('t', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('t', 't1', epsilon=e_repulsion, sigma=sigmat1_interaction, alpha=0.0)

# interactions of the t1 particles
lj.pair_coeff.set('t1', 't1', epsilon=e_repulsion, sigma=sigmat1_interaction, alpha=0.0)

# trajectory output
hoomd.dump.gsd(filename=run_gsd, period=gsd_freq, 
               group=hoomd.group.all(),phase=0)

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
