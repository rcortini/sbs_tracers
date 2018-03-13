from __future__ import print_function
import os
import numpy as np
import hoomd, hoomd.deprecated, hoomd.md
import mybiotools as mbt
import sys

# check for proper invocation
if len(sys.argv) < 6 :
    print("python sbs_tracers.py <phi> <e> <n> <init_seed> <integrate_seed>")
    sys.exit(1)

# initialize hoomd-blue
hoomd.context.initialize()

def particle_images(sim,frame_id) :
    """
    Get the image index of all particles in simulation, at the frame 'frame_id'
    """
    # get positions of all particles: define first the atom selection, then jump to
    # the user-requested trajectory frame, get the box dimensions (currently works
    # only for orthorhombic boxes, then calculate the image indices
    atoms = sim.u.select_atoms ('all')
    ts = sim.u.trajectory[frame_id]
    L = ts.dimensions[:3]
    pos = atoms.positions + L/2.
    return pos//L

def restore_images (images,system) :
    """
    Restore the images into the system defined by "system"
    """
    # get information on the image indices of each particle from the images file
    snapshot = system.take_snapshot()
    for i,im in enumerate (images) :
        snapshot.particles.image[i] = list (im)
    system.restore_snapshot (snapshot)

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
ntracers = 10            # number of tracer molecules
k = 330.0                # polymer stiffness constant
sigma = 1.0              # particle diameter, local units
e_repulsion = 1.0        # repulsion LJ epsilon
e_attraction = 10.0      # binder-to-binding sites affinity (kT)
r_cut = 3.0              # LJ interaction cutoff distance

###############################
# general parameters
###############################
name = "sbs_tracers-phi-%.2f-e-%.1f-%d"%(phi,e,sim)
xml_freq = 10000
dcd_freq = 10000
run_xml = "%s.xml"%name
run_dcd = "%s.dcd"%name
run_images = "%s-images.dat"%name
init_xml = "%s-init.xml"%name

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
# polymer is made of 'p' particles and 'a' particles:
# respectively: polymer and anchors
polymer_type = ['p']*N
np.random.seed(seed=init_seed)
bs = np.random.choice(N,size=nbs)
for i in bs :
    polymer_type[int(i)] = 'a'
polymer1 = dict(bond_len=b_polymer, type=polymer_type,
                bond="linear", count=1)

###############################
# init diffusing beads
###############################
# binders: 'b' particles
factors = dict(bond_len=b_polymer, type=['b'],
                bond="linear", count=nfactors)
# tracers: 't' particles
tracers = dict(bond_len=b_polymer, type=['t'],
                bond="linear", count=ntracers)

###############################
# system init
###############################
# check if trajectory file exists
if os.path.exists(run_dcd) :
    # if so, restore the system from last step of the simulation
    system = hoomd.deprecated.init.read_xml(run_xml)
    sim = mbt.hoomdsim(run_xml,run_dcd)
    run_images = particle_images(sim,-1)
    restore_images(run_images,system)
else :
    hoomd.deprecated.init.create_random_polymers(box=box, polymers=[polymer1,factors,tracers],
                                separation=dict(p=sep,
                                                a=sep,
                                                b=sep,
                                                t=sep), seed=init_seed)

# polymer bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set('polymer', k=k, r0=polymer_r0)

# polymer self interactions
nlist = hoomd.md.nlist.cell ()
lj = hoomd.md.pair.lj(r_cut,nlist)

# interactions of the p particles
lj.pair_coeff.set('p', 'p', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 'a', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 'b', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('p', 't', epsilon=e, sigma=sigma, alpha=1.0)

# interactions of the a particles
lj.pair_coeff.set('a', 'a', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('a', 'b', epsilon=e_attraction, sigma=sigma, alpha=1.0)
lj.pair_coeff.set('a', 't', epsilon=e, sigma=sigma, alpha=1.0)

# interactions of the b particles
lj.pair_coeff.set('b', 'b', epsilon=e_repulsion, sigma=sigma, alpha=0.0)
lj.pair_coeff.set('b', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)

# interactions of the t particles
lj.pair_coeff.set('t', 't', epsilon=e_repulsion, sigma=sigma, alpha=0.0)

# trajectory output
all = hoomd.group.all()
xml = hoomd.deprecated.dump.xml(all,filename=run_xml, vis=True, period=xml_freq, restart=True,
               phase=0)
dcd = hoomd.dump.dcd(filename=run_dcd, period=dcd_freq, unwrap_full=True, phase=0)

###############################
# RUN!
###############################

# integrate NVT for a bunch of time steps
hoomd.md.integrate.mode_standard(dt=dt)
hoomd.md.integrate.langevin(group=all, kT=kT, seed=integrate_seed)
hoomd.run_upto(run_length)
