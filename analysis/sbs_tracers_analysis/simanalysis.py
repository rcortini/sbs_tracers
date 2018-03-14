import numpy as np
from MDAnalysis.analysis.distances import distance_array
from .math import KL_divergence

def traj_nslice (u,teq,tsample) :
    """
    Returns the number of frames in the trajectory in universe u, using teq as
    equilibration time and tsample as sampling time
    """
    # get the number of frames in the slice (http://stackoverflow.com/a/7223557)
    traj_slice = u.trajectory[teq::tsample]
    return sum(1 for _ in traj_slice)

def tracers_analysis (sim,polymer_text,tracer_text,teq,tsample,t_threshold,p_threshold) :
    """
    This function does the complete analysis of the tracers in the simulation.
    It calculates the polymer contact matrix, traffic of tracers, Kullback-Leibler
    divergence between the two profiles as a function of time, and coverage of
    the tracers.
    
    "sim" is the "hoomdsim" class defined in core.

    User should supply "polymer_text" and "tracer_text", which define the
    polymer particles and tracer particles. Typically polymer_text is "type p or
    type a", and tracer_text is "type t".

    Sampling is done starting from time `teq`, and every `tsample` time steps.

    The threshold for the definition of contacts between the tracers and the
    polymer is `t_threshold`. The one between the polymer monomers is
    `p_threshold`.
    """
    # define DKL(t) vector
    nframes = traj_nslice(sim.u,teq,tsample)
    DKL_t = np.zeros(nframes)
    # define polymer and tracers
    polymer = sim.u.select_atoms(polymer_text)
    tracers = sim.u.select_atoms(tracer_text)
    N = polymer.n_atoms
    ntracers = tracers.n_atoms
    # init H and C vectors
    H = np.zeros((N,N),dtype=np.int32)
    C = np.zeros((N,ntracers),dtype=np.int32)
    # analyze all simulation frames as decided
    for i,ts in enumerate(sim.u.trajectory[teq::tsample]) :
        # calculate Hi-C at this time frame
        d = distance_array(polymer.positions,polymer.positions,box=ts.dimensions)
        H += (d<p_threshold)
        Rt = H.sum(axis=1)
        # calculate traffic at this time frame
        c = distance_array(polymer.positions,tracers.positions,box=ts.dimensions)
        C += (c<t_threshold)
        Ct = C.sum(axis=1)
        DKL_t[i] = KL_divergence(Ct,Rt)
    # coverage analysis
    C[C>1] = 1
    coverage = C.sum(axis=0).astype('float')/N
    return DKL_t,H,Ct.astype(np.int64),coverage

def msd_t(sim,particles_text,teq,tsample) :
    """
    Calculate the mean square displacement of the particles defined by
    'particles_text' in simulation `sim`, using sampling `tsample` and equilibration
    time `teq`. Returns the matrix corresponding to the mean square displacement
    of each particle, along with a matrix corresponding to the variance in the
    estimate of this quantity.

    This calculation is detailed in Qian, H., M. P. Sheetz, and E. L. Elson. 1991.
    ``Single Particle Tracking. Analysis of Diffusion and Flow in
    Two-Dimensional Systems.” Biophysical Journal 60 (4):910–21.
    """
    u = sim.u
    particles = u.select_atoms(particles_text)
    nparticles = particles.n_atoms
    nslice = traj_nslice (u,teq,tsample)
    # initialize the matrix containing all the positions
    # of the particles at all the sampling frames
    particles_pos = np.zeros ((nslice,nparticles,3))
    for i,ts in enumerate(u.trajectory[teq::tsample]) :
        particles_pos[i,:,:] = particles.positions
    # now initialize the Delta matrix, which contains the
    # squared differences between the particles' positions
    # at different time delays
    Nt = int(nslice/2)
    Delta = np.zeros((nparticles,Nt,Nt))
    for delay in xrange(1,Nt+1) :
        for t0 in xrange (Nt) :
            t1 = t0 + delay
            pos1 = particles_pos[t1,:,:]
            pos0 = particles_pos[t0,:,:]
            Delta[:,delay-1,t0] = np.sum((pos1-pos0)**2,axis=1)
    # return the matrices of MSD and its variance
    return np.mean(Delta,axis=2),np.var(Delta,axis=2)
