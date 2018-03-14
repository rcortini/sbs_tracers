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

    User should supply "polymer_text" and "tracer_text", which define the
    polymer particles and tracer particles. Typically polymer_text is "type p or
    type a", and tracer_text is "type t".
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
        # calculate ChIP-seq at this time frame
        c = distance_array(polymer.positions,tracers.positions,box=ts.dimensions)
        C += (c<t_threshold)
        Ct = C.sum(axis=1)
        DKL_t[i] = KL_divergence(Ct,Rt)
    # coverage analysis
    C[C>1] = 1
    coverage = C.sum(axis=0).astype('float')/N
    return DKL_t,H,Ct.astype(np.int64),coverage
