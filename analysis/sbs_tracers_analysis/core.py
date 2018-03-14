import os
import numpy as np
import MDAnalysis as mda

class hoomdsim :
    def __init__ (self,topology_file) :
        self.u = mda.Universe (topology_file)

def random_quaternion () :
    phi = 2*np.pi*np.random.random()
    z = -1.0 + 2.0*np.random.random()
    s = np.sqrt (1.-z*z)
    a = np.array([s*np.cos(phi),s*np.sin(phi),z])
    angle = np.pi*np.random.rand()
    l = np.sin(angle)
    return np.array([np.cos(angle),a[0]*l,a[1]*l,a[2]*l])
