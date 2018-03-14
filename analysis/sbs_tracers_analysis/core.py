from __future__ import print_function
import time, sys
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

def time_string () :
    return time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime ())

def error_message (program_name, message) :
    full_message = "%s %s: ERROR: %s"%(time_string (), program_name, message)
    print (full_message, file=sys.stderr)

def log_message (program_name, message) :
    full_message = "%s %s: INFO: %s"%(time_string (), program_name, message)
    print (full_message)

def warn_message (program_name, message) :
    full_message = "%s %s: WARNING: %s"%(time_string (), program_name, message)
    print (full_message)
