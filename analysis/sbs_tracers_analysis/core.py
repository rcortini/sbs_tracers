import os
import MDAnalysis as mda

class hoomdsim :
    def __init__ (self,topology_file) :
        self.u = mda.Universe (topology_file)
