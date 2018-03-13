import os
import MDAnalysis as mda

class hoomdsim :
    def __init__ (self,topology_file,dcd=None) :
        if dcd is not None :
            if not os.path.exists(dcd) and topology_file.endswith('.gsd'):
                u = mda.Universe (topology_file)
            else :
                u = mda.Universe (topology_file,dcd,format='DCD')
        else :
            u = mda.Universe (topology_file)
        self.u = u
