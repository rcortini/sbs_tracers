import numpy as np
import sbs_tracers_analysis as sbs
import os, sys, optparse

program_name = "msd"

def print_usage () :
    print "Usage: %s <gsd> [options]"%program_name
    print 
    print "OPTIONS:"
    print "--teq          <1000>                equilibration time"
    print "--tsample      <1>                   sampling time"
    print "--tracer_text  <name t>              string defining tracers"
    print 

# get parameters from the command line
parser = optparse.OptionParser()
parser.add_option('-t', '--tracer_text', 
                  dest="tracer_text", 
                  default="name t",
                  )
parser.add_option('--teq',
                  dest="teq",
                  default=1000,
                  type="int",
                  )
parser.add_option('--tsample',
                  dest="tsample",
                  default=1,
                  type="int",
                  )
options, remainder = parser.parse_args()

try :
    gsd = remainder[0]
except IndexError :
    sbs.error_message(program_name, "Incorrect usage")
    print_usage()
    sys.exit (1)

# init the HOOMD simulation instance from mybiotools
sim = sbs.hoomdsim(gsd)

# go for the calculation of the MSD
sbs.log_message(program_name, "Calculating MSD")
msd = sbs.msd_t(sim, options.tracer_text, options.teq, options.tsample)

# perform the gradient of the MSD
dinst = np.gradient(msd,axis=1) / 6.0

# save files
sbs.log_message(program_name, "Done. Saving file")
np.save('msd.npy', msd)
np.save('dinst.npy', dinst)
