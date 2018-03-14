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
    print "--msd_out      <msd>                 Hi-C output file name"
    print "--tracer_text  <name t>              string defining tracers"
    print 

# get parameters from the command line
parser = optparse.OptionParser()
parser.add_option('-m', '--msd_out', 
                  dest="msd_out", 
                  default="msd",
                  )
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
    mbt.error_message(program_name, "Incorrect usage")
    print_usage()
    sys.exit (1)

# init the HOOMD simulation instance from mybiotools
sim = mbt.hoomdsim(gsd)

# go for the calculation of the MSD
mbt.log_message(program_name, "Calculating MSD")
sbs.msd_t(sim, options.tracer_text, options.teq, options.tsample)

# save files
mbt.log_message(program_name, "Done. Saving file")
np.save(options.msd_out,sim.msd_t)
