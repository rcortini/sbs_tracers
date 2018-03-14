import numpy as np
import sbs_tracers_analysis as sbs
import os, sys, optparse

program_name = "tracers_analysis"

def print_usage () :
    print "Usage: %s <gsd> [options]"%program_name
    print 
    print "OPTIONS:"
    print "--threshold       <2.0>                 threshold for considering contacts"
    print "--t_threshold     <2.0>                 threshold for tracer contacts"
    print "--p_threshold     <2.0>                 threshold for polymer contacts"
    print "--teq             <5000>                equilibration time"
    print "--tsample         <50>                  sampling time"
    print "--DKL_t_out       <DKL_t>               DKL_t output file name"
    print "--hic_out         <hic>                 Hic output file name"
    print "--chipseq_out     <chipseq>             Chipseq output file name"
    print "--coverage_out    <coverage>            Coverage output file name"
    print "--polymer_text    <name A or name B>    string defining polymer"
    print "--tracer_text     <name D>              string defining tracers"
    print 

# get parameters from the command line
parser = optparse.OptionParser()
parser.add_option('-D', '--DKL_t_out', 
                  dest="DKL_t_out", 
                  default="DKL_t",
                  )
parser.add_option('-H', '--hic_out', 
                  dest="hic_out", 
                  default="hic",
                  )
parser.add_option('-C', '--chipseq_out', 
                  dest="chipseq_out", 
                  default="chipseq",
                  )
parser.add_option('-c', '--coverage_out', 
                  dest="coverage_out", 
                  default="coverage",
                  )
parser.add_option('-p', '--polymer_text', 
                  dest="polymer_text", 
                  default="type p or type a",
                  )
parser.add_option('-t', '--tracer_text', 
                  dest="tracer_text", 
                  default="type t",
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
parser.add_option('--threshold',
                  dest="threshold",
                  type="float",
                  )
parser.add_option('--t_threshold',
                  dest="t_threshold",
                  type="float",
                  )
parser.add_option('--p_threshold',
                  dest="p_threshold",
                  type="float",
                  )
options, remainder = parser.parse_args()

try :
    gsd = remainder[0]
except IndexError :
    sbs.error_message(program_name, "Incorrect usage")
    print_usage()
    sys.exit(1)

# now check that EITHER threshold is set OR (t_threshold AND p_threshold) are
# set
if options.threshold is None :
    if options.p_threshold is not None and options.t_threshold is not None :
        p_threshold = options.p_threshold
        t_threshold = options.t_threshold
    else :
        # error
        sbs.error_message(program_name,"Set either threshold or"
                          " p_threshold and t_threshold")
        sys.exit(1)
else :
    if options.p_threshold is not None or options.t_threshold is not None :
        # error
        sbs.error_message(program_name, "Set either threshold or"
                          " p_threshold and t_threshold")
        sys.exit(1)
    else :
        p_threshold = options.threshold
        t_threshold = options.threshold

# init the HOOMD simulation instance from mybiotools
sim = sbs.hoomdsim(gsd)

# go for the calculation of the DKL_t
sbs.log_message(program_name, "Starting tracer analysis")
DKL_t,H,C,coverage = sbs.tracers_analysis (sim,
                    options.polymer_text,
                    options.tracer_text,
                    options.teq,
                    options.tsample,
                    t_threshold,
                    p_threshold)

# save files
sbs.log_message(program_name, "Done. Saving files")
np.save(options.DKL_t_out,DKL_t)
np.save(options.hic_out,H)
np.save(options.chipseq_out,C)
np.save(options.coverage_out,coverage)
