import os
import sys
import glob

#####################################
# Path to the cooling table
#####################################
filebase = sys.argv[1]

########################################################################################################
# Get a list of all hdf5 tables in the cooling table path
########################################################################################################
runname_loc = [os.path.basename(x) for x in sorted(glob.glob('%s/*.hdf5'%(filebase)))]
runname = [os.path.splitext(val)[0] for val in runname_loc]
