print "Loading BOUT++ blob post processing routines"

# Load routines from separate files
import sys
import os

HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')


try:
    from blob_info import blob_info,turb_info
except:
    print "Sorry, no read_grid"
