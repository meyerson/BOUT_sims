#!/usr/bin/python
import sys, os
# sys.path.append('/home/cryosphere/BOUT/tools/pylib')
# sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
# sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
# sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')
HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')

sys.path.append('/usr/local/pylib')
sys.path.append(HOME+'/lib/python')
sys.path.append(HOME+'/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib/boutdata')
sys.path.append(BOUT_TOP+'/tools/pylib/boututils')
sys.path.append(BOUT_TOP+'/tools/pylib/post_bout')

from read_inp import metadata
import sys
import os
from boutdata import collect
from boututils import *
from post_bout import read_grid
import numpy as np
#from plot_CM import CM_mass, present
import pickle
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib 
matplotlib.use('pdf')

#output_grid.nc
def get_IC(file='output_grid.nc'):
    IC = read_grid(gridfile=file)
    

    meta={}
    
    for elem in IC.variables:
        print elem
        meta[elem] = np.array(IC.variables[elem][:])
  
    return meta


IC_rmp = get_IC()
print IC_rmp['nx']
print IC_rmp['ny']
#BG = get_IC(file='')

             
