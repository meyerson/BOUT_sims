#!/opt/apps/python/epd/7.2.2/bin/python
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


import matplotlib 
#matplotlib.use('pdf')
from boutdata import collect
from boututils import *
from post_bout import read_inp,parse_inp,read_cxx
import numpy as np

inp = read_inp(path='./',boutinp='BOUT_turb.inp')
inp = parse_inp(inp)

nx = np.int(inp['[mesh]']['nx'])
nz = np.int(inp['[main]']['MZ'])

path='/tmp/SOLblob/data_Ra1e4_turb'
#dukatpath = '/media/dukat'+path

#path = dukatpath

n = np.squeeze(collect('n',path=path,tind=[0,299]))
phi =  np.squeeze(collect('phi',path=path,tind=[0,299]))
u =  np.squeeze(collect('u',path=path,tind=[0,299]))

#nt,nx,ny = n.shape
#print n.shape
showdata(np.average(n[:,:,:],axis=2))

