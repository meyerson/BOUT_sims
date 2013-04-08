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
matplotlib.use('Agg')
from read_inp import metadata
import sys
import os
#from post_bout import normed as norm

path=sys.argv[1] #where to look for data
key=sys.argv[2] #some movie id
#key = '_XY'
#path = '/tmp/SOLblobXY/data_dirichlet_precon'


#key='movie'
cache=path+'_movie'

if os.path.exists(cache):
    os.rmdir(cache)
    
os.makedirs(cache)

print path
#meta = metadata(path=path)
#path='/scratch/01523/meyerson/ConvectSOL/data_convect_sol_3.18'



from boutdata import collect
from boututils import savemovie

import numpy as np


nz = np.squeeze(collect("MZ",xind=[0,0],path=path,info=False))
nx =  np.squeeze(collect("NXPE",xind=[0,0],path=path,info=False))*np.squeeze(collect("MXSUB",xind=[0,0],path=path,info=False)) #without gaurds

n = np.squeeze(collect("n",tind=[0,499],path=path,info=False))
u = np.squeeze(collect("u",tind=[0,499],path=path,info=False))
phi = np.squeeze(collect("phi",tind=[0,499],path=path,info=False))
alpha =  np.squeeze(collect("alpha",tind=[0,499],path=path,info=False))

#one movie per cpu
# savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",overcontour=Trueo,xO=xO,yO=yO,dx=dx,dy=dy)
print n.shape
nt,nx,ny = n.shape

dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
dy = np.squeeze(collect("dz",path=path,xind=[0,0]))
yO = -.5*(dy*ny)
xO = -.13 *(dx*nx)


savemovie(n[1:,:,:]+.1,data2=phi[1:,:,:],moviename='movie_n_phi'+key+'.avi',cache=cache+"/",
          overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy)

savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",
          overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy)

savemovie(n[1:,:,ny/2],moviename='movie_n_1D'+key+'.avi',cache=cache+"/", 
          overcontour=True,xO=xO, yO=yO,dx=dx,dy=dy)

#savemovie(n,data2=(n*(n>0)+10000000.*n*(n<0)),moviename='movie_n_neg'+key+'.avi',
#          cache=cache+"/",overcontour=True,xO=xO, yO=yO,dx=dx,dy=dy,norm=False,
#          nlevels = 3, removeZero = False)

#savemovie((n*(n>0)+10000000.*n*(n<0))[:,:,ny/2],moviename='movie_n_1D_neg'+key+'.avi',
#          cache=cache+"/",overcontour=True,xO=xO, yO=yO,dx=dx,dy=dy,norm=False)


#os.rmdir
 
