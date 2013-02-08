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


#import matplotlib
#matplotlib.use('Agg')
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
dx = 35.0/128
dy = 20.0/128
yO = -10.0
xO = -5.0

from boutdata import collect
from boututils import savemovie

import numpy as np
#gamma = collect("Gammax",yi5Cnd=[5,5],path=path)
#n = np.squeeze(collect("n",yind=[2,2],path=path))
#u = np.squeeze(collect("u",yind=[2,2],path=path))
#phi = np.squeeze(collect("phi",yind=[2,2],path=path))
n = np.squeeze(collect("n",path=path))
u = np.squeeze(collect("u",path=path))
phi = np.squeeze(collect("phi",path=path))

#one movie per cpu
# savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy)
print u.shape

#computer the center of mass

nt,nx,ny = n.shape

#savemovie(n[1:,:,:],data2=phi[1:,:,:],moviename='movie_n_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy,norm=False)


# savemovie(n[1:,:,:],data2=phi[1:,:,:],moviename='movie_n_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy)

#savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy,norm=False)

savemovie(u[1:,:,ny/2],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy,norm=False)

#os.rmdir
