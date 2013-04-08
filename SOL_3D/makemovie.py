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
from boututils import savemovie,showdata

import numpy as np
#gamma = collect("Gammax",yi5Cnd=[5,5],path=path)
#n = np.squeeze(collect("n",yind=[2,2],path=path))
#u = np.squeeze(collect("u",yind=[2,2],path=path))
#phi = np.squeeze(collect("phi",yind=[2,2],path=path))
# n = np.squeeze(collect("n",path=path))
# u = np.squeeze(collect("u",path=path))
# phi = np.squeeze(collect("phi",path=path))
          

          #one movie per cpu
          

nz = np.squeeze(collect("MZ",xind=[0,0],path=path,info=False))
nx =  np.squeeze(collect("NXPE",xind=[0,0],path=path,info=False))*np.squeeze(collect("MXSUB",xind=[0,0],path=path,info=False)) #without gaurds

#nz = 65
#nx=512
#nz = 129
print nx,nz

n = np.squeeze(collect("n",xind=[2,nx],yind=[16,16],tind=[0,36],path=path,info=False))
u = np.squeeze(collect("u",xind=[2,nx],yind=[16,16],tind=[0,36],path=path,info=False))
phi = np.squeeze(collect("phi",xind=[2,nx],yind=[16,16],tind=[0,36],path=path,info=False))

n_xy = np.squeeze(collect("n",zind=[64,64],tind=[1,400],path=path,info=False))
u_xy = np.squeeze(collect("u",zind=[64,64],tind=[1,400],path=path,info=False))
phi_xy = np.squeeze(collect("phi",zind=[64,64],tind=[1,400],path=path,info=False))


#one movie per cpu
savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",overcontour=True)
print n.shape
nt,nx,ny = n.shape

#try:
##    dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
#    dy = np.squeeze(collect("dz",path=path,xind=[0,0]))
#except:B
dx = 1.
dy = 1.

print dx,dy
yO = -.5*(dy*ny)
xO = -.16 *(dx*nx)


savemovie(n[1:,:,:]+.1,data2=phi[1:,:,:],moviename='movie_n_phi'+key+'.avi',cache=cache+"/",
          overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy)

#savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",
#          overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy)

#savemovie(n[1:,:,ny/2],moviename='movie_n_1D'+key+'.avi',cache=cache+"/",  
#          overcontour=True,xO=xO, yO=yO,dx=dx,dy=dy)

# savemovie(n_xy,data2=phi_xy,moviename='movie_nxy_phi'+key+'.avi',
#           cache=cache+"/",overcontour=True)

#savemovie(n,data2=(n*(n>0)+10000000.*n*(n<0)),moviename='movie_n_neg'+key+'.avi',
#          cache=cache+"/",overcontour=True,xO=xO, yO=yO,dx=dx,dy=dy,norm=False,
#          nlevels = 3, removeZero = False)

#savemovie((n*(n>0)+10000000.*n*(n<0))[:,:,ny/2],moviename='movie_n_1D_neg'+key+'.avi',
#          cache=cache+"/",overcontour=True,xO=xO, yO=yO,dx=dx,dy=dy,norm=False)


#os.rmdir
 
