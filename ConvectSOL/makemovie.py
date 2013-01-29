#! /share/home/01523/meyerson/local/bin/python
import sys
# sys.path.append('/home/cryosphere/BOUT/tools/pylib')
# sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
# sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
# sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')
# sys.path.append('/usr/local/pylib')
sys.path.append('/share/home/01523/meyerson/lib/python')
sys.path.append('/share/home/01523/meyerson/pylib')
sys.path.append('/share/home/01523/meyerson/BOUT/tools/pylib')
sys.path.append('/share/home/01523/meyerson/BOUT/tools/pylib/boutdata')
sys.path.append('/share/home/01523/meyerson/BOUT/tools/pylib/boututils')
sys.path.append('/share/home/01523/meyerson/BOUT/tools/pylib/post_bout')
#import matplotlib
#matplotlib.use('Agg')
from read_inp import metadata
import sys
import os
#from post_bout import normed as norm

path=sys.argv[1]
key=sys.argv[2]

#key='movie'
cache='/scratch/01523/meyerson/'+key

if os.path.exists(cache):
    os.rmdir(cache)
    
os.makedirs(cache)

print path
#meta = metadata(path=path)
#path='/scratch/01523/meyerson/ConvectSOL/data_convect_sol_3.18'
dx = 35.0/512
dy = 20.0/1025
yO = -10.0
xO = -5.0

from boutdata import collect
from boututils import savemovie
from post_bout import normed as norm
import numpy as np
#gamma = collect("Gammax",yi5Cnd=[5,5],path=path)
#n = np.squeeze(collect("n",yind=[2,2],path=path))
#u = np.squeeze(collect("u",yind=[2,2],path=path))
#phi = np.squeeze(collect("phi",yind=[2,2],path=path))
n = np.squeeze(collect("n",yind=[2,2],path=path))
u = np.squeeze(collect("u",yind=[2,2],path=path))
phi = np.squeeze(collect("phi",yind=[2,2],path=path))

#one movie per cpu
# savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy)

savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy,norm=False)


# savemovie(n[1:,:,:],data2=phi[1:,:,:],moviename='movie_n_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy)

savemovie(u[1:,:,:],data2=phi[1:,:,:],moviename='movie_u_phi'+key+'.avi',cache=cache+"/",overcontour=True,xO=xO,yO=yO,dx=dx,dy=dy,norm=False)

#os.rmdir
