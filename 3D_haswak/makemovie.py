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
import matplotlib
matplotlib.use('Agg')
from read_inp import metadata
import sys
import os


path=sys.argv[1]
key=sys.argv[2]


cache='/scratch/01523/meyerson/'+key

if os.path.exists(cache):
    os.rmdir(cache)
    
os.makedirs(cache)

print path
meta = metadata(path=path)

from boutdata import collect
from boututils import savemovie
gamma = collect("Gammax",yind=[5,5],path=path)
ni = collect("Ni",yind=[5,5],path=path)
rho = collect("rho",yind=[5,5],path=path)

#one movie per cpu
print ni.shape, gamma.shape
savemovie(gamma[:,:,0,:],ni[:,:,0,:],moviename='movie_'+key+'.avi',
          meta=meta,cache=cache+"/",overcontour=True)

savemovie(rho[:,:,0,:],ni[:,:,0,:],moviename='movie_phi_ni'+key+'.avi',
          meta=meta,cache=cache+"/",overcontour=True)

#os.rmdir(cache)

