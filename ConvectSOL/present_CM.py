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
from boututils import savemovie
import numpy as np
from plot_CM import CM_mass, present
import pickle
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib 
matplotlib.use('pdf')

path="/tmp/SOLblobXY/data_dirichlet_Ra1e6_fix"
n = np.squeeze(collect("n",path=path,tind=[0,299]))
#u = np.squeeze(collect("u",path=path))

print 'np.mean(n): ' ,np.mean(n) 

meta={'dx':.2343,'dy':.156,'y0':-10,'x0':-5.46,'dt':1e-1}

nt,nx,ny = n.shape

a = np.zeros((nt,nx,ny),np.float64)

#n = n*(n>1.0000000001e-3)

sim_data = []
z = CM_mass(n,meta=meta,label='BOUT++, Ra=1e6_fix_g')

f_db = open('local_XY_Re=1e6_fixg_db','w')
pickle.dump(z,f_db)
f_db.close()
sim_data.append(z)




#load older runs
filename = 'local_XY_Re=1e4_fixg_db'
f_db = open(filename,'r')
older = pickle.load(f_db)
f_db.close()
sim_data.append(older)

filename = 'local_XY_Re=1e4_topo_db'
f_db = open(filename,'r')
older = pickle.load(f_db)
f_db.close()
sim_data.append(older)


filename = 'local_XY_Re=1e4_BCfix_db'
f_db = open(filename,'r')
older = pickle.load(f_db)
f_db.close()
sim_data.append(older)

filename = 'craigs_data.txt_db'
f_db = open(filename,'r')
arcon_db = pickle.load(f_db)
f_db.close()

sim_data.append(arcon_db)

# filename = 'local_XY_Re=1e4_db'
# f_db = open(filename,'r')
# older = pickle.load(f_db)
# f_db.close()
# sim_data.append(older)

# filename = 'local_XY_Re=1e6_db'
# f_db = open(filename,'r')
# older = pickle.load(f_db)
# f_db.close()


sim_data.append(older)


#make a pdf



#present(sim_data,canvas=garcia)
#load Garcia results
import matplotlib.image as mpimg
garciaX = Image.open('GarciaX.png').transpose(Image.FLIP_TOP_BOTTOM)
garciaV = Image.open('GarciaV.png').transpose(Image.FLIP_TOP_BOTTOM)
garciaX = mpimg.imread('GarciaX.png')

# delta = 0.025
# x = np.arange(0, , delta)
# X, Y = np.meshgrid(x, y)
# Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
# Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)

pp = PdfPages('garcia.pdf')
figX = plt.figure()
allX = figX.add_subplot(1,1,1)
#fig.subplots_adjust(bottom=0.14)
#x = y = np.arange(-3.0, 3.0, delta)
#X, Y = np.meshgrid(x, y)

allX.imshow(garciaX,alpha = .1,extent=[0,10,0,10])
allX.axis('off')
# figX.savefig(pp, format='pdf')
# plt.close(figX)

present(sim_data,pp,xcanvas=allX)
figX.savefig(pp, format='pdf')
plt.close(figX)

figV = plt.figure()
allV = figV.add_subplot(1,1,1)
#fig.subplots_adjust(bottom=0.14)
allV.imshow(garciaV,alpha = .1)
allV.axis('off')
figV.savefig(pp, format='pdf')
plt.close(figV)

pp.close()

pp = PdfPages('cm.pdf')
present(sim_data,pp,compare_png_x=garciaX,compare_png_v=garciaV)
present(sim_data,pp)
pp.close()
