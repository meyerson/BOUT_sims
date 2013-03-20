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

#data_dir=`echo $PWD|sed "s/home1/scratch/g"`

 #where to look for data
sim_data =[]


def update():
    mu_list = ['1e-3','1e-2','1e-1']
    NP = 264
    

    for mu in mu_list:
        data_dir = '/scratch/01523/meyerson/BOUT_sims/ConvectSOL'
        label=str(NP)+'_mu='+mu+'_HD'
        sim_key = 'convect_sol_XZ_'+mu+'_bigger'

        path=data_dir+'/data_'+sim_key

        print path
        try:
            n = np.squeeze(collect("n",path=path,tind=[0,299]))
        except:
            n = np.squeeze(collect("n",path=path,tind=[0,200]))

        nt,nx,ny = n.shape

        dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
        dy = np.squeeze(collect("dz",path=path,xind=[0,0]))
        yO = -.5*(dy*ny)
        xO = -.17949 *(dx*nx)


        meta={'dx':dx,'dy':dy,'y0':yO,'x0':xO,'dt':1e-1}
        nt,nx,ny = n.shape
        
   
        z = CM_mass(n,meta=meta,label=label)
    
        f_db = open('TACC_'+label+'_db','w')
        pickle.dump(z,f_db)
        f_db.close()
        sim_data.append(z)
        
def render():
    older_runs = ['TACC_264_mu=1e-1_HD_db',
                  'TACC_264_mu=1e-2_HD_db','TACC_264_mu=1e-3_HD_db']

    for run in older_runs:
        f_db = open(run,'r')
        older = pickle.load(f_db)
        f_db.close()
        sim_data.append(older)

    import matplotlib.image as mpimg
    garciaX = Image.open('GarciaX.png').transpose(Image.FLIP_TOP_BOTTOM)
    garciaV = Image.open('GarciaV.png').transpose(Image.FLIP_TOP_BOTTOM)
    garciaX = mpimg.imread('GarciaX.png')
    garciaNamp = mpimg.imread('GarciaNamp.png')

    pp = PdfPages('garcia.pdf')
    figX = plt.figure()
    allX = figX.add_subplot(1,1,1)

    figV = plt.figure()
    allV = figV.add_subplot(1,1,1)
#fig.subplots_adjust(bottom=0.14)
    allV.imshow(garciaV,alpha = .1)
    allV.axis('off')
    figV.savefig(pp, format='pdf')
    plt.close(figV)

    pp.close()
    
    pp = PdfPages('cm.pdf')
    present(sim_data,pp,compare_png_x=garciaX,compare_png_v=garciaV,
            compare_png_max=garciaNamp)
    present(sim_data,pp)
    pp.close()

try:
    update_db = bool(sys.argv[1]) 
except:
    update_db = False

if update_db:
    update()

print update_db
render()
