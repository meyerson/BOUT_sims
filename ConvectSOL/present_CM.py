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


sim_data =[]


def update():
    mu_list = ['1e-3','1e-2','1e-1']
    NP = 264
    

    for mu in mu_list:
        data_dir = '/scratch/01523/meyerson/BOUT_sims/ConvectSOL'
        label=str(NP)+'_mu='+mu+'_trash'
        sim_key = 'convect_sol_XZ_'+mu+'_trash'

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
    older_runs = ['TACC_264_mu=1e-2_HD_db','TACC_264_mu=1e-3_HD_db']

    for run in older_runs:
        f_db = open(run,'r')
        older = pickle.load(f_db)
        f_db.close()
        sim_data.append(older)

    import matplotlib.image as mpimg
    garciaX = Image.open('GarciaX.png').transpose(Image.FLIP_TOP_BOTTOM)
    garciaV = Image.open('GarciaV.png')
    garciaX = mpimg.imread('GarciaX.png')
    garciaNamp = mpimg.imread('GarciaNamp.png')

    pp = PdfPages('garcia.pdf')
    figX = plt.figure()
    allX = figX.add_subplot(1,1,1)

    figV = plt.figure()
    allV = figV.add_subplot(1,1,1)

# sim_key='Ra1e4_d'
# # sim_key='neumann_Ra1e4_nue'
# # sim_key='dirichlet_Ra1e4_core'
# path="/tmp/SOLblobXZ/data_"+sim_key

# print path
# n = np.squeeze(collect("n",path=path,tind=[0,299]))
# print 'np.mean(n): ' ,np.mean(n) 

# #meta={'dx':1,'dy':1,'y0':-20,'x0':-23,'dt':1e-1}
# meta={'dx':.2343,'dy':.156,'y0':-10,'x0':-5.46,'dt':1e-1}
# nt,nx,ny = n.shape

# sim_data = []
# z = CM_mass(n,meta=meta,label='XZ_'+sim_key+'_dirich')

# f_db = open('local_XZ_'+sim_key+'_db','w')
# pickle.dump(z,f_db)
# f_db.close()
# sim_data.append(z)



sim_data = []
#load older runs
# older_runs = ['local_XY_Re=1e4_iso_db','local_XY_Re=1e4_x2_db',
#               'local_XY_Re=1e8_fixg_db','local_XY_Re=1e6_fixg_db',
#               'local_XY_Re=1e4_fixg_db','local_XY_Re=1e4_topo_db',
#               'local_XY_Re=1e4_BCfix_db', 'craigs_data.txt_db']
#older_runs = ['craigs_data.txt_db','local_XZ_Ra1e4_db',
#              'local_XY_Ra1e4_db']

#update()

older_runs = ['TACC_264_mu=1e-3_HD_db','TACC_264_mu=1e-2_HD_db']

# if refresh_db:
#     for elem in older_runs:
#         z = CM_mass(n,meta=meta,label=sim_key)
#         f_db = open('local_XY_'+sim_key+'_db','w')
#         pickle.dump(z,f_db)
#         f_db.close()

for run in older_runs:
    f_db = open(run,'r')
    older = pickle.load(f_db)
    f_db.close()
    sim_data.append(older)

#make a pdf



#present(sim_data,canvas=garcia)
#load Garcia results 
import matplotlib.image as mpimg
import csv

garciaX = Image.open('GarciaX.png').transpose(Image.FLIP_TOP_BOTTOM)
#garciaV = Image.open('GarciaV.png').transpose(Image.FLIP_TOP_BOTTOM)
garciaV = mpimg.imread('GarciaV.png')
garciaX = mpimg.imread('GarciaX.png')
garciaNamp = mpimg.imread('GarciaNamp.png')



xfiles = ['GarciaX_Ra1e2.csv','GarciaX_Ra1e4.csv','GarciaX_Ra1e6.csv']
vfiles = ['GarciaV_Ra1e2.csv','GarciaV_Ra1e4.csv','GarciaV_Ra1e6.csv']
amp_files = ['GarciaNamp_Ra1e2.csv','GarciaNamp_Ra1e4.csv','GarciaNamp_Ra1e6.csv']

with open('GarciaX_Ra1e6.csv', 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    spamreader.next()
    for row in spamreader:
        x,y = row[0].split(',')
        print x,y


#exit()
# pp = PdfPages('garcia.pdf')
# figX = plt.figure()
# allX = figX.add_subplot(1,1,1)
# #fig.subplots_adjust(bottom=0.14)
# #x = y = np.arange(-3.0, 3.0, delta)
# #X, Y = np.meshgrid(x, y)

# # allX.imshow(garciaX,alpha = .1,extent=[0,10,0,10])
# # allX.axis('off')
# # present(sim_data,pp,xcanvas=allX)
# # figX.savefig(pp, format='pdf')
# # plt.close(figX)

# figV = plt.figure()
# allV = figV.add_subplot(1,1,1)

# #fig.subplots_adjust(bottom=0.14)
# allV.imshow(garciaV,alpha = .1)
# allV.axis('off')
# figV.savefig(pp, format='pdf')
# plt.close(figV)

# pp.close()

pp = PdfPages('cm.pdf')
# present(sim_data,pp,compare_png_x=garciaX,compare_png_v=garciaV,
#         compare_png_max=garciaNamp)



present(sim_data,pp)
pp.close()
