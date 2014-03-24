import os, sys, inspect
import sqlite3 as sql
import pickle as pkl
from datetime import datetime

HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')
utc = datetime.utcnow()

import matplotlib
matplotlib.use('Agg')
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

sys.path.append('/usr/local/pylib')
sys.path.append(HOME+'/lib/python')
sys.path.append(HOME+'/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib/boutdata')
sys.path.append(BOUT_TOP+'/tools/pylib/boututils')
sys.path.append(BOUT_TOP+'/tools/pylib/post_bout')

cmd_folder = HOME+'/BOUT_sims/blob_py'

print cmd_folder

if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

from blob_info import blob_info, Blob2D
from blob_draw import BlobDraw
from blob_present import BlobPresent
from blob_movie import BlobMovie
from blob_db import BlobDB

import numpy as np
from boutdata import collect
from read_inp import read_inp, parse_inp
#prepare a list of directories
sim_key='Ra1e6'
path="/tmp/SOLblobXZ/data_"+sim_key

# inp = read_inp(path=path,boutinp='BOUT.inp')
# inp = parse_inp(inp)

# inp2 = {}
# for k, v in inp.items():
#      inp2[k] = v

# print inp2

#read the data and process
t1= 0
t2 = 200
try:
     inp = read_inp(path=path,boutinp='BOUT.inp')
     inp = parse_inp(inp)
     
     dx = np.double(inp['[mesh]']['dx'])
     #print 'dx',dx
     zmax = np.double(inp['[main]']['ZMAX'])

     print dx,zmax
     n = (np.squeeze(collect("n",path=path,tind =[t1,t2])))
     print dx,zmax
     n0 = (np.squeeze(collect("n0",path=path,tind =[t1,t2])))
     u = np.squeeze(collect("u",path=path,tind =[t1,t2]))
     phi = np.squeeze(collect("phi",path=path,tind =[t1,t2]))
     time = np.squeeze(collect("t_array",path=path,xind=[0,0]))
     #dx = np.squeeze(collect("dx",path=path))
     #zmax = np.squeeze(collect("ZMAX",path=path))
    
     nt,nx,ny = n.shape
     dy = (2.0*np.pi*zmax)/(ny)
     dky = 1.0/zmax

     print 'dky', dky
except:
     print "fail"

meta=[]
#blob_db = blob_info(n)

print dx.shape,dy
meta={'y0':0.0,'x0':0.0,'dx':dx,'dy':dy}

print 'movie'
print n.shape

#blob = BlobDB(n,meta=meta)
#blob = BlobMovie(n,meta=meta)
blob = Blob2D(n,meta=meta)

time = np.squeeze(collect("t_array",path=path,xind=[0,0]))
#how to quickly make a movie
from frame import Frame, FrameMovie

##pp = PdfPages('n0.pdf')
fig = plt.figure()
frm_n0 = Frame(n0,meta={'stationary':True,'dx':dx})
#frm_n0.render(fig,111)
frm_n0.ax = fig.add_subplot(111)
frm_n0.img = frm_n0.ax.plot(frm_n0[:,0])
#fig.savefig(pp,format='pdf')
fig.savefig('n0.eps')
plt.close(fig)
##pp.close()

data_c = phi
frm_data = Frame(blob.raw_data,meta={'mask':True,'dx':dx,'dy':dy,'title':'n'})
frm_amp = Frame(blob.amp,meta={'title':r'$n_{max}$'})


frm_fft = Frame((blob.power[:,0:30,0:60]),meta={'mask':True,'title':'power',
                                               'ylabel':r'$k_y \rho_s$',
                                              'xlabel':r'$k_x$'})


allk = dky*np.arange(ny)+(1e-8*dky)
mu = 1.0e-3
alpha = 3.0e-5
beta = 6.0e-4
Ln = 130.0/4.0 
n0 = 10.0
ii = complex(0,1)
soln = {}
soln['freq'] = []
soln['gamma'] = []
soln['gammamax'] = []
soln['freqmax'] = []

for i,k in enumerate(allk):
     M = np.zeros([2,2],dtype=complex)
     #density
     M[0,0] = -ii*mu*(k**2)
     M[0,1] = k*n0/Ln
     #potential
     M[1,0] = -beta/(n0*k)
     M[1,1] = -ii*(alpha + mu*k**4)/(k**2)
     #M = M.transpose()
     eigsys= np.linalg.eig(M)  
     gamma = (eigsys)[0].imag
     omega =(eigsys)[0].real
     eigvec = eigsys[1]
     #print 'k: ',k
     
     soln['gamma'].append(gamma)
     soln['gammamax'].append(max(gamma))
     where = ((gamma == gamma.max()))
     soln['freqmax'].append(omega[where])
     soln['freq'].append(omega)
     
     #print max(gamma)
gamma = (np.gradient(np.log(np.real(np.sqrt(blob.power))))[0])/(np.gradient(time)[0])

#gamma is nt x nx x ny
gamma = np.average(gamma[-10:-1,:,:],axis=0)

# gamma = np.sum(np.sqrt(blob.power[1:-1,:,:])*gamma[1:-1,:,:],axis=0)
# gamma = gamma/np.sum(np.sqrt(blob.power[1:-1,:,:]),axis=0)

# gamma = Frame(gamma[:,0,1:60],meta={'dx':dky,'xlabel':r'$k_y$',
#                                      'title':r'$\gamma$',
#                                      'ylabel':r'$\frac{\omega}{\omega_{ci}}$',
#                                      'overplot':soln['gammamax'][1:60],
#                                      'x0':dky})

#gamma_th = np.duplicate(soln['gammamax'])
#gamma_th = np.resize(gamma_th,ny)'overplot':0*np.arange(ny/10-1),

gamma_exp = Frame(gamma[0,1:ny/10:4],meta={'dx':4*dky,'xlabel':r'$k_y \rho_s$',
                                          'title':r'$\gamma_{linear}$',
                                          'ylabel':r'$\frac{\gamma}{ \omega_{ci}}$',
                                          'x0':dky,'shareax':False,'style':'ro',
                                          'ticksize':28,'stationary':True,
                                          'markersize':8,'linewidth':2,'style':'o',
                                          'overplot':0*np.arange(ny/8-1)})

gamma_th = Frame(np.array(soln['gammamax'][1:ny/10]),meta={'dx':dky,'x0':dky,'stationary':True,'nt':gamma_exp.nt,'yscale':'linear','title':r'$\gamma_{linear}$','fontsz':28,'ylabel':r'$\frac{\gamma}{ \omega_{ci}}$','xlabel':r'$k_y \rho_s$','ticksize':28,'style':'--','linewidth':15})

#let's create single frame
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
lin_formatter = ticker.ScalarFormatter()
from pylab import legend
lin_formatter.set_powerlimits((1, 1))
#plt.autoscale(axis='x',tight=True)
#self.ax.axis('tight')
#pp = PdfPages('gamma.pdf')
fig = plt.figure()
gamma_exp.ax = None
gamma_exp.t = 100
gamma_th.ax = None
gamma_exp.render(fig,111)
gamma_th.render(fig,111)
plt.tick_params(axis='both',direction='in',which='both',labelsize=20)
print dir( gamma_th.ax.yaxis.get_offset_text())
gamma_th.ax.yaxis.get_offset_text().set_size(20)
#exit()

gamma_exp.ax.yaxis.set_major_formatter(lin_formatter) 
plt.setp(gamma_th.img, color='b', linewidth=5.0,alpha=.7)
plt.setp(gamma_exp.img, color='r', linewidth=5.0,alpha=.7)
gamma_exp.ax.xaxis.set_label_coords(.65, -0.05)
plt.autoscale(axis='x',tight=True)
print 'gamma.img: ',gamma_exp.img
leg = plt.legend([gamma_exp.img,gamma_th.img],
                 ('BOUT++', 'analytic'),
                 'best', shadow=False, fancybox=True,
                 fontsize = 20)
leg.get_frame().set_alpha(0.6)
plt.tight_layout()
#fig.savefig(pp,format='pdf')
fig.savefig('gamma.eps')
plt.close(fig)

# # fig = plt.figure()


# # fig.savefig(pp,format='pdf')
# # plt.close(fig)

# #pp.close()

# #fancy overplot tricks - create the figure outside of MovieFrame
# fig = plt.figure()
# gamma_exp.t = 0
# gamma_exp.ax = fig.add_subplot(111)
# gamma_th.ax = gamma_exp.ax
# #gamma_th.ax = fig.add_subplot(111)
# #FrameMovie([gamma,gamma_th],fast=True,moviename='gamma',fps=6,fig=fig)




# #reset 
# gamma_exp.ax = None
# gamma_th.ax = None
# #gamma.t = 0
# #print gamma[0,:]
# #/elem['dt']]
# #fig.savefig(pp,format='pdf')
# fig.savefig('gamma.eps')
# #pp.close()
# fig = plt.figure()
# # gamma_exp.ax = fig.add_subplot(223)
# # gamma_th = Frame(np.array(soln['gammamax'][1:ny/3]),meta={'dx':dky,'x0':dky,'stationary':True,'nt':gamma_exp.nt,'yscale':'symlog'})
# # gamma_th.ax =fig.add_subplot(223)
# # gamma_exp.shareax = True
# # gamma_exp.yscale = 'symlog'

# print "power", blob.power.shape 
# sigma = blob.raw_data.std(axis=2)
# frm_data1D = Frame(np.log(np.sqrt(np.amax(blob.power[:,0,:],axis=1))))#,meta={'sigma':sigma,'t_array':time})
# frm_data1D.render(fig,111)
# frames= [frm_data,frm_fft,gamma_th,gamma,frm_amp]
# #FrameMovie(frames,fast=False,fig=fig)

# #fig.savefig(pp,format='pdf')
# fig.savefig('gamma.eps')
# pp.close()

# print blob.power.shape
#fig = plt.figure()


# #let's examine the boundary
# bc_left = Frame(blob.raw_data[:,0:20,:],meta={'mask':True})
# bc_right = Frame(blob.raw_data[:,-5:-1,:],meta={'mask':True})

# sig_lft = blob.raw_data[:,0:20,:].std(axis=2)
# bc1d_left = Frame(np.average(blob.raw_data[:,0:20,:],axis=2),meta={'sigma':sig_lft,'t_array':time})
# FrameMovie([bc_left,bc1d_left],fast=False,moviename='boundary',fps=2)
