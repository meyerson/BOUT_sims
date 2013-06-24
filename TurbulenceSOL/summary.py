#!/usr/bin/python2.7
import os, sys, inspect
import sqlite3 as sql
import pickle as pkl
from datetime import datetime

HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')
utc = datetime.utcnow()

from matplotlib.backends.backend_pdf import PdfPages
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
from turb_info import Turbulence as Trblnc2D


import numpy as np
from boutdata import collect

#prepare a list of directories
sim_key='Ra1e4_turb'
path="/tmp/SOLblob/data_"+sim_key



#read the data and process

n = np.squeeze(collect("n",path=path,tind =[1,550]))
#n0 = np.squeeze(collect("n0",path=path,tind =[1,550]))
u = np.squeeze(collect("u",path=path,tind =[1,550]))
phi = np.squeeze(collect("phi",path=path,tind =[1,550]))
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))
dx = np.squeeze(collect("dx",path=path))
zmax = np.squeeze(collect("ZMAX",path=path))
nt,nx,ny = n.shape
dy = (2.0*np.pi*zmax)/(ny)
dky = 1.0/zmax
print 'dky', dky
     
meta=[]
#blob_db = blob_info(n)

print dx.shape,dy
meta={'y0':0.0,'x0':0.0,'dx':dx[0],'dy':dy}

print 'movie'
print n.shape

#blob = BlobDB(n,meta=meta)
#blob = BlobMovie(n,meta=meta)
#blob = Blob2D(n,meta=meta)
blob = Trblnc2D(n,meta=meta)
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))
#how to quickly make a movie
from frame import Frame, FrameMovie

# pp = PdfPages('n0.pdf')
# fig = plt.figure()
# frm_n0 = Frame(n0,meta={'stationary':True,'dx':dx})
# #frm_n0.render(fig,111)
# frm_n0.ax = fig.add_subplot(111)
# frm_n0.img = frm_n0.ax.plot(frm_n0[:,0])
# fig.savefig(pp,format='pdf')
# plt.close(fig)
# pp.close()

data_c = phi
frm_data = Frame(blob.raw_data,meta={'mask':True,'dx':dx,'dy':dy,'title':'n'})
frm_amp = Frame(blob.amp,meta={'title':r'$n_{max}$'})


frm_fft = Frame((blob.power[:,0:30,0:60]),meta={'mask':True,'title':'power',
                                               'ylabel':r'$k_y$',
                                              'xlabel':r'$k_x$'})


allk = dky*np.arange(ny)+(1e-8*dky)
mu = 1.0e-2
alpha = 2.0e-5
beta = 6.0e-4
Ln = 25.5/4.0 
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

# gamma = Frame(gamma[:,0,1:60],meta={'dx':dky,'xlabel':r'$k_y$',
#                                      'title':r'$\gamma$',
#                                      'ylabel':r'$\frac{\omega}{\omega_{ci}}$',
#                                      'overplot':soln['gammamax'][1:60],
#                                      'x0':dky})

#gamma_th = np.duplicate(soln['gammamax'])
#gamma_th = np.resize(gamma_th,ny)
gamma = Frame(gamma[:,0,1:ny/6],meta={'dx':dky,'xlabel':r'$k_y$',
                          'title':r'$\gamma$',
                          'ylabel':r'$\frac{\omega}{\omega_{ci}}$',
                          'overplot':0*np.arange(ny/6-1),
                          'x0':dky,'shareax':False,'style':'ro',
                                      'ticksize':14})

gamma_th = Frame(np.array(soln['gammamax'][1:ny/6]),meta={'dx':dky,'x0':dky,'stationary':True,'nt':gamma.nt,'yscale':'symlog','title':r'$\gamma$','fontsz':18,'ylabel':r'$\frac{\omega}{\omega_{ci}}$','xlabel':r'$k_y$','ticksize':14})


import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
lin_formatter = ticker.ScalarFormatter()
from pylab import legend
lin_formatter.set_powerlimits((1, 1))
#plt.autoscale(axis='x',tight=True)
#self.ax.axis('tight')

#let's create single frame
pp = PdfPages('gamma.pdf')
fig = plt.figure()
gamma.ax = None
gamma.t = 20
gamma_th.ax = None
gamma.render(fig,111)
gamma_th.render(fig,111)

gamma.ax.yaxis.set_major_formatter(lin_formatter) 
plt.setp(gamma_th.img, color='b', linewidth=3.0,alpha=.7)
plt.setp(gamma.img, color='r', linewidth=2.0,alpha=.7)
plt.autoscale(axis='x',tight=True)
print 'gamma.img: ',gamma.img
leg = plt.legend([gamma.img,gamma_th.img],('BOUT++', 'analytic'),
                 'best', shadow=False, fancybox=True)
leg.get_frame().set_alpha(0.6)
fig.savefig(pp,format='pdf')
plt.close(fig)
pp.close()


#let's create a multipaneled summary of the results
pp = PdfPages('summary.pdf')
fig = plt.figure()
frm_amp.render(fig,221)

print 'shape: ', blob.lambdafit[0].shape

n_last = Frame(np.array(n[-1,:,ny/2]),meta={'dx':dx,'x0':0,'stationary':True,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14, 'overplot':blob.lambdafit[-1][:,ny/2]})
n_last.render(fig,222)

lam_history = Frame(np.array(blob.lam),meta={'stationary':False,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14})
lam_history.render(fig,223)

fig.savefig(pp,format='pdf')
plt.close(fig)
pp.close()

#fancy overplot tricks - create the figure outside of MovieFrame
fig = plt.figure()
gamma.t = 0
gamma.ax = fig.add_subplot(111)
gamma_th.ax = gamma.ax
#gamma_th.ax = fig.add_subplot(111)
#FrameMovie([gamma,gamma_th],fast=True,moviename='gamma',fps=6,fig=fig)




#reset 
gamma.ax = None
gamma_th.ax = None
#gamma.t = 0
#print gamma[0,:]
#/elem['dt']]
fig = plt.figure()
gamma.ax = fig.add_subplot(223)
gamma_th = Frame(np.array(soln['gammamax'][1:ny/3]),meta={'dx':dky,'x0':dky,'stationary':True,'nt':gamma.nt,'yscale':'symlog'})
gamma_th.ax =fig.add_subplot(223)
gamma.shareax = True
gamma.yscale = 'symlog'


sigma = blob.raw_data.std(axis=2)
frm_data1D = Frame(np.average(blob.raw_data,axis=2),meta={'sigma':sigma,'t_array':time})
frames= [frm_data,frm_fft,gamma_th,gamma,frm_amp]
#FrameMovie(frames,fast=False,fig=fig)


fig = plt.figure()


# #let's examine the boundary
# bc_left = Frame(blob.raw_data[:,0:20,:],meta={'mask':True})
# bc_right = Frame(blob.raw_data[:,-5:-1,:],meta={'mask':True})

# sig_lft = blob.raw_data[:,0:20,:].std(axis=2)
# bc1d_left = Frame(np.average(blob.raw_data[:,0:20,:],axis=2),meta={'sigma':sig_lft,'t_array':time})
# FrameMovie([bc_left,bc1d_left],fast=False,moviename='boundary',fps=2)
