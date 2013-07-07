#!/opt/apps/python/epd/7.2.2/bin/python
import os, sys, inspect
import sqlite3 as sql
import pickle as pkl
from datetime import datetime
import shutil


HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')
utc = datetime.utcnow()

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from cStringIO import StringIO
import sys

old_stdout = sys.stdout

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

# from blob_info import blob_info, Blob2D
# from blob_draw import BlobDraw
from turb_info import Turbulence as Trblnc2D


import numpy as np
from scipy.optimize import curve_fit

from boutdata import collect

from frame import Frame, FrameMovie

defaults ={} 
defaults['path'] = '/scratch/01523/meyerson/BOUT_sims/TurbulenceSOL/data_turb_mu1e-2_bcfalse_false'
defaults['tchunk'] =10
defaults['key'] = 'lamdahist'


arg_names = ['command', 'path', 'tchunk', 'key']
args = dict(zip(arg_names, sys.argv))
for x in defaults.keys():
     if x not in args:
          args[x] = defaults[x]
     
path = args['path']
tchunk = np.int(args['tchunk'])
key = args['key']     

save_path = path.replace('scratch','work')

if os.path.exists(save_path):
     #shutil.rmtree(save_path)
     for root, dirs,files in os.walk(save_path):
          for f in files:
               os.unlink(os.path.join(root, f))
               
# al =  [if x not in args else  for x in defaults.keys()]
# print defaults,va

#read the data and process

sys.stdout = mystdout = StringIO()
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))
dx = np.squeeze(collect("dx",path=path))
zmax = np.squeeze(collect("ZMAX",path=path))
#xpiece = np.squeeze(collect("x",path=path))
#ypiece = np.squeeze(collect("y",path=path))
NXPE = np.squeeze(collect("NXPE",path=path))
n = np.squeeze(collect("n",path=path,tind =[0,0]))
nx,ny = n.shape
n = np.squeeze(collect("n",path=path,xind =[0,0],zind=[0,0],yind=[0,0]))
nt = len(n)
dy = (2.0*np.pi*zmax)/ny
dx = np.mean(dx)
sys.stdout = old_stdout

x0=0
y0=0

ymin =y0
xmin =x0
xmax =nx*dx + xmin
ymax =ny*dy + ymin
print ymax,xmax,ymin,xmin,dx,dy


pos = np.mgrid[xmin:xmax:dx,ymin:ymax:dy]
#pos = np.mgrid[xmin:xmax-dx:nx*complex(0,1),ymin:ymax-dy:ny*complex(0,1)]
#pos_i = np.mgrid[0:nx-1:nx*complex(0,1),0:ny-1:ny*complex(0,1)]
                


def expfall(x,y0,l):
     return y0*np.exp(-x/l)    

def fit_lambda(y,x,appendto=None,p0=None):
     try:
          return curve_fit(expfall, x.flatten(), y.flatten(),p0=p0)
     except:
          return np.zeros(2),np.zeros([2,2])
print nx,ny,nt

stop = True
t = 0
lam = []
lam_rough =[]
lambdafit = []
offset = []
edgecut = np.round(nx/10.)
xstart = np.int(edgecut)
xstop= np.int(nx - edgecut)
time  = []
nmax  =[]
print xstart,xstop
while  t < np.round(.99*nt):
     sys.stdout = mystdout = StringIO()
     n = np.squeeze(collect("n",path=path,xind=[xstart,xstop-1],tind =[t,t]))
     phi = np.squeeze(collect("phi",path=path,xind=[xstart,xstop-1],tind =[t,t]))
     t_len = n.shape[0]
     sys.stdout = old_stdout
     time.append(t)
     t=t+tchunk
     

     print 't: ',t
     
     #for i in xrange(t_len):
     nave = np.mean(n,axis = 1)
          
     print pos.shape,pos[0][xstart:xstop,:].shape,n.shape,nave[-1]
     #popt, pcov= fit_lambda(n-nave[-1],pos[0][xstart:xstop,:])
     
     est_lam = (pos[0][xstop,5]-pos[0][xstart,5])/(np.log(nave[0]/nave[-1]))
     p0=[nave[0]-nave[-1],est_lam]
     popt, pcov= fit_lambda(n - nave[-1],pos[0][xstart:xstop,:],p0=p0)
     #est_lam = (pos[0][xstop,5]-pos[0][xstart,5])/(np.log(nave[0]/nave[-1]))
     #popt, pcov= curve_fit(expfall,pos[0][xstart:xstop,:],n-nave[-1],p0=[nave[0]-nave[-1],est_lam])
     print est_lam,popt,pcov,nave[0]-nave[-1]
     lam.append(popt[1])
     lam_rough.append(est_lam)
     lambdafit.append(popt[0]*np.exp(-pos[0][xstart:xstop,:]/popt[1]))
     offset.append(popt[0])
     nmax.append(n.max())

     print n.shape,phi.shape
     print 'time:', t


     sys.stdout = mystdout = StringIO()
     pp = PdfPages(save_path+'/'+key+'lam.pdf')
     fig = plt.figure() 
     lam_history = Frame(lam,meta={'dx':tchunk,'stationary':False,'fontsz':18,'ylabel':'',
                                   'xlabel':r'$t$','ticksize':14,'title':r'$\lambda$','xlabel':r't'}) 
     
     lam_history.ax = fig.add_subplot(111)
     lam_history.ax.set_ylim([0.0,100.*np.round((lam_history[-1]+100.)/100.)])
     lam_history.render(fig,111)
     
     print lam_history.x
     fig.savefig(pp,format='pdf')
     plt.close(fig)
     pp.close()

     pp = PdfPages(save_path+'/'+key+'compare_lam.pdf')
     fig = plt.figure()
     lam_history.ax = None
     lam_history.yscale='linear'
     
     #lam_history = Frame(lam,meta={'stationary':False,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$',
     #                              'ticksize':14,'yscale':'linear'})
     lam_rough_history = Frame(lam_rough,meta={'stationary':False,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$',
                                               'ticksize':14,'yscale':'linear','dx':tchunk})
     lam_history.render(fig,111)
     lam_rough_history.render(fig,111)
     
     fig.savefig(pp,format='pdf')
     plt.close(fig)
     pp.close()
                                        

     pp = PdfPages(save_path+'/'+key+'loglam.pdf')
     fig = plt.figure()
     lam_history.ax = None
     lam_history = Frame(lam,meta={'stationary':False,'dx':tchunk,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14,'yscale':'symlog'}) 
     lam_history.render(fig,111)
     fig.savefig(pp,format='pdf')
     plt.close(fig)
     pp.close()

     pp = PdfPages(save_path+'/'+key+'lamda_summary.pdf')
     fig = plt.figure() 
     lam_history.ax = None
     #lam_history = Frame(lam,meta={'stationary':False,'dx':tchunk'title':r'$\lambda_n$','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14,'yscale':'symlog'}) 
     lam_history.render(fig,211)
     nmax_history = Frame(nmax,meta={'stationary':False,'dx':tchunk,'title':r'$n_{max}$','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14,'yscale':'linear'}) 
     nmax_history.render(fig,212)
     fig.savefig(pp,format='pdf')
     plt.close(fig)
     pp.close()
     
     sys.stdout = old_stdout

     
     


n = n[:,:,:]
phi = phi[1:,:,:]

dy = (2.0*np.pi*zmax)/(ny)
dky = 1.0/zmax
print 'dky', dky
     
meta=[]
#blob_db = blob_info(n)

print dx.shape,dy,dx
meta={'y0':0.0,'x0':0.0,'dx':dx,'dy':dy}

print 'movie'
print n.shape

#blob = BlobDB(n,meta=meta)
#blob = BlobMovie(n,meta=meta)
#blob = Blob2D(n,meta=meta)
blob = Trblnc2D(n,meta=meta)
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))
#how to quickly make a movie




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
pp = PdfPages('lamda.pdf')
fig = plt.figure() 
lam_history = Frame(lam,meta={'stationary':False,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14}) 
lam_history.render(fig,111)
fig.savefig(pp,format='pdf')
plt.close(fig)
pp.close()


pp = PdfPages('gamma.pdf')
fig = plt.figure()
gamma.ax = None
gamma.t = 1
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

print blob.lam
#lam_history = Frame(np.array(blob.lam),meta={'stationary':False,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14})
lam_history = Frame(lam,meta={'stationary':False,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14})

lam_history.render(fig,223)

n_fit = Frame(blob.lambdafit[-1][:,ny/2],meta={'dx':dx,'x0':0,'stationary':True,'title':'','fontsz':18,'ylabel':'','xlabel':r'$t$','ticksize':14})

n_fit.render(fig,224)
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
