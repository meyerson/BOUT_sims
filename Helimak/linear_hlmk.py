#!/opt/apps/python/epd/7.2.2/bin/python
import sys, os, gc
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




from read_inp import read_inp, parse_inp
import sys,os,inspect,shutil,subprocess
import argparse
import multiprocessing
from scipy.optimize import curve_fit ,fmin#,minimize
from scipy import signal
from frame import Frame, FrameMovie

#from turb_info2 import field_info as sim
cmd_folder = HOME+'/BOUT_sims/blob_py'

print cmd_folder

if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

parser = argparse.ArgumentParser()

parser.add_argument("path", type=str,
                    help="data location")

parser.add_argument("key", type=str,
                    help="some ID string")

parser.add_argument("tstart", type=int,
                    help="start time",nargs='?', default=0)

parser.add_argument("tstop", type=int,
                    help="stop time",nargs='?', default=100)
 
parser.add_argument("tchunk", type=int,
                    help="time chunk",nargs='?', default=10)

parser.add_argument('--debug',dest='debug',action='store_true')
parser.set_defaults(debug=False)

args = parser.parse_args()

path = args.path
key = args.key
tstart = args.tstart
tstop = args.tstop
tchunk = args.tchunk
debug = args.debug

save_path = path.replace('scratch','work')+'movie'

import matplotlib
# if debug:
#      matplotlib.use('pdf')
#      from matplotlib.backends.backend_pdf import PdfPages
# else:
#      #matplotlib.use('Agg')

matplotlib.use('pdf')
from matplotlib.backends.backend_pdf import PdfPages

#from matplotlib.backends.backend_pf import PdfPages
import matplotlib.pyplot as plt

if os.path.exists(save_path):
     #shutil.rmtree(save_path)
     for root, dirs,files in os.walk(save_path):
          for f in files:
               os.unlink(os.path.join(root, f))
else:
     os.makedirs(save_path)

#key='movie'
cache=path+'_movie'

nfile = path+ 'nfile.dat'
ufile = path+ 'ufile.dat'
phifile = path+ 'phifile.dat'
Akfiile = path+'Akfile.dat'
Tefile =  path+'Tefile.dat'
from boutdata import collect
from boutdata import collect2 
from collect2 import collect2 as collect

from boututils import savemovie

import numpy as np

inp = read_inp(path=path)
meta = parse_inp(inp)
#print meta
nz = np.double(meta['[main]']['MZ'])
nxpe = np.double(meta['[main]']['NXPE'])
nype = np.double(meta['[main]']['NYPE'])
mxg = np.double(meta['[main]']['MXG'])
dx =  np.double(meta['[mesh]']['dx'])
zmax = np.double(meta['[main]']['ZMAX'])
#mxsub = meta['[main]']['']
#nz = np.squeeze(collect("MZ",xind=[0,0],path=path,info=False))
#nxpe=  np.squeeze(collect("NXPE",xind=[0,0],path=path,info=False))
mxsub = 15
#mxsub = np.squeeze(collect("MXSUB",xind=[0,0],path=path,info=True)) #without gaurds
#mxg = np.squeeze(collect("MXG",xind=[0,0],path=path,info=False))

nx = mxsub*nxpe + 2*mxg
ny = nz
#dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
#dy = np.squeeze(collect("dz",path=path,xind=[0,0]))
#zmax = np.squeeze(collect("ZMAX",path=path))
dy = zmax/nz
yO = -.5*(dy*ny)
xO = 0.0
#dy = (2.0*np.pi*zmax)/(ny)
dky = 1.0/zmax
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))#[tstart:tstop+1]
dt = time[1]-time[0]
# print time
# exit()
a = np.squeeze(collect("alpha",path=path,info=False))
a_smooth =a
# a_smooth = np.squeeze(collect("alpha_smooth",path=path))
# mask = np.squeeze(collect("alpha_mask",path=path))
#beta = 5.0e-4
beta = np.double(meta['[physics]']['beta'])
B0 = np.squeeze(collect("B0",path=path,info=False))
x0=0
y0=0

ymin =y0
xmin =x0
xmax =nx*dx + xmin
ymax =ny*dy + ymin
pos = np.mgrid[xmin:xmax:dx,ymin:ymax:dy]

def get_data(start,stop):
     
     #n = np.exp(np.squeeze(collect("n",tind=[start,stop],path=path,info=False)))
     n = (np.squeeze(collect("n",tind=[start,stop],path=path,info=False)))
     #phi = (np.squeeze(collect("phi",tind=[start,stop],path=path,info=False)))
     n_mmap = np.memmap(nfile,dtype=n.dtype.name,mode='w+',shape=n.shape)
     n_mmap[:] = n[:]
     
     print 'n_mmap.shape :',n_mmap.shape,n.shape
     del n
     gc.collect()
     u = np.squeeze(collect("u",tind=[start,stop],path=path,info=False))
     u_mmap = np.memmap(ufile,dtype=u.dtype.name,mode='w+',shape=u.shape)
     u_mmap[:] = u[:]
     del u

     gc.collect()
     #A_k = np.fft.rfft(u_mmap,axis=2)
     A_k = np.fft.rfft(n_mmap,axis=2)
     #power = fft_u.conj()*fft_u
     #A_k = np.real(np.sqrt(power))
     
     #del fft_u,power
     
     gc.collect()
     phi = np.squeeze(collect("phi",tind=[start,stop],path=path,info=False))
     phi_mmap = np.memmap(phifile,dtype=phi.dtype.name,mode='w+',shape=phi.shape)
     phi_mmap[:] = phi[:]
     
     del phi


     gc.collect()
     
     T = (np.squeeze(collect("Te",tind=[start,stop],path=path,info=False)))
     T_mmap = np.memmap(Tefile,dtype=T.dtype.name,mode='w+',shape=T.shape)
     T_mmap[:] = T[:]
    
     del T
     gc.collect()
     return n_mmap,u_mmap,A_k,phi_mmap,T_mmap

def expfall(x,y0,l):
     return y0*np.exp(-x/l) 

def atan(params,*args):
     ydata = args[0]
     x = args[1]

     y0 = params[0]
     l = params[1]
     start = params[2]

     return y0*np.arttan(-x/l)+np.pi/2

def linearfall(x,y0,l):
    #print y0,l,np.min(x)
    return np.log(y0) - x/l

def gauss_kern(size, sizey=None):
     """ Returns a normalized 2D gauss kernel array for convolutions """
     size = int(size)
     if not sizey:
          sizey = size
     else:
          sizey = int(sizey)
     x, y = np.mgrid[-5:5+1, -sizey:sizey+1]
     g = np.exp(-(x**2/float(5)+y**2/float(sizey)))
     return g / g.sum()


def expfall2(params,*args):
     #print 'args', len(args),args.__class__,len(params)
     #print 'params', params

     ydata = args[0]
     x = args[1]
     mx = args[2]


     y0 = params[0]
     l = params[1]
     start = params[2]
     
    # print mx*start
     #popt, pcov= fit_lambda(x[start:-1],x[xstart:-1],p0=p0)
     #print start,ydata.shape,x.shape
     #print ydata
     nnx = len(x[int(start*mx):0])
     
     if start < 1 and start > 0:
          return np.sum((y0*np.exp(-x[np.int(start*mx):0]/l) -ydata[np.int(start*mx):0])**2)/nnx**2
     else:
          return 1e10

def fit_lambda(y,x,appendto=None,p0=None):
     #try:
     
     return curve_fit(expfall, x.flatten(), y.flatten(),p0=p0)
     #except:
     #     return np.zeros(2),np.zeros([2,2])
def fit_lambda2(y,x,mx,p0):
     print p0, type(p0),type(p0[0])
     #x0 = np.asfarray(p0.flatten())

     res = fmin(expfall2,p0,args=(y,x,mx))
     return res
     
     #return curve_fit(expfall2, x.flatten(), y.flatten(),p0=p0)

t1 = tstart
t2 = t1+tchunk

if debug:
     print debug
     n,u,Ak,phi,Te = get_data(t1,t2)
     pp = PdfPages('debug.pdf')
     fig = plt.figure()
     print a_smooth.shape,path
     vy = -((np.gradient(phi)[2])/dy)
     debug_frm = Frame(np.mean(np.mean(vy,axis=2),axis=0),meta={'dx':dx,
                                                                   'stationary':True,'xlabel':'x['+r'$\rho_s$'+']',
                                                                   'ylabel':r'$v_y[C_s]$','yscale2':1,'ylabel2':''})
     debug_frm.render(fig,111)
     fig.savefig(pp, format='pdf')
     pp.close()
     exit()
print 'enter main loop'
while t2<=tstop:
     
     print 'about to get data'
     n,u,Ak,phi,Te = get_data(t1,t2)
     
     print 'shape:', Ak.shape



     nt,nx,ny = n.shape
     time = np.squeeze(collect("t_array",path=path,xind=[0,0]))[t1:t2+1]
     nx_sol = np.round(.4*nx) 
          
     #-.17949 *(dx*nx)
     data_c = phi
     
     power = (Ak.conj()*Ak)

     gamma = np.gradient(np.sqrt(power))[0]/dt
     w = np.gradient(np.angle(Ak))[0]/dt
     
     print 'gamma: ', gamma.shape
     t1 = t1+tchunk
     t2 = t2+tchunk

     # exit()
 

pp = PdfPages('linear_hlmk.pdf')  
debug_frm = Frame(np.mean(np.mean(gamma[:,-nx/4:-nx/4+10,0:30],axis=0),axis=0),
                  meta={'dx':dky,
                        'stationary':True,'xlabel':r'$k_y \rho_s$',
                        'ylabel':r'$\frac{\gamma}{\omega_{ci}}$','yscale2':1,'ylabel2':''})
fig = plt.figure() 
debug_frm.render(fig,111)
fig.savefig(pp,format='pdf')

debug_frm = Frame(np.mean(np.mean(w[:,-nx/4:-nx/4+10,0:30],axis=0),axis=0),
                  meta={'dx':dky,
                        'stationary':True,'xlabel':r'$k_y \rho_s$',
                        'ylabel':r'$\frac{\gamma}{\omega_{ci}}$','yscale2':1,'ylabel2':''})
fig = plt.figure() 
debug_frm.render(fig,111)
fig.savefig(pp,format='pdf')

debug_frm = Frame(np.mean(np.mean((gamma/w)[:,-nx/4:-nx/4+10,0:30],axis=0),axis=0),
                  meta={'dx':dky,
                        'stationary':True,'xlabel':r'$k_y \rho_s$',
                        'ylabel':r'$\frac{\gamma}{\omega_{ci}}$','yscale2':1,'ylabel2':''})
fig = plt.figure() 
debug_frm.render(fig,111)
fig.savefig(pp,format='pdf')

def analytic(ny,dky,mu = 1.0e-2,alpha = 2.0e-3,beta = 1.0e-2,Ln = 30.0,
             L_phi = 2*30.0,L_t = 2*30.0,Lambda = 5.0):
     
     allk = .1*dky*np.arange(10*ny)+(1e-8*dky)
     # mu = 1.0e-2
     # alpha = 2.0e-3
     # beta = 1.0e-2
     # Ln = 30.0
     # L_phi = 2*Ln
     # L_t = 2*Ln
     # Lambda = 5.0
     n0 = 10.0
     ii = complex(0,1)
     soln = {}
     soln['freq'] = []
     soln['gamma'] = []
     soln['gammamax'] = []
     soln['freqmax'] = []


     #k_r = 0 for now
     for i,k in enumerate(allk):
          M = np.zeros([3,3],dtype=complex)
          #density
          M[0,0] = -ii*mu*(k**2)-ii*alpha-Lambda*k/L_phi - beta*k
          M[0,1] = k/Ln+ ii*alpha + beta*k
          M[0,2] = -ii* alpha/2 - ii*alpha*Lambda*-beta*k
          #potential
          M[1,0] = -beta/k
          M[1,1] = -ii*(alpha + mu*k**4)/(k**2) 
          M[1,1] = M[1,1] + Lambda*(1.0/(k*L_phi**3) + k/L_phi)
          M[1,2] = ii*alpha*Lambda/k**2 - beta/k
          #temprature
          M[2,0] = -2.0*beta*k/3.0
          M[2,1] = ii*57*alpha/50 + 2.0*beta*k/3.0 + k/L_t
          M[2,2] = -ii*alpha - (57/50)*ii*alpha*Lambda - (7.0/3.0)*beta*k 
          M[2,2] = M[2,2] - ii*mu*k**2 - Lambda*k/L_phi

          #M = M.transpose()
          eigsys= np.linalg.eig(M)  
          gamma = (eigsys)[0].imag
          omega =(eigsys)[0].real
          eigvec = eigsys[1]
          #print 'k: ',k
          
          soln['gamma'].append(gamma)
          soln['gammamax'].append(max(gamma))
          where = ((gamma == gamma.max()))
          #print where, omega,omega[where],gamma
          #exit()
          soln['freqmax'].append(omega[where])
          soln['freq'].append(omega)
     return soln

soln  = analytic(ny,dky)
gamma_th = Frame(np.array(soln['gammamax'][0:ny/6]),meta={'dx':.1*dky,'x0':0,'stationary':True,'yscale':'symlog','title':r'$\gamma$','fontsz':18,'ylabel':r'$\frac{\omega}{\omega_{ci}}$','xlabel':r'$k_y$','ticksize':14})
fig = plt.figure() 
gamma_th.render(fig,111)
fig.savefig(pp,format='pdf')

# print np.array(soln['freqmax'][1:ny/6]).shape, np.array(soln['gammamax'][1:ny/6]).shape
# exit()
omega_th= Frame(np.squeeze(np.array(soln['freqmax'][1:ny/6])),meta={'dx':.1*dky,'x0':.1*dky,'stationary':True,'yscale':'symlog','title':r'$\omega$','fontsz':18,'ylabel':r'$\frac{\omega}{\omega_{ci}}$','xlabel':r'$k_y$','ticksize':14})
fig = plt.figure() 
omega_th.render(fig,111)
fig.savefig(pp,format='pdf')


soln  = analytic(ny,dky,Ln = 20,L_phi = -40,L_t = -40)
gamma_th = Frame(np.array(soln['gammamax'][0:ny/6]),meta={'dx':.1*dky,'x0':0,'stationary':True,'yscale':'symlog','title':r'$\gamma$','fontsz':18,'ylabel':r'$\frac{\omega}{\omega_{ci}}$','xlabel':r'$k_y$','ticksize':14})
fig = plt.figure() 
gamma_th.render(fig,111)
fig.savefig(pp,format='pdf')

# print np.array(soln['freqmax'][1:ny/6]).shape, np.array(soln['gammamax'][1:ny/6]).shape
# exit()
omega_th= Frame(np.squeeze(np.array(soln['freqmax'][1:ny/6])),meta={'dx':.1*dky,'x0':.1*dky,'stationary':True,'yscale':'symlog','title':r'$\omega$','fontsz':18,'ylabel':r'$\frac{\omega}{\omega_{ci}}$','xlabel':r'$k_y$','ticksize':14})
fig = plt.figure() 
omega_th.render(fig,111)
fig.savefig(pp,format='pdf')


#soln.plot

pp.close()
exit()


