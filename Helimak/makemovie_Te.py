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


import matplotlib
matplotlib.use('Agg')
#from matplotlib.backends.backend_pf import PdfPages
import matplotlib.pyplot as plt

from read_inp import read_inp, parse_inp
import sys,os,inspect,shutil,subprocess
import argparse
import multiprocessing
from scipy.optimize import curve_fit ,fmin#,minimize
from scipy import signal
from frame import Frame, FrameMovie
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


args = parser.parse_args()

path = args.path
key = args.key
tstart = args.tstart
tstop = args.tstop
tchunk = args.tchunk

save_path = path.replace('scratch','work')+'movie'


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
from collect2 import collect

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
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))[tstart:tstop+1]
a = np.squeeze(collect("alpha",path=path,info=True))
a_smooth =a
# a_smooth = np.squeeze(collect("alpha_smooth",path=path))
# mask = np.squeeze(collect("alpha_mask",path=path))
beta = 5.0e-4

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
     fft_u = np.fft.rfft(u_mmap)
     power = fft_u.conj()*fft_u
     A_k = np.real(np.sqrt(power))
     
     del fft_u,power
     
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

while t2<=tstop:

     n,u,Ak,phi,Te = get_data(t1,t2)
     
     print n.shape

     nt,nx,ny = n.shape
     time = np.squeeze(collect("t_array",path=path,xind=[0,0]))[t1:t2+1]
     nx_sol = np.round(.4*nx) 
          
     #-.17949 *(dx*nx)
     data_c = phi
     print 'a',a,np.sum((np.array(np.isfinite(a))).astype(int)-1)
     #ddtu = np.squeeze(collect("t_array",path=path,xind=[0,0]))[t1:t2+1]

    # for i in np.arange(2.*nx/3.,nx-1):
    #      a[i,:] = a[i,0]
    # a = (np.repeat(a,nt)).reshape(nx,ny,nt)
     print a.shape
     #a = np.transpose(a,(2,0,1))

     frm_n = Frame(n,meta={'dx':dx,'dy':dy,'title':'w/out chaos','cmap':'hot',
                           'xlabel':'x['+r'$\rho_s$'+']',
                           'ylabel':'y['+r'$\rho_s$'+']',
                           'fontsz':20,'interpolation':'linear','grid':False,
                           'linewidth':1,'contour_color':'black',
                           't_array':time,'x0':dx*250.0 })

     n_DC = np.swapaxes(n,1,0)
     n_DC = n_DC.reshape(nx,nt*ny)
     print n_DC.shape
     n_std = n_DC.std(axis=1)
     n_DC  = n_DC.mean(axis=1)                        
    
     
     n_DC = np.repeat(n_DC,nt*ny)
     n_DC = n_DC.reshape(nx,nt,ny)
     n_DC = np.swapaxes(n_DC,1,0)
     
     n_std = np.repeat(n_std,nt*ny)
     n_std = n_std.reshape(nx,nt,ny)
     n_std = np.swapaxes(n_std,1,0)

     n_AC = n - n_DC
     n_AC_norm = n_AC/n_std
     
 
     frm_n_AC = Frame(n_AC,
                         meta={'dx':dx,'dy':dy,'title':'w/out chaos','cmap':'hot',
                               'xlabel':'x['+r'$\rho_s$'+']',
                               'ylabel':'y['+r'$\rho_s$'+']',
                               'fontsz':20,'interpolation':'linear','grid':False,
                               'linewidth':1,'contour_color':'black',
                               't_array':time,'x0':0})

     #blobs = np.exp(n)*(np.gradient(phi)[2])
     blobs = u
     import copy
    # sigma = blobs.std(axis=2)
     blobs_data1D = Frame(np.mean(np.mean(blobs,axis=2),axis=0),meta={'t_array':time,'dx':dx})

    # blobs = n* np.gradient(phi)[2]
     b_DC = np.swapaxes(blobs,1,0)
     b_DC = b_DC.reshape(nx,nt*ny)
     print b_DC.shape
     b_std = b_DC.std(axis=1)
     b_DC  = b_DC.mean(axis=1)                        
    
     
     b_DC = np.repeat(b_DC,nt*ny)
     b_DC = b_DC.reshape(nx,nt,ny)
     b_DC = np.swapaxes(b_DC,1,0)
     
     blobs = blobs - b_DC

     sigma = np.swapaxes(blobs,1,0)
     print sigma.shape
     sigma = sigma.reshape(nx,nt*ny)
     print sigma.shape
     sigma = sigma.std(axis=1)
     sigma = np.repeat(sigma,nt*ny)
     sigma = sigma.reshape(nx,nt,ny)
     sigma = np.swapaxes(sigma,1,0)
     print sigma.shape
     
     
     sigma = sigma + np.mean(sigma)*(sigma<(np.mean(sigma)*1e-3))
     blobs = blobs/(sigma)

     blobs = blobs * (abs(blobs) >1.5)

     #blobs = blobs * (sigma<(np.mean(sigma)*1e-2))


     frm_blob = Frame(n_AC,meta={'dx':dx,'dy':dy,'title':'blobs',
                                  'xlabel':'x['+r'$\rho_s$'+']',
                                  'ylabel':'y['+r'$\rho_s$'+']',
                                  'fontsz':20,'interpolation':'linear','grid':False,
                                  'linewidth':.5,'contour_color':'red',
                                  't_array':time,'x0':dx*250.0 })

     frm_exp_data = Frame(np.exp(n),meta={'mask':True,'dx':dx,'dy':dy,'cmap':'hot'})
     
     frm_log_data = Frame(np.log(np.abs(n)),meta={'dx':dx,'dy':dy,'cmap':'hot'})
     
     
     frm_phi_data = Frame(phi,meta={'dx':dx,'dy':dy,'cmap':'hot'})

     frm_u_data = Frame(u,meta={'dx':dx,'dy':dy,'cmap':'hot'})
     frm_du_data = Frame(np.gradient(u)[0],meta={'dx':dx,'dy':dy,'cmap':'hot'})
     
     #we can include as many overplot as we want - just grab the canvas and draw whatever
     #if you are going to make movies based on stationary include nt
    
  
     # a_contour = Frame(a,meta={'stationary':True,'dx':dx,'dy':dy,'contour_only':True,'alpha':.2,'colors':'blue','grid':False,'x0':0})
     # # alpha_contour = Frame(mask,meta={'stationary':True,'dx':dx,'dy':dy,'contour_only':True,'alpha':.1,'colors':'k'})
     # a_contour.nt = frm_n.nt

     # for t in range(frm_data.nt):
     #      phi[t,:,:]-np.mean(phi[t,:,:])

     phi_contour = Frame(phi,meta={'stationary':False,'dx':dx,'contour_only':True,'alpha':.5,'colors':'red'})
     phi_contour.nt = frm_n.nt

     #frm_data_SOL = Frame(n[:,nx_sol:-1,:],meta={'mask':True,'dx':dx,'x0':dx*nx_sol})
     #frm_data = Frame(a,meta={'data_c':a,'mask':True,'dx':dx})
     print n.shape
     amp = abs(n).max(1).max(1)   
     frm_amp = Frame(amp)

     
     dky = 1.0/zmax
     allk = dky*np.arange(ny)+(1e-8*dky)
     mu = 1.0e-2
     #alpha = 3.0e-5
     # beta = 6.0e-4
     # Ln = 130.0/4.0
     
     n0 = n[0,:,:].mean(axis=1)
     #Ln  = n0/gradient(n0)
     ii = complex(0,1)
     soln = {}
     soln['freq'] = []
     soln['gamma'] = []
     soln['gammamax'] = []
     soln['freqmax'] = []
     
     # for i,k in enumerate(allk):
     #      M = np.zeros([2,2],dtype=complex)
     #      #density
     #      M[0,0] = -ii*mu*(k**2)
     #      M[0,1] = k*n0/Ln
     #      #potential
     #      M[1,0] = -beta/(n0*k)
     #      M[1,1] = -ii*(alpha + mu*k**4)/(k**2)
     #      #M = M.transpose()
     #      eigsys= np.linalg.eig(M)  
     #      gamma = (eigsys)[0].imag
     #      omega =(eigsys)[0].real
     #      eigvec = eigsys[1]
     #      #print 'k: ',k
          
     #      soln['gamma'].append(gamma)
     #      soln['gammamax'].append(max(gamma))
     #      where = ((gamma == gamma.max()))
     #      soln['freqmax'].append(omega[where])
     #      soln['freq'].append(omega)
     
     
     a_m = np.power(beta/a_smooth[:,0]**2,.20)
     a_par = np.power(beta/a_smooth[:,0]**2,1.0/3.0)
     
     a_D = beta/(a_smooth[:,0] * mu)

     a_mu = np.power(mu/a_smooth[:,0],.25)

     a_L = np.power((beta/a_smooth[:,0]**2),1./3.)

     


     print a_m
     
     frm_Ak = Frame(Ak[:,:,0:60],meta={'dy':dky,'dx':dx,
                                        'overplot':[2.*np.pi/a_m,2.*np.pi/a_L,
                                                    2.*np.pi/a_mu,2.*np.pi/a_D]})
     #FrameMovie([[frm_data,alpha_contour]],fast=True,moviename=save_path+'/'+'n_phi'+key+str(t2),fps = 10,encoder='ffmpeg')
     #FrameMovie([frm_Ak],fast=True,moviename=save_path+'/'+'u_k_phi'+key+str(t2),fps = 10,encoder='ffmpeg')

     frm_Ak.reset()
     frm_n.reset()
     # a_contour.reset()
     # a_contour.nt = frm_n.nt
     # a_contour.dx = frm_n.dx
     sigma = n.std(axis=2)
     sigma_exp = (np.exp(n)).std(axis=2)
     frm_data1D = Frame(np.average(n,axis=2),meta={'sigma':sigma,'t_array':time,'dx':dx})

     sigma = Te.std(axis=2)
     frm_Te1D = Frame(np.average(Te,axis=2),
                      meta={'sigma':sigma,'error_color':'pink',
                            't_array':time,'dx':dx})
     frm_Te2D = Frame(Te,meta={'t_array':time,'dx':dx,'dy':dy})
     
     frm_exp_data1D = Frame(np.average(np.exp(n),axis=2),meta={'sigma':sigma_exp,'t_array':time,'dx':dx})
     frm_log_data1D = Frame(np.average(np.log(np.abs(n)),axis=2),meta={'sigma':(np.log(n)).std(axis=2),'t_array':time,'dx':dx})

     

     
     # sigma = n.std(axis=2)
     # frm_data1D = Frame(np.average(n,axis=2),meta={'sigma':sigma,'t_array':time,'dx':dx})


     sigma = phi.std(axis=2)
     phi_data1D = Frame(np.average(phi/4.7,axis=2),meta={'sigma':sigma/4.7,'t_array':time,'dx':dx})
     u_data1D = Frame(np.average(u,axis=2),meta={'sigma':sigma,'t_array':time,'dx':dx})
     

     nave  = np.average(np.average(n,axis=2),axis=0)
     a_ave = np.average(a_smooth,axis=1)
     
     

    # n_fit = popt[0]*np.exp(-pos[0][xstart:xstop,5]/popt[1])
     # n_fit = Frame(n_fit,meta={'dx':dx,'x0':pos[0][xstart,5],'stationary':True})
     print 'new_frm?'
    # new_frm = phi_data1D/4.7
    # new_frm.array_finalize(phi_data1D) #shouldn't have to do this
     
     frames= [frm_exp_data1D ,frm_n_AC,[frm_Te1D,phi_data1D],frm_Te2D]
     #frames= [frm_data1D,[frm_data,phi_contour],frm_log_data1D,frm_log_data]

     
      
     frm_n.t = 0
     # frm_Ak.t = 0
     # frm_Ak.reset()
     # frm_data.reset()
     # alpha_contour.reset()
     #FrameMovie([[frm_blob_AC,dw_contour]],fast=True,moviename=save_path+'/'+key+str(t2),fps = 10,encoder='ffmpeg')
     
     frm_n.t = 0
     frm_Ak.t = 0
     frm_Ak.reset()
     frm_n.reset()
     
     FrameMovie(frames,fast=True,moviename=save_path+'/'+key+str(t2),fps = 10,encoder='ffmpeg')
     #print time, n_fit.shape,popt,pcov,nave[0:40],popt
     
     frm_n.t = 0
     frm_Ak.t = 0
     t1 = t1+tchunk
     t2 = t2+tchunk


#ls -rt $PWD/movie*mp4 | awk '{ print "file " "'\''" $1 "'\''"}'

movienames = [key]#,'n_phi'+key,'u_k_phi'+key]

#from subprocess import call
# for name in movienames:
#      print name, save_path
#      command = ('makemovlist.sh',save_path+'/',name)
#      subprocess.check_call(command)



#ls -rt $PWD/movie*mp4 | awk '{ print "file " "'\''" $1 "'\''"}'
