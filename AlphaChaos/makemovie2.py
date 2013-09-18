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
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

from read_inp import metadata
import sys,os,inspect,shutil,subprocess
import argparse
import multiprocessing
from scipy.optimize import curve_fit ,fmin#,minimize

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

from boutdata import collect
from boututils import savemovie

import numpy as np



nz = np.squeeze(collect("MZ",xind=[0,0],path=path,info=False))
nxpe=  np.squeeze(collect("NXPE",xind=[0,0],path=path,info=False))

mxsub = np.squeeze(collect("MXSUB",xind=[0,0],path=path,info=False)) #without gaurds
mxg = np.squeeze(collect("MXG",xind=[0,0],path=path,info=False))

nx = mxsub*nxpe + 2*mxg
ny = nz
dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
dy = np.squeeze(collect("dz",path=path,xind=[0,0]))
zmax = np.squeeze(collect("ZMAX",path=path))
yO = -.5*(dy*ny)
xO = 0.0
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))[tstart:tstop+1]
a = np.squeeze(collect("alpha",path=path))
a_smooth = np.squeeze(collect("alpha_smooth",path=path))
mask = np.squeeze(collect("alpha_mask",path=path))
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
     
     return n_mmap,u_mmap,A_k,phi_mmap

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

     n,u,Ak,phi = get_data(t1,t2)
     
     print n.shape

     #nt,nx,ny = n.shape
     time = np.squeeze(collect("t_array",path=path,xind=[0,0]))[t1:t2+1]
     nx_sol = np.round(.4*nx) 
          
     #-.17949 *(dx*nx)
     data_c = phi
     print 'a',a,np.sum((np.array(np.isfinite(a))).astype(int)-1)
    

    # for i in np.arange(2.*nx/3.,nx-1):
    #      a[i,:] = a[i,0]
    # a = (np.repeat(a,nt)).reshape(nx,ny,nt)
     print a.shape
     #a = np.transpose(a,(2,0,1))

     frm_data = Frame(n,meta={'mask':True,'dx':dx,'cmap':'hot'})
     print u.shape
     #u_data =  Frame(u[:,0:10,:],meta={'mask':True,'dx':dx,'cmap':'hot'})
     u_data =  Frame(u,meta={'mask':True,'dx':dx,'cmap':'hot'})
     
     phi_data =  Frame(u,meta={'mask':True,'dx':dx,'cmap':'hot'})

     frm_exp_data = Frame(np.exp(n),meta={'mask':True,'dx':dx,'cmap':'hot'})
     
     frm_log_data = Frame(np.log(np.abs(n)),meta={'mask':True,'dx':dx,'cmap':'hot'})
     

     
     #we can include as many overplot as we want - just grab the canvas and draw whatever
     #if you are going to make movies based on stationary include nt
     alpha_contour = Frame(abs(mask),meta={'stationary':True,'dx':dx,'contour_only':True,'alpha':.3,'colors':'blue'})
     alpha_contour.nt = frm_data.nt

     a_contour = Frame(a,meta={'stationary':True,'dx':dx,'contour_only':True,'alpha':.3,'colors':'blue'})
     a_contour.nt = frm_data.nt
  
     phi_contour = Frame(phi,meta={'stationary':False,'dx':dx,'contour_only':True,'alpha':.5,'colors':'blue'})
     phi_contour.nt = frm_data.nt

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
     frm_data.reset()
     alpha_contour.reset()
     alpha_contour.nt = frm_data.nt
     alpha_contour.dx = frm_data.dx
     sigma = n.std(axis=2)
     sigma_exp = (np.exp(n)).std(axis=2)
     frm_data1D = Frame(np.average(n,axis=2),meta={'sigma':sigma,'t_array':time,'dx':dx})
     frm_exp_data1D = Frame(np.average(np.exp(n),axis=2),meta={'sigma':sigma_exp,'t_array':time,'dx':dx})
     frm_log_data1D = Frame(np.average(np.log(np.abs(n)),axis=2),meta={'sigma':(np.log(n)).std(axis=2),'t_array':time,'dx':dx})

     

     
     # sigma = n.std(axis=2)
     # frm_data1D = Frame(np.average(n,axis=2),meta={'sigma':sigma,'t_array':time,'dx':dx})


     sigma = phi.std(axis=2)
     phi_data1D = Frame(np.average(phi,axis=2),meta={'sigma':sigma,'t_array':time,'dx':dx})

     nave  = np.average(np.average(n,axis=2),axis=0)
     a_ave = np.average(a_smooth,axis=1)
     
     

     #xstart= np.int(np.round(np.mean(np.where(abs(a_ave-.1*a_ave.max()) < .05* a_ave.mean()))))
     #xstart= np.int(np.round(np.mean(np.where(abs(nave-.5*nave.max()) < .1* nave.mean()))))
     #xstart = np.int(nx*.25)
     #xstop = xstart+ nx/3.
     # xstart= np.int(np.round(np.mean(np.where(abs(nave-.6*nave.max()) < .1* nave.mean()))))
#      cond_mean = nave[np.where(nave>0)].mean()
#      #print cond_mean,nave.mean(),np.where(((1*(nave<cond_mean) +1.*(abs(nave-.1*cond_mean)<.1*cond_mean))==2 ))
#      xstop = np.int(np.mean(np.where((1*(nave>0) +1.*(abs(nave-.1*cond_mean)<.1*cond_mean))==2 )))


#      if xstop > nx:
#           xstop = nx -nx/10 ;
     
#      print nx,xstart,xstop,nave.shape,pos[0].shape,nx,dx,nx*dx,np.mgrid[xmin:nx*dx:dx,ymin:ymax:dy].shape
#      est_lam = (pos[0][xstop,5]-pos[0][xstart,5])/(np.log(nave[xstart]/nave[xstop]))
#      p0=[nave[xstart],est_lam]#
# #print 
#      popt, pcov= curve_fit(linearfall,pos[0][xstart:xstop,5],np.log(nave[xstart:xstop]),p0=p0)
#      linear_est = popt[1]
#      popt, pcov= fit_lambda(nave[xstart:xstop],pos[0][xstart:xstop,5],p0=p0)
     #p0 =[nave[xstart],est_lam,xstart/(nx-xstop)]
    
     #print nave[::xstop],nave.shape
     #print nave.shape,pos.shape
     #brutal force
     # fval_old = 1e10
     # pmin = p0
     # for x0 in [.4, .8,.9,1.0,1.1,1.2,1.5]:
     #      print 'x0', x0
     #      #p0 =[nave[xstart],est_lam,x0*xstart/(nx-xstop)]
     #      #res= fit_lambda2(nave[0:xstop],pos[0][0:xstop,5],nx-xstop,p0)
     #      xxstop = x0*xstart+ nx/2.
     #      if xxstop > nx:
     #           xxstop = nx -nx/10 ;
     #      popt, pcov= fit_lambda(nave[x0*xstart:xxstop],pos[0][x0*xstart:xxstop,5],p0=p0)
          
     #      fval = np.sum((expfall(pos[0][x0*xstart:xxstop,5],popt[0],popt[1]) - nave[x0*xstart:xxstop])**2)/(xxstop-x0*xstart)**2

     #      print 'fval: ',fval
     #      if fval<fval_old:
     #           fval_old = fval
     #           pmin = popt
     # print pmin
     
     #print pmin
     # popt, pcov= fit_lambda2(nave[0:xstop],pos[0][0:xstop,5],
     #                         p0=[nave[np.int(nx/2.0)],est_lam,np.int(nx/2.0)])
     # #print 'min parameters: ',popt,res[0]
     # n_fit = popt[0]*np.exp(-pos[0][xstart:xstop,5]/popt[1])
     # n_fit = Frame(n_fit,meta={'dx':dx,'x0':pos[0][xstart,5],'stationary':True})
     
     #frames= [frm_exp_data1D,frm_exp_data,frm_data1D,[frm_data,phi_contour]]
     #frames= [frm_data1D,frm_data,frm_exp_data1D,[frm_exp_data,phi_contour]]
     frames= [frm_data1D,[frm_data,a_contour],frm_exp_data1D,phi_data]
     #frames= [frm_data1D,[frm_data,alpha_contour],frm_exp_data1D,[phi_data,a_contour]]


     frm_data.t = 0
     frm_Ak.t = 0
     frm_Ak.reset()
     frm_data.reset()
     alpha_contour.reset()

     FrameMovie(frames,fast=True,moviename=save_path+'/'+key+str(t2),fps = 10,encoder='ffmpeg')
     #print time, n_fit.shape,popt,pcov,nave[0:40],popt
     
     frm_data.t = 0
     frm_Ak.t = 0
     t1 = t1+tchunk
     t2 = t2+tchunk


#ls -rt $PWD/movie*mp4 | awk '{ print "file " "'\''" $1 "'\''"}'

movienames = [key]#,'n_phi'+key,'u_k_phi'+key]

#from subprocess import call
for name in movienames:
     print name, save_path
     command = ('makemovlist.sh',save_path+'/',name)
     subprocess.check_call(command)



#ls -rt $PWD/movie*mp4 | awk '{ print "file " "'\''" $1 "'\''"}'
