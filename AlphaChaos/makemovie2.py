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
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from read_inp import metadata
import sys,os,inspect,shutil,subprocess
import argparse
import multiprocessing


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
nx =  np.squeeze(collect("NXPE",xind=[0,0],path=path,info=False))

mxsub = np.squeeze(collect("MXSUB",xind=[0,0],path=path,info=False)) #without gaurds

ny = nz
dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
dy = np.squeeze(collect("dz",path=path,xind=[0,0]))
zmax = np.squeeze(collect("ZMAX",path=path))
yO = -.5*(dy*ny)
xO = 0.0
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))[tstart:tstop+1]
a = np.squeeze(collect("alpha",path=path))
a_smooth = np.squeeze(collect("alpha_smooth",path=path))
beta = 5.0e-4


def get_data(start,stop):
     
     n = np.squeeze(collect("n",tind=[start,stop],path=path,info=False))
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


t1 = tstart
t2 = t1+tchunk

while t2<=tstop:

     n,u,Ak,phi = get_data(t1,t2)
     print n.shape

     nt,nx,ny = n.shape
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
     
     #we can include as many overplot as we want - just grab the canvas and draw whatever
     #if you are going to make movies based on stationary include nt
     alpha_contour = Frame(abs(np.log(a)),meta={'stationary':True,'dx':dx,'contour_only':True,'alpha':.5,'colors':'blue'})
     alpha_contour.nt = frm_data.nt
  


     frm_data_SOL = Frame(n[:,nx_sol:-1,:],meta={'mask':True,'dx':dx,'x0':dx*nx_sol})
     #frm_data = Frame(a,meta={'data_c':a,'mask':True,'dx':dx})
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
     FrameMovie([[frm_data,alpha_contour]],fast=True,moviename=save_path+'/'+'n_phi'+key+str(t2),fps = 10,encoder='ffmpeg')
     FrameMovie([frm_Ak],fast=True,moviename=save_path+'/'+'u_k_phi'+key+str(t2),fps = 10,encoder='ffmpeg')

     frm_Ak.reset()
     frm_data.reset()
     alpha_contour.reset()
     alpha_contour.nt = frm_data.nt
     alpha_contour.dx = frm_data.dx
     sigma = n.std(axis=2)
     frm_data1D = Frame(np.average(n,axis=2),meta={'sigma':sigma,'t_array':time,'dx':dx})
     

     frames= [[frm_data,alpha_contour],frm_Ak,frm_data1D,frm_data_SOL]
     #frames= [[frm_data,alpha_contour],frm_data1D]
     #frames= [[alpha_contour,frm_data],frm_data1D]
     #frames= [[alpha_contour,alpha_contour],frm_data1D]
     print time
     frm_data.t = 0
     frm_Ak.t = 0
     frm_Ak.reset()
     frm_data.reset()
     alpha_contour.reset()

     FrameMovie(frames,fast=True,moviename=save_path+'/'+key+str(t2),fps = 10,encoder='ffmpeg')
     
     
     frm_data.t = 0
     frm_Ak.t = 0
     t1 = t1+tchunk
     t2 = t2+tchunk


#ls -rt $PWD/movie*mp4 | awk '{ print "file " "'\''" $1 "'\''"}'

movienames = [key,'n_phi'+key,'u_k_phi'+key]

#from subprocess import call
for name in movienames:
     print name, save_path
     #command = ('makemovlist.sh',save_path+'/',name)
     #subprocess.check_call(command)



#ls -rt $PWD/movie*mp4 | awk '{ print "file " "'\''" $1 "'\''"}'
