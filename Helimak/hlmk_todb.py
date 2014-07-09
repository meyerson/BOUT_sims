#!/opt/apps/python/epd/7.2.2/bin/python
import sys, os, gc
from datetime import datetime
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


from scipy.interpolate import UnivariateSpline

import matplotlib
matplotlib.use('Agg')

from boutdata import collect_alt as collect

from read_inp import read_inp, parse_inp
import sys,os,inspect,shutil,subprocess
import argparse
import multiprocessing
from scipy.optimize import curve_fit ,fmin#,minimize
from scipy import signal
from frame import Frame, FrameMovie
cmd_folder = HOME+'/BOUT_sims/blob_py'

from pymongo import Connection, MongoClient
from pymongo import errors as mongoErr
from bson.binary import Binary
import cPickle as pickle
import zlib

def zdumps(obj):
  return zlib.compress(pickle.dumps(obj,pickle.HIGHEST_PROTOCOL),9)
def zloads(zstr):
  return pickle.loads(zlib.decompress(zstr)) 
import hashlib,types
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

parser.add_argument("-s","--server", type=str,
                    help="server url",default='dataplay.no-ip.org')

parser.add_argument('--debug',dest='debug',action='store_true')
parser.set_defaults(debug=False)

args = parser.parse_args()

path = args.path
key = args.key
tstart = args.tstart
tstop = args.tstop
tchunk = args.tchunk
debug = args.debug
server = args.server


save_path = path.replace('scratch','work')+'movie'

import matplotlib
if debug:
     matplotlib.use('pdf')
     from matplotlib.backends.backend_pdf import PdfPages
else:
     matplotlib.use('Agg')

from boutdata import collect_alt as collect

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
# from boutdata import collect
# from boutdata import collect2 
# from collect2 import collect2 as collect


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
dy = 2*np.pi*zmax/nz
yO = -.5*(dy*ny)
xO = 0.0
time = np.squeeze(collect("t_array",path=path,xind=[0,0]))#[tstart:tstop+1]
dt = time[1]-time[0]
try:
     bias_phi = np.squeeze(collect("bias_phi",path=path))#[tstart:tstop+1]  
except:
     bias_phi = 0
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


def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]


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
     ydata = args[0]
     x = args[1]
     mx = args[2]


     y0 = params[0]
     l = params[1]
     start = params[2]

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

c = MongoClient(host=server)
db = c.hlmk
db.authenticate("lonestar","britoebito")
comments = db.comments
hlmk_blob = {"author": "Dmitry",
            "text": "My 2nd blog post!",
            "tags": ["fusion", "in", "5 days"],
            "date": datetime.utcnow()}
comments.insert(hlmk_blob)
n_coll = db.n_profiles
phi_coll = db.phi_profiles
Te_coll = db.Te_profiles
V_coll = db.V_profiles
meta_coll = db.meta


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

ban_types = [type(np.array([123]))]

def is_numeric_paranoid(obj):
    try:
        obj+obj, obj-obj, obj*obj, obj**obj, obj/obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        return True

def serialize(input_dict):
     ser_dict = {}
     
          # elif value.size == 1:
          #      sim_blob[key] = value

     for key, value in input_dict.iteritems():
          if type(value) == type(np.array([12])):
               print key, value.nbytes,value.size
               if value.ndim == 1:
                    value = value.tolist()


          if is_numeric_paranoid(value):
               if 'all' in dir(value):
                    if value.ndim == 1 and value.size == 1 or value.ndim == 0:
                         ser_dict[key] = np.asscalar(value)
                    else:
                         print key, ' to opaque binary zlib obj',type(value)
                         ser_dict[key] = Binary(zdumps(value))
               else:
                    ser_dict[key] = value
          elif type(value) is types.ListType:
                if len(value) == 0:
                    ser_dict[key] = value
                elif type(value[0]) == type(np.array([12.23])):
                    print key, ' to opaque binary zlib obj',type(value),type(value[0])
                    ser_dict[key] = Binary(zdumps(value))
                else:
                    if type(value[0]) == type(np.int64(1)):
                        value = [ int(x) for x in value ]
                    ser_dict[key] = value
          elif type(value) is types.DictType:
               ser_dict[key] = value
          else:
               #print key, 'serializing', value, type(value), value.size()
               #ser_dict[key] = Binary(pickle.dumps(value,pickle.HIGHEST_PROTOCOL))
               print key, ' to opaque binary zlib obj',type(value)
               ser_dict[key] = Binary(zdumps(value))

     return ser_dict

# def to_json(dict_obj):
#      for key, value in dict_obj.iteritems():
#           try:
#                print key, value.__class__, type(value),type(np.array([12]))
#                dict_obj[key] = np.asscalar(value)
#           except:
#                print key
          

#      return dict_obj

# log_f=open(path+'/BOUT.log.0', 'r')
# f.close()
while t2<=tstop:
     #basic for join in our little database

     idstring = path+str(t1)+str(t2)
     idstring = hashlib.md5(idstring).hexdigest()

     n,u,Ak,phi,Te = get_data(t1,t2)
     
     print n.shape

     nt,nx,ny = n.shape

     
     time = np.squeeze(collect("t_array",path=path,xind=[0,0]))[t1:t2+1]
     phi_bias = np.squeeze(collect("bias_phi",zind=[0,0],path=path))

     
     #eventually one would want to dump the entire contents of BOUT.inp here
     meta_dict = {"dx":dx,"dy":dy,"t1":t1,"t2":t2,"dt":dt,"nx":nx,
                  "ny":ny,"nt":nt}

     [n_dict,phi_dict,Te_dict,u_dict] = map(lambda nd_obj: {"xt":np.ravel(np.mean(nd_obj,axis=2)),"max":np.ravel(np.max(nd_obj,axis=2)),"min":np.ravel(np.min(nd_obj,axis=2)),"sig":np.ravel(np.std(nd_obj,axis=2))},[np.exp(n),phi,Te,u])

     profile_dicts = [n_dict,phi_dict,Te_dict,u_dict]
     profile_coll = [n_coll,phi_coll,Te_coll,u_coll]
     

     # n_dict = {"n_xt":np.ravel(np.mean(np.exp(n),axis=2)),"n_max":np.ravel(np.max(np.exp(n),axis=2)),"n_min":np.ravel(np.min(np.exp(n),axis=2)),"n_sig":np.ravel(np.std(np.exp(n),axis=2))}

     print n_dict["n_xt"].shape
     

     try:
          z = serialize(meta_dict).copy()
          z.update({"_id":idstring})
          meta_coll.insert(z)
     except:
          meta_coll.update({"_id":idstring},serialize(meta_dict))
      
     for coll_elem,dict_elem in zip(profile_dicts):
          try:
               z = serialize(dict_elem).copy()
               z.update({"_id":idstring})
               coll_elem.insert(dict_elem)
          except:
               coll_elem.update({"_id":idstring},serialize(dict_elem))
     
     
     nx_sol = np.round(.4*nx) 
          

     data_c = phi
     print 'a',a,np.sum((np.array(np.isfinite(a))).astype(int)-1)
     #ddtu = np.squeeze(collect("t_array",path=path,xind=[0,0]))[t1:t2+1]

    # for i in np.arange(2.*nx/3.,nx-1):
    #      a[i,:] = a[i,0]
    # a = (np.repeat(a,nt)).reshape(nx,ny,nt)
     print a.shape
     #a = np.transpose(a,(2,0,1))

     frm_n = Frame(np.exp(n),meta={'dx':dx,'dy':dy,'title':r'$n_{AC}$','cmap':'hot',
                           'xlabel':'x['+r'$\rho_s$'+']',
                           'ylabel':'y['+r'$\rho_s$'+']',
                           'interpolation':'linear','grid':False,
                           'linewidth':1,'contour_color':'black',
                           't_array':time,'x0':dx*250.0 })

     n_v = np.gradient(np.exp(n))
     
     cond = n_v[1]<1e-12
     n_v[1] = n_v[1]+(cond)*np.min(n_v[1][np.where(n_v[1]>1e-12)])
     cond = n_v[2]<1e-12
     n_v[2] = n_v[2]+(cond)*np.min(n_v[2][np.where(n_v[2]>1e-12)])

     n_v = [n_v[0]/n_v[1],n_v[0]/n_v[2],n_v[0]/np.sqrt(n_v[1]**2+ n_v[2]**2)]

     n_v_display = map(lambda v: np.sign(v)*np.log(np.exp(-(v/np.mean(np.abs(v)))**2)+np.abs(v)),n_v)

     frm_dn = Frame(n_v_display[1],meta={'dx':dx,'dy':dy,'title':r'$n_{AC}$','cmap':'hot',
                                 'xlabel':'x['+r'$\rho_s$'+']',
                                 'ylabel':'y['+r'$\rho_s$'+']',
                                 'interpolation':'linear','grid':False,
                                 'linewidth':1,'contour_color':'black',
                                 't_array':time,'x0':dx*250.0,'stationary':False })

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
                         meta={'dx':dx,'dy':dy,'title':r'$n_{AC}$','cmap':'hot',
                               'xlabel':'x['+r'$\rho_s$'+']',
                               'ylabel':'y['+r'$\rho_s$'+']',
                               'fontsz':10,'interpolation':'linear','grid':False,
                               'linewidth':1,'contour_color':'black',
                               't_array':time,'x0':0})


     vyEB = -((np.gradient(phi)[1])/dx)/B0

     dky = 1.0/zmax
     allk = dky*np.arange(ny/8.0)+(1e-8*dky)
  
     dens = np.exp(n)
     dens_fft = np.fft.rfft(dens)
     dens_pow = dens_fft.conj()*dens_fft
    
     k_max = [((np.where(col == np.max(col)))[0])+1 for col in np.mean(dens_pow,axis=0)[:,1:]]
    
     
     dens_acorr = np.real(np.fft.irfft(dens_pow))
 
     omega = np.gradient(np.angle(dens_fft))[0]

  
     nt,nx,nky = dens_fft.shape

     omega_r = np.real(omega[:,:,0:ny/8.0])
  
     omega_r = (omega_r  - np.sign(omega_r)*np.pi)*(np.abs(omega_r)>np.pi/2)+omega_r*(np.abs(omega_r)<np.pi/2)
    
     vy_phase = -1.0*np.array(map(lambda omega_r_x: omega_r_x/(dt*allk),omega_r))

     weight = np.real(dens_pow[:,:,1])
     ave_vy_phase = np.sum(vy_phase[:,:,1] * weight,axis=0)/np.sum(weight,axis=0)
     
     #we can also wash out the noise in the time resolved version
     #with spectral power weighted splines
     vy_phase_s = []
     for jw,vy in enumerate(vy_phase[:,:,1]):
          s = UnivariateSpline(dx*np.arange(nx), vy, w = weight[jw,:], s=.1,k=1)
          vy_phase_s.append(s(dx*np.arange(nx)))

     vy_phase_s = np.array(vy_phase_s)

     vy_phase_frm =  Frame(vy_phase_s,meta={'t_array':time,'dx':dx,'yscale2':5e3,'ylabel2':r'$\frac{m}{s}$'})
     vy_phase_stat =  Frame(ave_vy_phase,meta={'stationary':True,
                                                 'dx':dx,'x0':0,
                                                 'yscale2':5e3,'ylabel2':r'$\frac{m}{s}$'})


     pow_frm =  Frame(np.sqrt(np.real(dens_pow[:,:,1]/dens_pow[:,:,2])),
                      meta={'stationary':False,
                            'dx':dx,'x0':0,'yscale':'log',
                            'yscale2':1,'ylabel2':''})

    
     mean_vyExB = np.mean(np.exp(n)*vyEB,axis=2)/np.mean(np.exp(n),axis=2)
     import copy
     sigma = vyEB.std(axis=2)
 

     vy_1D =  Frame(mean_vyExB,meta={'t_array':time,'dx':dx,'yscale2':5e3,
                                     'ylabel2':r'$\frac{m}{s}$',
                                     'sigma':sigma})
     vy_1D_static = Frame(np.mean(mean_vyExB,axis=0),meta={'dx':dx,
                                                                   'stationary':True,'xlabel':'x['+r'$\rho_s$'+']',
                                                                   'ylabel':r'$v_y[C_s]$','yscale2':5e3,'ylabel2':r'$\frac{m}{s}$'})
                              


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
     
  

     phi_contour = Frame(phi,meta={'stationary':False,'dx':dx,'contour_only':True,'alpha':.5,'colors':'red'})
     phi_contour.nt = frm_n.nt

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
  
     
     a_m = np.power(beta/a_smooth[:,0]**2,.20)
     a_par = np.power(beta/a_smooth[:,0]**2,1.0/3.0)
     
     a_D = beta/(a_smooth[:,0] * mu)

     a_mu = np.power(mu/a_smooth[:,0],.25)

     a_L = np.power((beta/a_smooth[:,0]**2),1./3.)

     


     print a_m
     
     frm_Ak = Frame(Ak[:,:,0:60],meta={'dy':dky,'dx':dx,
                                        'overplot':[2.*np.pi/a_m,2.*np.pi/a_L,
                                                    2.*np.pi/a_mu,2.*np.pi/a_D]})
  
     frm_Ak.reset()
     frm_n.reset()

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

     sigma = phi.std(axis=2)
     phi_data1D = Frame(np.average(phi/4.7,axis=2),meta={'sigma':sigma/4.7,'t_array':time,
                                                         'dx':dx,'overplot':phi_bias/4.7})
     sigma = u.std(axis=2)
     u_data1D = Frame(np.average(u/4.7,axis=2),meta={'sigma':sigma,'t_array':time,'dx':dx})
     

     nave  = np.average(np.average(n,axis=2),axis=0)
     a_ave = np.average(a_smooth,axis=1)
   
     print 'new_frm?'
  
     #frames= [vy_phase_frm ,frm_dn,[frm_Te1D,phi_data1D],[vy_1D_static,vy_1D]]
     frames= [vy_phase_stat,frm_n,[frm_Te1D,phi_data1D],[vy_1D_static,vy_phase_stat,vy_1D]]

      
     frm_n.t = 0
    
     
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
for name in movienames:
     print name, save_path
     command = ('makemovlist.sh',save_path+'/',name)
     subprocess.check_call(command)



#ls -rt $PWD/movie*mp4 | awk '{ print "file " "'\''" $1 "'\''"}'
