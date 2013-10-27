import numpy as np
import math

import matplotlib.pyplot as plt
import hashlib
from pymongo import Connection, MongoClient
from pymongo import errors as mongoErr
from datetime import datetime
from bson import BSON

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
#from scipy.signal import argrelextrema  
from scipy.optimize import curve_fit
import sys, os, gc,types
from cStringIO import StringIO
from bson.binary import Binary
import cPickle as pickle

from boutdata import collect2 as collect
import subprocess 
from scipy import interpolate
from scipy import signal
from scipy import stats
from pympler.asizeof import asizeof,flatsize,asizesof
from read_inp import read_inp, parse_inp

copy_types = [types.StringType,types.LongType,types.FloatType,types.IntType,type( datetime.utcnow())]
ban_types = [type(np.array([123]))]
old_stdout = sys.stdout



def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    """ the sigma of gaussian is just size, and the gaussian is computed 3 sd in each dir"""
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-2*size:2*size+1, -2*sizey:2*sizey+1]
        #g = exp(-(x**2/float(size)+y**2/float(sizey)))
    g = np.exp(-(x**2/float(size)**2 + y**2/float(sizey)**2 ))
    return g / g.sum()  

def is_numeric_paranoid(obj):
    try:
        obj+obj, obj-obj, obj*obj, obj**obj, obj/obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        return True

def wave_kern(size,sizey = None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    """ the sigma of gaussian is just size, and the gaussian is computed 3 sd in each dir"""
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-3*size:3*size+1, -3*sizey:3*sizey+1]
        #g = exp(-(x**2/float(size)+y**2/float(sizey)))
    w_0 = (2.0*np.pi)/(2.0)
    g = np.pi**(-1./2.)*np.exp(-(x**2/float(size)**2 + y**2/float(sizey)**2 ))
    g = g/g.sum()
    #g = g*np.exp(w_0*complex(0,1)*x/size)*np.exp(w_0*complex(0,1)*y/size)
    g = g*(np.exp(w_0*complex(0,1)*np.sqrt(x**2+y**2)/np.sqrt(size**2 + sizey**2)) - 
           np.exp(-(1./2.) * w_0**2) )
    #g = g*(np.exp(w_0*complex(0,1)*np.sqrt(x**2+y**2)/np.sqrt(size**2 + sizey**2)))
    #g = 1.0*(np.cos(w_0*np.sqrt(x**2+y**2)/np.sqrt(size**2 + sizey**2)))
    #g = g * np.cos((np.pi*x)/(size))*np.cos((np.pi*y)/(sizey))
    #g = g * np.cos(np.sqrt((2*np.pi*x)**2 + (2*np.pi*y)**2)/np.sqrt(size**2+sizey**2))
    return g# / g.sum()
    

def blob_info(data,meta=None,label=None):

    print 'in basic_info'
    dims = data.shape
    ndims = len(dims)

    # dc = data.mean(1).mean(1).mean(1) # there MUST be a way to indicate all axis at once
    # amp = abs(data).max(1).max(1).max(1)
    
    
    nt,nx,ny = data.shape

    
    if meta is not None:
        dt = meta['dt']
        dx = meta['dx']
        dy = meta['dy']
        x0 = meta['x0']
        y0 = meta['y0']
    else:
        dt,dx,dy = (1,1,1)
        x0,y0 = (0,0)
        
    ymin =y0
    xmin =x0
    xmax =nx*dx + xmin
    ymax =ny*dy + ymin
    print ymax,xmax,ymin,xmin
    

    pos = np.mgrid[xmin:xmax-dx:nx*complex(0,1),ymin:ymax-dy:ny*complex(0,1)]
    pos_i = np.mgrid[0:nx-1:nx*complex(0,1),0:ny-1:ny*complex(0,1)]

    x = []#np.zeros(nt)
    y = []
    max_val =[]

    #moments
    print pos[0,::]
    xmoment={1:[],2:[]}
    ymoment={1:[],2:[]}
    for t in xrange(nt):
        #centered momennts, 
        moments(data[t,:,:],pos[0,:,:],appendto=xmoment)
        moments(data[t,:,:],pos[1,:,:],appendto=ymoment)
      
        max_val.append(np.max(data[t,:,:]))
        print t,':',np.max(data[t,:,:])

        
    #wavelets
        print
    freq_vs_t = wvlt(data[:,:,ny/2])
        
    c_moments =  {'x':xmoment,'y':ymoment,'t':dt*np.arange(nt),'label':label,
          'dt':dt,'max':max_val}
    
    wavelet_info = {'cwt':freq_vs_t}

    
    

    #create a blob obj
    #return Blob(cm)
    #print cm
    # print ymin
    # print ymax
    # print xmin
    # print xmax
    return c_moments

def expfall(x,y0,l):
    #print y0,l,np.min(x)
    return y0*np.exp(-x/l)

def linearfall(x,y0,l):
    #print y0,l,np.min(x)
    return y0 - x/l

def fit_lambda(y,x):
    return curve_fit(linearfall, x.flatten(), y.flatten())

def linearfit(x,y0):
    return y0 + y1*x

#popt,pcov,estimatex

#calculates the first n central moments of variable given some probabality 
#density function, for 
def moments(density,variable,nmom=2,appendto=None):
    if nmom < 1:
        return 0
    mu ={}
    mu['1'] = (np.sum(density * variable)/np.sum(density))
    
    for n in range(2,nmom+1):
        mu[str(n)] = (np.sum(density * (variable-mu['1'])**n )/np.sum(density))
        
    
    if appendto is not None:
        for n in range(1,nmom+1):
            #print len(appendto['1'])
            appendto[str(n)].append(mu[str(n)])

    return mu

def wvlt(data,J=10,dt=1,s0=None,dj=None):
    try:
        import kPyWavelet as wavelet #cwt 
        import pywt
    except:
        print 'no wavelet tools'
        return 0
    mother = wavelet.Morlet()
   
    print data.ndim

    if data.ndim == 1:
        nt = data.shape
        std1 = data.std()
        # returns a a listass with containing [wave, scales, freqs, coi, fft, fftfreqs] 
        cwt1 = wavelet.cwt(data / std1, dt, s0=s0,dj=dj,
                           wavelet=mother,J =J)

    if data.ndim == 2:
        nt,nx = data.shape
        
        cwt1 = []
        for t in xrange(nt):
            std1 = data[t,:].std()
        # returns a a list with containing [wave, scales, freqs, coi, fft, fftfreqs] 
            cwt1.append(wavelet.cwt(data[t,:] / std1, dt, wavelet=mother,s0 = nscale))
    #sig1 = wavelet.significance(1.0, dt, cwt1[1], 0, data1['alpha'], 
   # wavelet=mother)
    return cwt1



def get_data(path,field="n",fun=None,start=0,stop=50,*args):
     
     data = np.squeeze(collect("n",tind=[start,stop],path=path,info=False))
     #print data.shape,start,stop,data.dtype.name
     #data_mmap = np.memmap(data,dtype=data.dtype.name,mode='w+',shape=data.shape)
     #data_mmap[:] = data[:]
     
     #print 'n_mmap.shape :',n_mmap.shape,n.shape
     #del data
     gc.collect()
     
     if fun is None:
         return data
     else:
         return fun(data)
 
class field_info(object):
    def __init__(self,path,field="n",meta=None,fast_center=True,get_Xc=True,
                 get_lambda=True,debug=False,fast = True):
    
        self.path=path
        self.field=field
        #print 'read the inp'
        f=open(path+'/BOUT.log.0', 'r')
        self.md5 = hashlib.md5(f.read()).hexdigest()
        f.close()


        
        defaults = {'dx':1,'x0':0,'dy':1,'y0':0,'log_n':False}   

        for key,val in defaults.items():
            if not hasattr(self,key):
                print 'setting: ',key,val
                setattr(self, key, val)
                
        #read all metadata from BOUT.inp
        #print 'read the inp'
        inp = read_inp(path=path,boutinp='BOUT.inp')
        inp = parse_inp(inp) #inp file
        
        
        for k, v in inp.items():
                setattr(self, k, v)
          
        for k,v in (getattr(self,'[main]')).items():
            print k,v
            try:
                setattr(self, str.lower(k), float(v))
            except:
                setattr(self, str.lower(k), v)


        for k,v in (getattr(self,'[physics]')).items():
            print k,v
            try:
                setattr(self, str.lower(k), float(v))
            except:
                try:
                    setattr(self, str.bool(k), v)
                except:
                    setattr(self, str.lower(k), v)
       

        if meta is not None:
            for k, v in meta.items():
                setattr(self, k, v)

      
        self.mz = self.mz-1
        self.nz  = self.mz       #self.nz = self.np.squeeze(collect("MZ",xind=[0,0],path=path,info=False))-1
        nxpe = self.nxpe# np.squeeze(collect("NXPE",xind=[0,0],path=path,info=False))
        mxg = self.mxg #np.squeeze(collect("MXG",xind=[0,0],path=path,info=False))
        
        self.mxsub = np.squeeze(collect("MXSUB",xind=[0,0],path=path,info=False)) #gaurds
        self.nx = int(self.mxsub*nxpe + 2*mxg) #including gaurds
        nx,ny,nz = self.nx,self.nz,self.nz
        
        self.dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
        self.dy = np.squeeze(collect("dz",path=path,xind=[0,0]))
        self.zmax = np.squeeze(collect("ZMAX",path=path))

        self.time = np.array(np.squeeze(collect("t_array",path=path,xind=[0,0])))
        nt = self.nt = len(self.time)

        try:
            self.nave = np.squeeze(collect("nave",path=path))
            self.nave = self.nave/nt
            have_nave = True
        except:
            have_nave = False

        try:
            self.G = np.squeeze(collect("gamma",path=path))
            self.G = self.G/nt
            have_flux = True
        except:
            have_flux = False
           
        self.debug = debug

        if debug:        
            t_chunk = 20
            t_stop  = np.max(self.nt)
            t1 = np.int(self.nt*.90)
        else:
            t_chunk = 30
            t_stop  = np.max(self.nt)
            t1 = 0

        self.t_chunk  = t_chunk  
        self.t1 = t1
        self.t_stop = t_stop
        tarray = np.array(np.arange(nt - nt/2))
        
        t2 = np.min([t1+t_chunk,t_stop])

        self.xmoment = {'1':[],'2':[]}
        self.ymoment = {'1':[],'2':[]}
        self.max_val = []
        
        
        y0,x0 = (self.y0,self.x0)
        dx,dy = (self.dx,self.dy)

        
        ymin = y0
        xmin = x0
        xmax = nx*dx + xmin
        ymax = ny*dy + ymin

        self.kx_max = nx
        self.ky_max = ny
        kxmax,kymax = nx,ny
        kxmin = 0.
        kymin = 0.
   
        dky = self.dky = 1.0/self.zmax
        dkx = self.dkx = 2.*np.pi/self.dx
        self.kx = np.arange(0,nx)*self.dkx
        self.ky = np.arange(0,ny)*self.dky

        #print nx,nx*self.mxsub
        self.pos_i = np.mgrid[0:nx:1,0:ny:1]
        self.x = np.mgrid[xmin:xmax:self.dx]
        self.y = np.mgrid[xmin:xmax:self.dy]
        
       
        sys.stdout = mystdout = StringIO()
        a = np.squeeze(collect("alpha",path=path))
        sys.stdout = old_stdout
        

        #self.a_ave 
        a_ave = self.a_ave = np.average(a,axis=1)
        # popt, pcov= curve_fit(fit_tanh,xnew,a_ave,p0=p0)
        # w = run['w'] =  popt[0]*dx

        
        try:  
            xstart = np.min(np.where((a_ave-.1*a_ave.max()>0)) )
        
            self.xstart = xstart
            
            xstop = np.int(xstart+ nx/2.)
            if xstop > nx:
                xstop = nx -1;
        except:
            xstart= np.int(nx/3.)            
            xstop = np.int(xstart+ nx/2.)

        self.fast = fast
        #will return 'real n' not 'log n' statistics and profiles
        moving_min_n, moving_max_n, moving_min, moving_max, nave_net, nrms_net ,nmin_net,nmax_net = self.compute_profile() 
        pdf_x = []
        df_x = []
        
        
        self.nbins = nbins = 300
        for x in xrange(nx):
            #print moving_min_n[x]
            pdf_x.append(np.mgrid[moving_min_n[x]-.01:moving_max_n[x] +.01:complex(0,nbins)])
            df_x.append(np.mgrid[moving_min[x]-.01:moving_max[x] +.01:complex(0,nbins)])
            
        
        
        self.pdf_x = pdf_x
        self.df_x = df_x
        pdf_y=nx*[None]
        df_y=nx*[None]
        pdf_y_net = 0
        df_y_net = 0
        # self.moving_min_n= moving_min_n
        # self.moving_max_n= moving_max_n
        # self.moving_min= moving_min
        # self.moving_max= moving_max


        for i in range(nx):
            #print "some unique object %d" % ( i, )
            pdf_y[i]= 0*pdf_x[i]
            df_y[i]= 0*pdf_x[i]

        
        
        xstart = np.min(