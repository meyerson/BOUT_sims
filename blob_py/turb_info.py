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
def mean_and_std(data):
    dshape = data.shape
    
    DC = np.swapaxes(data,1,0)
    DC = DC.reshape([dshape[1], np.prod(dshape) / dshape[1]])
    
    std = DC.std(axis=1)
    DC  = DC.mean(axis=1)


    DCC = np.repeat(DC,np.prod(dshape) / dshape[1])
    DCC = DCC.reshape([dshape[1], dshape[0], dshape[2]])
    print 'DCC: ', DCC.shape
    DCC = np.swapaxes(DCC,1,0)

    sigma2 =  np.repeat(std,np.prod(dshape) / dshape[1])
    sigma2 = sigma2.reshape([dshape[1], dshape[0], dshape[2]])
    sigma2 = np.swapaxes(sigma2,1,0)
    
    AC = data - DCC
    sigma2 = sigma2 + np.mean(sigma2)*(sigma2<(np.mean(sigma2)*1e-3))
    AC_norm = AC/sigma2

    return DC, std, AC, AC_norm

     
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
    #print pos[0,::]
    xmoment={1:[],2:[]}
    ymoment={1:[],2:[]}
    for t in xrange(nt):
        #centered momennts, 
        moments(data[t,:,:],pos[0,:,:],appendto=xmoment)
        moments(data[t,:,:],pos[1,:,:],appendto=ymoment)
      
        max_val.append(np.max(data[t,:,:]))
        #print t,':',np.max(data[t,:,:])

        
    #wavelets
        #print
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
                #print 'setting: ',key,val
                setattr(self, key, val)
                
        #read all metadata from BOUT.inp
        #print 'read the inp'
        inp = read_inp(path=path,boutinp='BOUT.inp')
        inp = parse_inp(inp) #inp file
        
        
        for k, v in inp.items():
                setattr(self, k, v)
          
        for k,v in (getattr(self,'[main]')).items():
            #print k,v
            try:
                setattr(self, str.lower(k), float(v))
            except:
                setattr(self, str.lower(k), v)


        for k,v in (getattr(self,'[physics]')).items():
            #print k,v
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
            havex_flux = False
           
        self.debug = debug

        if debug:        
            t_chunk = 20
            t_stop  = np.max(self.nt)
            t1 = np.int(self.nt*.90)
        else:
            t_chunk = 30
            t_stop  = np.max(self.nt)
            t1 = 0 #self.time[0]

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
        #will return 'real n' not 'log n' statistics and profiles, spasly sampled in x, dense in t
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

        flux_pdf_y=nx*[None]
        flux_df_y =nx*[None]
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
            flux_pdf_y[i]= 0*pdf_x[i]
            flux_df_y[i]= 0*pdf_x[i]

        
        
        xstart = np.min(np.where(a_ave - 0.1 * a_ave.max() > 0))
        xstop = np.int(xstart + nx / 2.0)
        x1 = np.int(xstart)
        x2 = np.int(xstart + nx/5.0)
        # x1 = np.min(np.where((a_ave-.1*a_ave.max()>0)) )
        # x2 = np.min(np.where((a_ave-.9*a_ave.max()>0)) )

        df_x_start = np.array(df_x[x1:x2]).min()
        df_x_stop = np.array(df_x[x1:x2]).max()
        pdf_x_start = np.array(pdf_x[x1:x2]).min()
        pdf_x_stop = np.array(pdf_x[x1:x2]).max()
        self.df_x_start = df_x_start
        self.df_x_stop = df_x_stop
        self.pdf_x_start = pdf_x_start
        self.pdf_x_stop = pdf_x_stop
        if xstop > nx:
            xstop = nx - 1
            
        j = 0
        print 't2,t_stop:', t2, t_stop
        self.nave = []
        self.flux = []
        self.v = []
        while t2 < t_stop:
            print 'computing statistical measures'
            (nave, dn), (flux, flux_max, d_flux), (vx,vy),data = self.dens_and_flux(t1, t2)
            self.save_linlam(np.log(nave), np.log(dn))
            self.nave.append([nave, dn])
            self.flux.append([flux, flux_max, d_flux])
            self.v.append([vx,vy])
            #print len(self.flux)
            ntt, nxx, nyy = data.shape
            #print 'data.shape: ', data.shape, t1, t_stop
            if t1 > np.int(2*t_stop)/3:
                j = j + 1
                for x in xrange(nx):
                    if not fast:
                        p = np.polyfit(np.arange(ntt), np.squeeze(data[:, x, :]).mean(axis=1), 1)
                        nave2d = p[1] + p[0] * np.arange(ntt)
                        
                        p = np.polyfit(np.arange(ntt), np.squeeze(flux[:, x, :]).mean(axis=1), 1)
                        flux2d_x = p[1] + p[0] * np.arange(ntt)
                        if self.log_n:
                            cond = np.exp(data[:, x, :]) - np.exp(np.reshape(np.repeat(nave2d, ny), [ntt, ny]))
                            cond_flux =  np.exp(data[:, x, :])*vy - np.reshape(np.repeat(flux2d_x, ny), [ntt, ny])
                        else:
                            cond = data[:, x, :] - np.reshape(np.repeat(nave2d, ny), [ntt, ny])
                    elif self.log_n:
                        #print nave.shape,data.shape
                        cond = np.exp(data[:, x, :]) - nave[x]
                    else:
                        cond = data[:, x, :] - nave[x]
                    nrms = cond.std()
                    datamax = data[:, x, :].max()
                    if np.isfinite(datamax):
                        dfadd = np.squeeze(np.histogram(cond, bins=nbins, normed=False, range=(df_x[x].min(), df_x[x].max()))[0])
                        pdfadd = np.squeeze(np.histogram(cond / nrms, bins=nbins, normed=False, range=(pdf_x[x].min(), pdf_x[x].max()))[0])
                        if np.isfinite(np.sum(pdfadd)):
                            pdf_y[x] = pdf_y[x] + pdfadd
                            df_y[x] = df_y[x] + dfadd
                            if abs(x - xstart) < nx / 20.0:
                                pdf_y_net = pdf_y_net + np.squeeze(np.histogram(cond / nrms, bins=nbins, normed=False, range=(pdf_x_start, pdf_x_stop))[0])
                                #print x1, x2
                                df_y_net = df_y_net + np.squeeze(np.histogram(cond, bins=nbins, normed=False, range=(df_x_start, df_x_stop))[0])

            self.pdf_y = pdf_y
            self.pdf_y_net = pdf_y_net
            self.df_y = df_y
            self.df_y_net = df_y_net
            self.t.append(t1 + 0.5 * (t2 - t1))
            t1 = t2 + 1
            t2 = np.min([t1 + t_chunk - 1, t_stop + 1])

        self.nave = np.array(self.nave)
        self.flux = np.array(self.flux)
        self.pdf_stats = []
        self.df_stats = []
        print 'almost done'
        for x in xrange(nx):
            if np.mean(nrms_net[:, x]) != 0:
                self.pdf_y[x] = pdf_y[x] / sum(pdf_y[x])
                self.df_y[x] = df_y[x] / sum(df_y[x])
            else:
                self.pdf_x[x] = pdf_x[x - 1]
                
        if hasattr(self, 'pdf_y_net'):
            self.pdf_y_net = pdf_y_net / np.float(sum(pdf_y_net))
            self.df_y_net = df_y_net / np.float(sum(df_y_net))

    def save_v(self):
        return 0
    
    def save_pow(self):
        return 0
    

    def save_linlam(self, nave, dn):
        a_ave = self.a_ave
        nx = self.nx
        xstart = self.xstart
        xstop = np.int(xstart + nx / 2.0)
        if xstop > nx:
            xstop = nx - 1
        if self.log_n:
            est_lam = (self.x[xstop] - self.x[xstart]) / (nave[xstart] - nave[xstop])
            p0 = [nave[xstart], est_lam]
            popt, pcov = curve_fit(linearfall, self.x[xstart:xstop], nave[xstart:xstop], p0=p0)
        else:
            est_lam = (self.x[xstop] - self.x[xstart]) / np.log(nave[xstart] / nave[xstop])
            p0 = [np.log(nave[xstart]), est_lam]
            popt, pcov = curve_fit(linearfall, self.x[xstart:xstop], np.log(nave[xstart:xstop]), p0=p0)
        self.linlam.append(popt[1])

    def dens_and_flux(self, t1, t2):
        sys.stdout = mystdout = StringIO()
        path = self.path
        dy = self.dy
        dx = self.dx
        # ny = self.ny
        # nx = self.nx
        
        if self.fast:
            data = np.squeeze(collect('n', tind=[t1, t1 + 1], path=path, info=False))
            ntt, nxx, nyy = data.shape
            vx = np.squeeze(np.gradient(np.squeeze(collect('phi', tind=[t1, t1 + 1], path=path, info=False)))[2]) / dy
            vy = np.squeeze(np.gradient(np.squeeze(collect('phi', tind=[t1, t1 + 1], path=path, info=False)))[1]) / dx
        else:
            try:
                data = np.squeeze(collect('n', tind=[t1, t2], path=path, info=False))
                vx = np.squeeze(np.gradient(np.squeeze(collect('phi', tind=[t1, t2], path=path, info=False)))[2]) / dy
                vy = np.squeeze(np.gradient(np.squeeze(collect('phi', tind=[t1, t2], path=path, info=False)))[1]) / dx
            except:
                data = np.squeeze(collect('n', tind=[t1, t2 - 1], path=path, info=False))
                vx = np.squeeze(np.gradient(np.squeeze(collect('phi', tind=[t1, t2 - 1], path=path, info=False)))[2]) / dy
                vy = np.squeeze(np.gradient(np.squeeze(collect('phi', tind=[t1, t2 - 1], path=path, info=False)))[1]) / dx

        if self.log_n:
            flux = np.exp(data) * vx
        else:
            flux = data * vx
        sys.stdout = old_stdout
        if self.log_n:
            n_ave = np.exp(data).mean(axis=0)
            n_ave = n_ave.mean(axis=1)
        else:
            n_ave = n.mean(axis=0)
            n_ave = n_ave.mean(axis=1)


        # b_DC, b_std,b_AC,b_AC_norm = mean_and_std(flux)
        # blob_fft = np.fft.rfft(b_AC_norm)
        # blob_power  = blob_fft.conj()*blob_fft
        # blob_R = np.real(np.fft.ifft(blob_power,axis=2)) #n 
        # temp = np.mean(blob_R[:,:,0:np.int(ny/8)],axis=0)
        # norm = np.amax(temp,axis=1)
        # norm = (np.repeat(norm,np.int(ny/3))).reshape(nx,np.int(ny/3))
        # temp = temp/norm
        # norm = (np.repeat(norm,np.int(ny/8))).reshape(nx,np.int(ny/8))
        # temp = temp/norm

        # z = np.array([((np.where(col > np.exp(-.5)*np.max(col)))[0]) for col in temp[:,:]])
        # for i,elem in enumerate(z):
        #     if len(elem) >0:
        #         z[i] = 2.0*(1./np.sqrt(2))*dy*(len(elem) + (-np.exp(-.5)+temp[i,max(elem)])/(temp[i,max(elem)]- temp[i,max(elem)+1]))


        

        vx = np.mean(np.mean(vx,axis=0),axis=1)
        vy = np.mean(np.mean(vy,axis=0),axis=1)
        
        dshape = data.shape
        if self.log_n:
            temp = np.transpose(data, [1, 0, 2])
        else:
            temp = np.transpose(data, [1, 0, 2])
        temp = np.reshape(temp, [dshape[1], np.prod(dshape) / dshape[1]])
        n_sigma = temp.std(axis=1)

        temp = np.transpose(flux, [1, 0, 2])
        temp = np.reshape(temp, [dshape[1], np.prod(dshape) / dshape[1]])

        flux_sigma = temp.std(axis=1)
        flux_max = temp.max(axis=1)
        flux_ave = flux.mean(axis=0)
        flux_ave = flux.mean(axis=1) #nx long, ave along y and t

        flux_hist = flux-flux_ave
        n_hist = n - n_ave

        return ((n, n_sigma), (flux, flux_sigma, flux_max),(vx,vy), data)

    def compute_profile(self):
        self.linlam = []
        self.lam = []
        self.t = []
        self.pdf_x = []
        self.pdf_y = []
        self.df_x = []
        self.df_y = []
        nrms = []
        nave = []
        nmax = []
        nmin = []
        coarse_x = []
        xchunk = 10
        nx = self.nx
        nt = self.nt
        path = self.path
        xpos = np.arange(xchunk) * nx / xchunk
        xpos = list(xpos)
        xpos.append(nx - 1)
        self.nave = []
        for x in np.array(xpos):
            #print 'x: ', x / nx, nt
            coarse_x.append(x)
            #print x, round(nt / 2)
            sys.stdout = mystdout = StringIO()
            data = np.squeeze(collect('n', xind=[int(x), int(x) + 5], tind=[nt / 2, nt - 2], path=path, info=False))
            if self.log_n:
                data = np.exp(data)
            sys.stdout = old_stdout
            dshape = data.shape
            temp = np.reshape(data, [dshape[0], np.prod(dshape) / dshape[0]])
            nrms.append(temp.std(axis=1))
            try:
                nave.append(data.mean(axis=1).mean(axis=1))
                nmax.append(data.max(axis=1).max(axis=1))
                nmin.append(data.min(axis=1).min(axis=1))
            except:
                nave.append(data.mean(axis=1))
                nmax.append(data.max(axis=1))
                nmin.append(data.min(axis=1))

        nave = np.array(nave)
        nrms = np.array(nrms)
        nmax = np.array(nmax)
        nmin = np.array(nmin)
        print 'nave.shape: ', nave.shape, nmax.shape,
        print len(coarse_x), len(np.arange(nt - nt / 2))
        tarray = np.array(np.arange(nt - nt / 2 - 1))
        self.nave_net = np.transpose(nave)
        self.nrms = np.transpose(nrms)
        self.nmin = np.transpose(nmin)
        self.nmax = np.transpose(nmax)
        self.coarse_x = coarse_x
        self.moving_min = np.min(self.nmin - self.nave_net, axis=0)
        self.moving_max = np.max(self.nmax - self.nave_net, axis=0)
        self.moving_min_n = np.min((self.nmin - self.nave_net) / self.nrms, axis=0)
        self.moving_max_n = np.max((self.nmax - self.nave_net) / self.nrms, axis=0)
        nmin_net = []
        nmax_net = []
        nrms_net = []
        nave_net = []
        coarse_x = self.coarse_x
        xnew = np.arange(0, nx)
        for tt in tarray:
            f = interpolate.interp1d(coarse_x, nave[:, tt], kind='linear')
            nave_net.append(f(xnew))
            f = interpolate.interp1d(coarse_x, nrms[:, tt], kind='linear')
            nrms_net.append(f(xnew))
            f = interpolate.interp1d(coarse_x, nmin[:, tt], kind='linear')
            nmin_net.append(f(xnew))
            f = interpolate.interp1d(coarse_x, nmax[:, tt], kind='linear')
            nmax_net.append(f(xnew))

        nave_net = np.array(nave_net)
        nrms_net = np.array(nrms_net)
        nmin_net = np.array(nmin_net)
        nmax_net = np.array(nmax_net)
        moving_min_n = np.min((nmin_net - nave_net) / nrms_net, axis=0)
        moving_max_n = np.max((nmax_net - nave_net) / nrms_net, axis=0)
        moving_min = np.min(nmin_net - nave_net, axis=0)
        moving_max = np.max(nmax_net - nave_net, axis=0)
        return (moving_min_n,
         moving_max_n,
         moving_min,
         moving_max,
         nave_net,
         nrms_net,
         nmin_net,
         nmax_net)

    def serialize(self, input_dict):
        ser_dict = {}
        for key, value in input_dict.iteritems():
            #print '#######################'
            if is_numeric_paranoid(value):
                if 'all' in dir(value):
                    if value.ndim == 1 and value.size == 1 or value.ndim == 0:
                        ser_dict[key] = np.asscalar(value)
                    else:
                        ser_dict[key] = Binary(pickle.dumps(value, protocol=2))
                else:
                    ser_dict[key] = value
            elif type(value) in copy_types:
                ser_dict[key] = value
            elif type(value) is types.ListType:
                #print 'list: ', key, len(value)
                if len(value) == 0:
                    ser_dict[key] = value
                elif type(value[0]) == type(np.array([12.23])):
                    #print value[0].size, len(value), type(value[0])
                    #print key, type(value[0]), type(value[0]) == type(np.array([12.23]))
                    ser_dict[key] = Binary(pickle.dumps(value, protocol=2))
                else:
                    #print len(value), type(value[0])
                    #print type(np.array([12.23]))
                    if type(value[0]) == type(np.int64(1)):
                        value = [ int(x) for x in value ]
                    ser_dict[key] = value
            elif type(value) is types.DictType:
                ser_dict[key] = value
            else:
                #print key, 'serializing', value, type(value), value.size()
                ser_dict[key] = Binary(pickle.dumps(value, protocol=2))

        return ser_dict

    def to_db(self, server = 'beavis.ph.utexas.edu'):
        print 'in to_db'
        sim_blob = {'author': 'Dmitry',
         'tags': ['fusion', 'in', '50 years'],
         'date': datetime.utcnow()}
        for key, value in self.__dict__.iteritems():
            try:
                print key, value.__class__, sys.getsizeof(value), asizeof(value)
            except:
                print key

            if type(value) not in ban_types:
                sim_blob[key] = value
            elif type(value) == type(np.array([12])):
                print key, value.nbytes
                if value.ndim == 1:
                    sim_blob[key] = value.tolist()
                elif value.size == 1:
                    sim_blob[key] = value
                else:
                    print 'adding to the db obj', key, type(value), value.shape, value.size
                    if value.size< 1600000:
                        sim_blob[key] = value
                    else:
                        print 'nevermind, ', key, ' is too big'

        print 'before: ', asizeof(sim_blob), ' ', sys.getsizeof(sim_blob), ' ', flatsize(sim_blob)
        cutlist = ['pos_i','kx','ky','pdf_x','df_x','v']
        #cutlist = ['pos_i','kx','ky','pdf_x','df_x','v']
        #cutlist = ['t_stop','dx']

        keep_list = ['nz','nx','alpha_c']
        
        for k in cutlist:
            try:
                print 'k: ',(sim_blob[k])
                sim_blob.pop(k, None)
            except:
                print 'not found'
        #new_blob = {}
    # for k in keep_list:
        #     new_blob.pop(k, None)

        # print 'look for bugs'
        # for k in sim_blob:
        #     print 'trying: ',k#,sim_blob[k]
        #     new_blob[k] = sim_blob[k]
        #     ser_dict = self.serialize(new_blob)
        #     #print sys.getsizeof(new_dict), asizeof(ser_dict), ' ', asizesof(sim_blob)
        #     print BSON.encode(ser_dict).__sizeof__()

        # exit()

        #sim_blob = new_blob

        print 'after: ', asizeof(sim_blob), ' ', sys.getsizeof(sim_blob), ' ', flatsize(sim_blob)
        print 'coarse_x' in sim_blob
        print 'df_x' in sim_blob
        #print type(sim_blob['coarse_x']), sim_blob['coarse_x'][0].__class__
        #exit()
    
        ser_dict = self.serialize(sim_blob)
        # for elem in ser_dict:
        #     print elem,': ',type(ser_dict[elem])m
        #print ser_dict
        print sys.getsizeof(ser_dict), asizeof(ser_dict), ' ', asizesof(sim_blob)
        print BSON.encode(ser_dict).__sizeof__()
        #exit()
        c = MongoClient(host=server)
        db = c.new_database
        print 'db info'
        print db.command('collstats', 'alpha_runs')
        alpha_runs = db.alpha_runs
        alpha_runs.ensure_index('md5', unique=True, dropDups=True)
        #db.things.ensureIndex({'source_references.key' : 1}, {unique : true, dropDups : true})

        try:
            #print ser_dict.keys()
            alpha_runs.insert(ser_dict, db)
        except mongoErr.DuplicateKeyError:
            alpha_runs.remove({'path': ser_dict['path']})
            alpha_runs.remove({'md5': ser_dict['md5']})
            
            #alpha_runs.remove({ "$and" :[{'path': ser_dict['path']},{'md5': ser_dict['md5']}]})
            alpha_runs.insert(ser_dict, db)
            print 'Duplicate run not adding to db, but updating'
