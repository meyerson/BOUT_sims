import numpy as np
import math

import matplotlib.pyplot as plt
import hashlib
from pymongo import Connection, MongoClient
from pymongo import errors as mongoErr
from datetime import datetime

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
            t1 = np.int(self.nt*.95)
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
        

        
        for x in xrange(nx):
            #print moving_min_n[x]
            pdf_x.append(np.mgrid[moving_min_n[x]-.01:moving_max_n[x] +.01:complex(0,300)])
            df_x.append(np.mgrid[moving_min[x]-.01:moving_max[x] +.01:complex(0,300)])
            
        
        
        self.pdf_x = pdf_x
        self.df_x = df_x
        pdf_y=nx*[None]
        df_y=nx*[None]
        pdf_y_net = 0
        df_y_net = 0

        

        for i in range(nx):
            #print "some unique object %d" % ( i, )
            pdf_y[i]= 0*pdf_x[i]
            df_y[i]= 0*df_x[i]

        
        
        xstart = np.min(np.where((a_ave-.1*a_ave.max()>0)) )
        xstop = np.int(xstart+ nx/2.)

        if xstop > nx:
            xstop = nx -1;                      
        
        #########################################    
        #now lets actually get the statistical info 
        ###############    ########    ##########

       
    
        j = 0 #course time chunk counter
        print 't2,t_stop:',t2,t_stop
        self.nave = []
        self.flux = []

        while t2<t_stop:
            print 'computing statistical measures'
            #try:
            #everything but the raw "data" variable are based on n not log(n)
            (nave,dn),(flux,flux_max,d_flux),data = self.dens_and_flux(t1,t2) #coarse in t, fine in x 
            self.save_linlam(np.log(nave),np.log(dn))
            # sys.stdout = mystdout = StringIO()
            self.nave.append([nave,dn])
            self.flux.append([flux,flux_max,d_flux])

            print "data.shape: ", data.shape
            if t1 > t_stop/2:
                j = j + 1
                for x in xrange(nx):
                    if not fast:
                        p = np.polyfit(np.arange(t2-t1+1),np.squeeze(data[:,x,:]).mean(axis=1),1)
                        nave2d = p[1]+p[0]*np.arange(t2-t1+1)
                        if self.log_n:
                            cond = np.exp(data[:,x,:]) - np.exp(np.reshape(np.repeat(nave2d,ny),[t2-t1+1,ny])) #nt * ny - nt*ny 
                        else:
                            cond = data[:,x,:] - np.reshape(np.repeat(nave2d,ny),[t2-t1+1,ny])
                        
                    else:
                        if self.log_n:
                            cond = np.exp(data[:,x,:]) - nave #nt * ny - nt*ny
                        else:
                            cond = data[:,x,:] - nave #nt *

                    nrms = cond.std()
                    datamax = data[:,x,:].max()
                    
                    #local_pdf_x

                    #if  datamax is not 0 and np.isfinite(datamax):
                    if np.isfinite(datamax):  
                       
                        #print data[:,x,:].shape,nave[x]#pdf[x].shape
                        dfadd  = np.squeeze((np.histogram(cond,bins = 300,normed=False,
                                                            range=(df_x[x].min(),df_x[x].max()))[0]))
                        
                        pdfadd =  np.squeeze((np.histogram(cond/nrms,bins = 300,normed=False,
                                                           range=(pdf_x[x].min(),pdf_x[x].max()))[0]))
                        # pdfadd =  np.squeeze((np.histogram(cond/nrms,bins = 300,normed=False,
                        #                                    range=(-2.0,2.0))[0]))
                      
                        if np.isfinite(np.sum(pdfadd)):
                            pdf_y[x] = pdf_y[x] + pdfadd
                            df_y[x] = pdf_y[x] + dfadd
                            if (abs(x-xstart)<nx/20.):
                                pdf_y_net = pdf_y_net+pdfadd
                                df_y_net = df_y_net + pdfadd
                                #print pdf_y[x]
                                # print np.min(pdf_x[np.int(xstart)])
                                # exit()
                    else:
                        ##pdf_x = 0*pdf_x
                        ##pdf_yy = 0*pdf_yy
                        pass
               
            self.pdf_y = pdf_y
            self.pdf_y_net = pdf_y_net
            self.df_y = pdf_y
            self.df_y_net = pdf_y_net

            self.t.append(t1 + .5*(t2-t1))

            t1 = t2+1
            t2 = np.min([t1+t_chunk-1,t_stop+1])

            


        self.nave = np.array(self.nave)
        self.flux = np.array(self.flux)

        print 'almost done'
        for x in xrange(nx):    
            #print x,pdf_x[x].mean(),nave_net[x].mean(),np.mean(pdf_x[x])
            if np.mean(nrms_net[x,:]) != 0:
                #self.pdf_x[x] = (pdf_x[x] - nave_net[x])/(nrms_net[x])
                self.pdf_y[x] = pdf_y[x]/sum(pdf_y[x])
            else:
                self.pdf_x[x] = pdf_x[x-1]

        if hasattr(self,'pdf_y_net'):
            print 'sum(pdf_y_net): ',sum(pdf_y_net)
            self.pdf_y_net = pdf_y_net/np.float(sum(pdf_y_net))

            self.df_y = self.pdf_y
        
        print self.pdf_x[self.xstart]
        #exit()    
    def save_linlam(self,nave,dn):
        a_ave = self.a_ave
        nx = self.nx
        #xstart = np.min(np.where((a_ave-.1*a_ave.max()>0)) )
        #xstop = np.int(xstart+ nx/2.)
        xstart = self.xstart
        #xstop = self.xstop
        xstop = np.int(xstart+ nx/2.)
        if xstop > nx:
            xstop = nx -1;
        
        if self.log_n:
            est_lam = (self.x[xstop]-self.x[xstart])/(nave[xstart] - nave[xstop])    
            p0=[nave[xstart],est_lam]
            popt, pcov= curve_fit(linearfall,self.x[xstart:xstop],nave[xstart:xstop],p0=p0)
            
        else:
            est_lam = (self.x[xstop]-self.x[xstart])/(np.log(nave[xstart]/nave[xstop]))
            p0=[np.log(nave[xstart]),est_lam]
            popt, pcov= curve_fit(linearfall,self.x[xstart:xstop],np.log(nave[xstart:xstop]),p0=p0)
        self.linlam.append(popt[1])
           

    def dens_and_flux(self,t1,t2): #fine in x, smooth in t profiles
        sys.stdout = mystdout = StringIO()
        path = self.path
        dy = self.dy
        if self.fast:
            data = np.squeeze(collect("n",tind=[t1,t1+1],path=path,info=False))
            v = np.squeeze(np.gradient(np.squeeze(collect("phi",tind=[t1,t1+1],path=path,info=False)))[2])/dy
        else:
            try:
                data = np.squeeze(collect("n",tind=[t1,t2],path=path,info=False))
                v = np.squeeze(np.gradient(np.squeeze(collect("phi",tind=[t1,t2],path=path,info=False)))[2])/dy
            except:
                data = np.squeeze(collect("n",tind=[t1,t2-1],path=path,info=False))
                v = np.squeeze(np.gradient(np.squeeze(collect("phi",tind=[t1,t2-1],path=path,info=False)))[2])/dy

        if self.log_n:
            flux = np.exp(data)*v
        else:
            flux = data*v
        sys.stdout = old_stdout
        if self.log_n:
            n = (np.exp(data).mean(axis = 0)) #ave in time, still x,z
            n = (n.mean(axis = 1))
           # n = np.log(n) #yes I take the log again
        else:
            n = (data.mean(axis = 0)) #ave in time, still x,z
            n = (n.mean(axis = 1))

        dshape =  data.shape
        
        if self.log_n:
            temp = np.transpose(data,[1,0,2]) #seems ok given that we will fit an exp. to this
        else:
            temp = np.transpose(data,[1,0,2])
        
        temp = np.reshape(temp,[dshape[1],np.prod(dshape)/dshape[1]]) 
        n_sigma = temp.std(axis=1)

        
        temp = np.transpose(flux,[1,0,2])
        temp = np.reshape(temp,[dshape[1],np.prod(dshape)/dshape[1]])
        flux_sigma = temp.std(axis=1)
        flux_max = temp.max(axis=1)
        flux = flux.mean(axis = 0)
        flux = flux.mean(axis = 1)

        

        return (n,n_sigma),(flux,flux_sigma,flux_max),data


    def compute_profile(self): #coarse in x, fine in t
        self.linlam  = []
        
        self.lam =[]
        self.t=[]
        self.pdf_x =[]
        self.pdf_y = []
        self.df_x =[]
        self.df_y = []

        nrms = []
        nave = []
        nmax = []
        nmin = []
        coarse_x =[]
       


        #get mean profiles for t > nt/2
        xchunk = 10
        nx = self.nx
        nt = self.nt
        path = self.path
        xpos = np.arange(xchunk)*nx/xchunk
        xpos = list(xpos)
        xpos.append(nx-1)

        #build 
            
        self.nave = [] 

        for x in np.array(xpos):
            print "x: ", x/nx,nt
            coarse_x.append(x)
            print x,round(nt/2)
            sys.stdout = mystdout = StringIO()
            data =  np.squeeze(collect("n",xind=[int(x),int(x)+5],tind=[nt/2,nt-1],path=path,info=False))
            # if self.log_n:
            #     data = (np.squeeze(collect("n",xind=[int(x),int(x)+5],tind=[nt/2,nt-1],path=path,info=False)))
            if self.log_n:
                data = np.exp(data)

            sys.stdout = old_stdout
            #nrms.append(data.std(axis))
            
           
            dshape =  data.shape
            temp = np.reshape(data,[dshape[0],np.prod(dshape)/dshape[0]])
           
            nrms.append(temp.std(axis=1))

            
            try:
                nave.append((data.mean(axis=1)).mean(axis=1))
                nmax.append((data.max(axis=1)).max(axis=1))
                nmin.append((data.min(axis=1)).min(axis=1))
           
            except:
                nave.append(data.mean(axis=1))
                nmax.append(data.max(axis=1))
                nmin.append(data.min(axis=1))
                #nrms.append(data.std(axis=1))
                
        nave = np.array(nave) # a time resolved ave
        nrms = np.array(nrms)
        nmax = np.array(nmax)
        nmin = np.array(nmin)

        print "nave.shape: ",nave.shape,nmax.shape,
        print len(coarse_x),len((np.arange(nt - nt/2)))

        
        tarray = np.array(np.arange(nt - nt/2))
       
        self.nave_net  = np.transpose(nave)
        self.nrms  = np.transpose(nrms)
        self.nmin  = np.transpose(nmin)
        self.nmax  = np.transpose(nmax)
        self.coarse_x = coarse_x #nx x nt

        nmin_net =[] 
        nmax_net = []
        nrms_net= []
        nave_net = []
        coarse_x = self.coarse_x

        xnew = np.arange(0,nx)
        for tt in tarray:
            #print tt
            
            f = interpolate.interp1d(coarse_x,nave[:,tt],kind='linear')
            nave_net.append(f(xnew))
            
            f = interpolate.interp1d(coarse_x,nrms[:,tt],kind='linear')
            nrms_net.append(f(xnew))

            f = interpolate.interp1d(coarse_x,nmin[:,tt],kind='linear')
            nmin_net.append(f(xnew))
            
            f = interpolate.interp1d(coarse_x,nmax[:,tt],kind='linear')
            nmax_net.append(f(xnew))

        nave_net = np.array(nave_net)
        nrms_net = np.array(nrms_net)
        nmin_net = np.array(nmin_net)
        nmax_net = np.array(nmax_net)
        
        moving_min_n = np.min((nmin_net - nave_net)/nrms_net,axis=0)
        moving_max_n = np.max((nmax_net - nave_net)/nrms_net,axis=0)
        # print moving_min_n
        # exit()
        moving_min = np.min((nmin_net - nave_net),axis=0)
        moving_max = np.max((nmax_net - nave_net),axis=0)
        
        return moving_min_n, moving_max_n, moving_min, moving_max, nave_net, nrms_net ,nmin_net,nmax_net

    def serialize(self,input_dict):
        ser_dict = {}
    
        #print 'pdf_x', input_dict['pdf_y'][-1]
        for key,value in input_dict.iteritems():
            print '#######################'
            print key,type(value)
      #if type(value) not in copy_types:
            if is_numeric_paranoid(value):# or type(value) in copy_types:
                print key, 'np.array',type(value)#, value.shape
                if 'all' in dir(value):
                    #print key,value.size,value.ndim
                    #print key, value.ndim,type(value), value.size
                    if (value.ndim ==1 and value.size ==1) or (value.ndim == 0):
                        #print value
                        ser_dict[key] = np.asscalar(value)
                    else:
                        print 'packing into a binary'
                        print value.size
                        ser_dict[key] = Binary(pickle.dumps(value,protocol=2))
                else:
                    ser_dict[key] = value

            elif type(value) in copy_types:
                ser_dict[key] = value

            elif type(value) is types.ListType:
                print "list: ", len(value)
                if len(value) == 0:
                    ser_dict[key] = value
                elif  type(value[0]) == type(np.array([12.23])):
                    print value[0].size,len(value)
                    #print key, type(value[0]),type(value[0]) == type(np.array([12.23]))
                    ser_dict[key] =  Binary(pickle.dumps(value,protocol=2))
                else:
                    ser_dict[key] = value

            elif type(value) is types.DictType:
                ser_dict[key] = value

            else:
                print key, 'serializing',value,type(value),value.size()
                ser_dict[key] = Binary(pickle.dumps(value,protocol=2))

        # for elem in ser_dict:
        #    print elem, ser_dict[elem].__class__
           
        # exit()
        #print sim_blob.keys()
        #print 'self.pdf_x: ',self.pdf_x[-1]
        

        return ser_dict
    
    def to_db(self,server='beavis.ph.utexas.edu'):

        #print self.nave_net.shape
        #exit()
       # print dir(self)
        #exit()
        print self.pdf_x[self.xstart]
        #exit()
        sim_blob={"author": "Dmitry",
                 "tags": ["fusion", "in", "50 years"],
                 "date": datetime.utcnow()}
        for key,value in self.__dict__.iteritems():
            #print key, type(value)
            #print 'take data from sim to a big dictionary'
            if type(value) not in ban_types:
                #print key, type(value)
                sim_blob[key] = value
            else:
                if type(value) == type(np.array([12])):
                    #print key, value.shape
                    if value.ndim == 1:
                        sim_blob[key] = value.tolist()
                        #print 'convert ', key
                    elif value.size == 1:
                        #print key
                        sim_blob[key] = value
                    else:          
                        print 'adding to the db obj',key, type(value),value.shape, value.size
                        if value.size < 1600000:
                            sim_blob[key] = value
                        else:
                            print 'nevermind, ', key, ' is too big'
                            
            
      
        #a very primitive way to look for bad keys

        #cutlist = ['zmax', 'linlam', 't_chunk', 'nave', 'y0', 'author', 'fast','nx', 'field', 'nz','location', 'nt', 'ymoment', 'lam', 'tags', 'xmoment', 'max_val','ky_max', 'dx','dy', 'date', 'path','x0', 'mxsub', 'kx_max','dkx', 'dky','y', 'x', 'ky', 'kx','nave_net','nrms','nmin','pos_i','t','df_y','df_x','coarse_x' ]
        #cutlist = ['pdf_x','pdf_y','coarse_x','df_x','df_y']
        #print 'nave',sim_blob['nave'].size,sim_blob['nave'][0].size,sim_blob['nave'].shape
        #print 'flux',sim_blob['flux'].__class__,sim_blob['flux'][0].size,sim_blob['flux'].shape
        #exit()
        #sim_blob['flux'] = sim_blob['nave']
        cutlist = ['df_x','df_y','coarse_x']
        #cutlist = ['coarse_x' ]
        for k in cutlist:
            print k,type(sim_blob[k])#,type(sim_blob[k][0]),is_numeric_paranoid(sim_blob[k])
            sim_blob.pop(k, None)
        
        #print sim_blob
        print 'sim blob: ',sim_blob.keys()
        
        #print self.pdf_x
        #print sim_blob['pdf_x']

        #exit()
       # for elem in sim_blob:
        #    print elem, sim_blob[elem].__class__
        #     try:
        #         print sim_blob[elem].shape
        #     except:
        #         pass
        # exit()  
        #sim_blob = {'test':1}
        #exit()
        ser_dict = self.serialize(sim_blob)
        #ser_dict = self.serialize({'test':1})
        #exit()
        c = MongoClient(host=server)
        #db = c.test_database # 
        db = c.new_database
        #c.drop_database(db)
        alpha_runs = db.alpha_runs
        alpha_runs.ensure_index('md5',unique=True,dropDups=True)#,{'unique': True, 'dropDups': True})

        try:
            print ser_dict.keys()
            #for elem in ser_dict:
            #    print elem, ser_dict[elem].__class__
            alpha_runs.insert(ser_dict,db)
        except mongoErr.DuplicateKeyError:
            
            #print sim_blob.keys()
            #print 'the hash: ',ser_dict['md5'],sim_blob['md5'],ser_dict.keys()

            ser_dict.pop("_id",None) #rip off the id 

            alpha_runs.update({'md5':ser_dict['md5']}, {"$set":ser_dict})#, upsert=False)
            print 'Duplicate run not adding to db, but updating'
            

  
class Turbulence(object):
    def __init__(self,data,meta=None,pp=None,fast_center=True,get_Xc=True,
                 get_lambda=True):
        
        self.raw_data = data
        self.md5 = hashlib.md5(data).hexdigest() #unique blob IDl

        self.fft = np.fft.fft2(data)
        self.power = self.fft.conj()*self.fft
        self.acorr = np.fft.ifft2(self.power)

        self.shape = data.shape
        self.ndim = data.ndim
        
        nt,nx,ny = data.shape
        self.nt, self.nx, self.ny = nt,nx,ny
        
        if meta is not None:
            for k, v in meta.items():
                setattr(self, k, v)
        
        defaults = {'dx':1,'x0':0,'dy':1,'y0':0}   

        for key,val in defaults.items():
            if not hasattr(self,key):
                print 'setting: ',key,val
                setattr(self, key, val)      
      
        #amplitude
        self.amp = abs(data).max(1).max(1)    
        
#compute the first several moments
        self.xmoment = {1:[],2:[]}
        self.ymoment = {1:[],2:[]}
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
        dkx, dky = 1.0, 1.0

        self.k =  np.mgrid[kxmin:kxmax:dkx,kymin:kymax:dky]
        self.kx,self.ky = self.k
        
        #print xmin,xmax,ymin,ymax
        #self.pos = np.mgrid[xmin:xmax:nx*complex(0,1),ymin:ymax:ny*complex(0,1)]
        self.pos = np.mgrid[xmin:xmax:dx,ymin:ymax:dy]
        self.x, self.y = self.pos

        self.center_data = 0.0*data
        self.pos_i = np.mgrid[0:nx:1,0:ny:1]
        #self.pos_i = np.mgrid[0:nx-1:nx*complex(0,1),0:ny-1:ny*complex(0,1)]
        old_points = np.transpose(np.array([self.pos_i[0,...].flatten(),self.pos_i[1,...].flatten()]))
        
        grad = np.gradient(data)
        self.grad = np.sqrt(grad[1]**2 + grad[2]**2)
        
        
        for t in xrange(nt):
        
            moments(data[t,:,:],self.pos_i[0,:,:],appendto=self.xmoment)
            moments(data[t,:,:],self.pos_i[1,:,:],appendto=self.ymoment)
            self.max_val.append(np.max(data[t,:,:]))
            #print "xmomnt: ",self.xmoment[1][-1]
            
            #we can find the index displacement needed to center the density
            #field
            if np.isfinite(self.xmoment[1][-1]):
                xshift = nx/2.0 - self.xmoment[1][-1]
            else:
                xshift = nx/2.0
            
            if np.isfinite(self.ymoment[1][-1]):
                yshift = ny/2.0 - self.ymoment[1][-1]
            else:
                yshift = ny/2.0
 
          
            #
           
            #phase = -2*np.pi*(xshift)*pos_i[0,...]/nx #- 2*np.pi*(yshift)*pos_i[1,...]/ny 
            
            #fftphaseshift = np.exp(np.vectorize(complex)(0,phase))
            
            
            #self.center_data[t,:,:]= np.real(np.fft.ifft2(fftphaseshift*(np.fft.fft2(data[t,:,:]))))
            
            #fancy interp method to go the blob cm frame
            if fast_center:
                self.center_data[t,:,:]= np.roll(data[t,:,:],np.int(xshift),axis=0)
                self.center_data[t,:,:]= np.roll(self.center_data[t,:,:],np.int(yshift),axis=1)    
            else:
                pos_newx = np.mod(self.pos_i[0,...]-xshift,nx-1).flatten()
                pos_newy = np.mod(self.pos_i[1,...]-yshift,ny-1).flatten()
                new_points = np.transpose(np.array([pos_newx,pos_newy]))
           
                centered = interpolate.griddata(old_points,
                                           data[t,...].flatten(), 
                                           new_points,method='linear',fill_value = .1)
                self.center_data[t,:,:] = centered.reshape(nx,ny)

            
        print data.shape,self.nx*self.ny
        U,s,V = np.linalg.svd(self.center_data.reshape(self.nt,self.nx*self.ny),
                              full_matrices=False)

        self.svd = {'V':V.reshape(self.nt,self.nx,self.ny),
                    's':s,'U':np.transpose(U),'R':[],'P':[],'cutoff':0}
        
        #find the cutoff that will keep the field w/in .1%
        tol = (np.cumsum(self.svd['s']**2)/np.sum(self.svd['s']**2)-.999)**2
        max_indx = argrelextrema(tol, np.less,
                                 order = 2,mode='wrap')

        self.svd['cutoff'] = max_indx
        self.svd['tol'] = tol


        #let computere autocorrelation of each V (topos)
        
        #for v_i in xrange(max_indx[0]):
        for v_i in xrange(3):
            fftv = np.fft.fft2(self.svd['V'][v_i,:,:])
            power = fftv.conj() * fftv
            self.svd['R'].append(np.fft.ifft2(power))#autocorrelation
            self.svd['P'].append(np.abs(power))
        
        #print U.shape, s.shape,V.reshape(self.nt,self.nx,self.ny)
        
            
        
        dt = 1.0
 
