import numpy as np
import math
import numpy as np
import matplotlib.pyplot as plt
import hashlib

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
#from scipy.signal import argrelextrema  
from scipy.optimize import curve_fit
import sys, os, gc


from boutdata import collect
import subprocess 
from scipy import interpolate
from scipy import signal

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

def expfall(x,y0,l,c):
    return y0*np.exp(-x/l) + c

def fit_lambda(y,x,appendto=None):
    return curve_fit(expfall, x.flatten(), y.flatten())
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
        # returns a a list with containing [wave, scales, freqs, coi, fft, fftfreqs] 
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

def expfall(x,y0,l):
     return y0*np.exp(-x/l)  

class field_info(object):
    def __init__(self,path,field="n",meta=None,fast_center=True,get_Xc=True,
                 get_lambda=True,debug=False):
        self.path=path
        self.field=field
        
        f=open(path+'/BOUT.log.0', 'r')
        self.md5 = hashlib.md5(f.read()).hexdigest()
        f.close()
        

        if meta is not None:
            for k, v in meta.items():
                setattr(self, k, v)
        
        self.nz = np.squeeze(collect("MZ",xind=[0,0],path=path,info=False))-1
        nxpe =  np.squeeze(collect("NXPE",xind=[0,0],path=path,info=False))
        mxg = np.squeeze(collect("MXG",xind=[0,0],path=path,info=False))
        
        
        self.mxsub = np.squeeze(collect("MXSUB",xind=[0,0],path=path,info=False)) #gaurds
        self.nx = self.mxsub*nxpe + 2*mxg #including gaurds
        nx,ny,nz = self.nx,self.nz,self.nz
        
        self.dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
        self.dy = np.squeeze(collect("dz",path=path,xind=[0,0]))
        self.zmax = np.squeeze(collect("ZMAX",path=path))

        self.time = np.squeeze(collect("t_array",path=path,xind=[0,0]))
        self.nt = len(self.time)

        defaults = {'x0':0,'y0':0}   

        for key,val in defaults.items():
            if not hasattr(self,key):
                setattr(self, key, val)

        if debug:        
            t_chunk = 10
            t_stop  = 100#np.max(self.nt)
        else:
            t_chunk = 40
            t_stop  = np.max(self.nt)

        self.t_chunk  = t_chunk  
        t1 = 0
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
        #print nx,ny,nz,nxpe,mxg, self.mxsub , 2*mxg 
        self.k =  np.mgrid[kxmin:kxmax:dkx,kymin:kymax:dky]
        self.kx = np.arange(0,nx)*self.dkx
        self.ky = np.arange(0,ny)*self.dky

        #print nx,nx*self.mxsub
        self.pos_i = np.mgrid[0:nx:1,0:ny:1]

        #print 't: ',t_stop,t1,t2
        
          
        #print pos.shape,pos[0][xstart:xstop,:].shape,n.shape,nave[-1]
        #popt, pcov= fit_lambda(n-nave[-1],pos[0][xstart:xstop,:])
        self.x = np.mgrid[xmin:xmax:self.dx]
        
        try:
            a = np.squeeze(collect("alpha",path=path))
            a_ave = self.a_ave = np.average(a,axis=1)
            xstart= np.int(np.round(np.mean(np.where(abs(a_ave-.3*a_ave.max()) < .1* a_ave.mean()))))
            xstop = np.int(xstart+ nx/2.)
            if xstop > nx:
                xstop = nx -1;
        except:
            xstart= np.int(nx/3.)            
            xstop = np.int(xstart+ nx/2.)



        self.lam  = []
        while t2<t_stop:
            data = np.squeeze(collect("n",tind=[t1,t2],path=path,info=False))
            
            nave = (data.mean(axis = 0)).mean(axis=1)[xstart:xstop]

            est_lam = (self.x[xstop]-self.x[xstart])/(np.log(nave[0]/nave[-1]))
            p0=[nave[0],est_lam]
            popt, pcov= curve_fit(expfall,self.x[xstart:xstop],nave,p0=p0)
           
            self.lam.append([t1 + .5*(t2-t1),popt[0]])

            print popt, p0
            # print data.shape,nx
            nt = data.shape[0]
            for t in xrange(nt):
                #print data.shape, self.pos_i.shape
                moments(data[t,:,:],self.pos_i[0,:,:],appendto=self.xmoment)
                moments(data[t,:,:],self.pos_i[1,:,:],appendto=self.ymoment) 
                #print t
            t1 = t2+1
            t2 = np.min([t1+t_chunk-1,t_stop+1])
            print t1,t2,self.dx

     
                
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
        nscales = 40
        maxScale = (nt/dt)/2.0
        minScale = 10*dt
        
        dj = (np.log(maxScale) - np.log(minScale))/(nscales* np.log(2))
        

        if (get_lambda):
            self.lambdafit = []
            self.lam = []
            #lamda_x={1:[],2:[]}
            for t in xrange(nt):
                popt, pcov= fit_lambda(data[t,:,:],self.pos[0,:,:])
                self.lam.append(popt[1])
                self.lambdafit.append(popt[0]*np.exp(-self.pos[0,:,:]/popt[1])
                                     + popt[2])
                
            #print t,':',lamda_x[-1])

        self.fftdata = np.fft.fft2(data)
        self.power = self.fftdata.conj() * self.fftdata
        
        
        
       # self.kxy = np.mgrid[kxmin:xmax:nx*complex(0,1),ymin:ymax:ny*complex(0,1)]
        
        kz_max,kx_max = ny,nx

        kz = (2.0*np.pi/(1.0*dy))*np.linspace(0,kz_max-1,kz_max)
        kz = np.repeat(kz,kx_max)
        kz = np.transpose(kz.reshape(kz_max,kx_max))
        
        kx =(2.0*np.pi/(1.0*dx))*np.linspace(0,kx_max-1,kx_max)
        kx = np.repeat(kx,kz_max)
        kx = kx.reshape(kx_max,kz_max)
    
        k = np.sqrt(kx**2 + kz**2)
 

    def plot_wave(self):
        pp = PdfPages('wave.pdf')
        fig = plt.figure()
        canvas = fig.add_subplot(111)
        wave =  self.freq_vs_t_kx[0]
        power = np.real((abs(wave)) ** 2)   
        canvas.imshow((power+.1),aspect='auto',interpolation='none')
        fig.savefig(pp,format='pdf')

        fig = plt.figure()
        canvas = fig.add_subplot(221)
        nt = np.round(self.nt*(4./5.))
        ymid = np.round(self.ny/2)
        xslice = self.center_data[nt,:,ymid]
        dt = 1.0
        nscales = 60
        maxScale = (self.nx/dt)/2.0
        minScale = 2*dt
        dj = (np.log(maxScale) - np.log(minScale))/(nscales* np.log(2))

        canvas.plot(xslice)
        
        canvas = fig.add_subplot(222)
        xscale_vs_y = wvlt(xslice,s0=minScale,dj=dj,J=nscales)
        wave =  xscale_vs_y[0]
        power = np.real((abs(wave)) ** 2)   
        canvas.imshow((power+.1),aspect='auto',interpolation='none')
        
        canvas = fig.add_subplot(223)
        xslice = self.center_data[0,:,ymid]
        xscale_vs_y = wvlt(xslice,s0=minScale,dj=dj,J=nscales)
        canvas.plot(xslice)
        canvas = fig.add_subplot(224)
        xscale_vs_y = wvlt(xslice,s0=minScale,dj=dj,J=nscales)
        wave =  xscale_vs_y[0]
        power = np.real((abs(wave)) ** 2)   
        canvas.imshow((power+.1),aspect='auto',interpolation='none')
        fig.savefig(pp,format='pdf')
        fig = plt.figure()


        canvas = fig.add_subplot(221)
        nt = np.round(self.nt*(4./5.))
        xmid = np.round(self.nx/2)
        yslice = self.center_data[nt,xmid,:]
        dt = 1.0
        nscales = 60
        maxScale = (self.ny/dt)/2.0
        minScale = 10*dt
        dj = (np.log(maxScale) - np.log(minScale))/(nscales* np.log(2))

        canvas.plot(yslice)
        
        canvas = fig.add_subplot(222)
        yscale_vs_y = wvlt(yslice,s0=minScale,dj=dj,J=nscales)
        wave =  yscale_vs_y[0]
        power = np.real((abs(wave)) ** 2)   
        canvas.imshow((power+.1),aspect='auto',interpolation='none')
        
        canvas = fig.add_subplot(223)
        yslice = self.center_data[0,xmid,:]
        yscale_vs_y = wvlt(yslice,s0=minScale,dj=dj,J=nscales)
        canvas.plot(yslice)
        canvas = fig.add_subplot(224)
        yscale_vs_y = wvlt(yslice,s0=minScale,dj=dj,J=nscales)
        wave =  yscale_vs_y[0]
        power = np.real((abs(wave)) ** 2)   
        canvas.imshow((power+.1),aspect='auto',interpolation='none')
        
        #canvas.imshow(self.raw_data[nt,:,:],aspect='auto',interpolation='none')
      
        plt.close(fig)
        fig.savefig(pp,format='pdf')


        pp.close()
        return 0

    def plot_grad(self):
        nt = self.nt
    
        pp = PdfPages('grad.pdf')
        fig = plt.figure()
        canvas = fig.add_subplot(111)
        # wave =  self.freq_vs_t_kx[0]
        # power = np.real((abs(wave)) ** 2)   
        #canvas.imshow((self.center_data[self.nt/2,self.xmoment[1][np.round(nt/2)]:,
         #                              self.ymoment[1][np.round(nt/2)]:]),aspect='auto',interpolate='none')
        ymin = self.ymoment[1][nt/2] - 2*np.sqrt(self.ymoment[2][nt/2])
        ymax = self.ymoment[1][nt/2] + 2*np.sqrt(self.ymoment[2][nt/2])
        xmin = self.xmoment[1][nt/2] - 2*np.sqrt(self.xmoment[2][nt/2])
        xmax = self.xmoment[1][nt/2] + 2*np.sqrt(self.xmoment[2][nt/2])

        canvas.imshow((self.raw_data[self.nt/2,xmin:xmax,ymin:ymax]),aspect='auto')
        
        #

        plt.close(fig)
        fig.savefig(pp,format='pdf')
        
        import kPyWavelet as wavelet #cwt 
        morlet = wavelet.Morlet(3)
        morl = morlet.psi(np.mgrid[-10:11])
        fig = plt.figure()
        canvas = fig.add_subplot(221)

        scalex = 4 #in pixels
        kern = wave_kern(scalex)
        kern2 = wave_kern(2*scalex)
        canvas.imshow(np.real(kern2))
        canvas = fig.add_subplot(222)
        canvas.plot(np.real(kern[:,3*scalex+1]))
        #canvas = fig.add_subplot(222)
        canvas.plot(np.real(kern2[:,6*scalex+1]))
        plt.close(fig)

        fig.savefig(pp,format='pdf')
        
        # fig = plt.figure()
        # canvas = fig.add_subplot(111)
        
        # result3 = signal.fftconvolve(self.raw_data[nt/2,:,:],kern)
        # canvas.imshow(np.real(result3),aspect='auto')
        # plt.close(fig)
        # fig.savefig(pp,format='pdf')

        fig = plt.figure()
        canvas = fig.add_subplot(111)
        sizes = [2,5,10,20,40,60]
        proj = []
        
        for size in sizes: 
            kern = wave_kern(size)
            result3 = signal.fftconvolve(self.raw_data[nt/2,:,:],kern)
            norm = signal.fftconvolve(self.raw_data[nt/2,:,:]*0 + 1.0,kern)
            result = (result3/norm)
            result = np.sqrt(result.conj() * result).sum()
            
         
        canvas.plot(sizes,proj)
   
        plt.close(fig)
        fig.savefig(pp,format='pdf')

        pp.close()
        return 0
    
    def plot_svd(self):
        
        pp = PdfPages('svd.pdf')
        fig = plt.figure()
        canvas = fig.add_subplot(121)
        canvas.imshow(np.log(self.center_data[-1,:,:]+.1),aspect='auto')
        # fig.savefig(pp,format='pdf')
        
        # fig = plt.figure()
        canvas = fig.add_subplot(122)
        fftimg = np.fft.fft2(self.center_data[-1,:,:])
        
        canvas.imshow(np.log(np.real(fftimg*fftimg.conj())[0:self.nx/8,0:self.ny/16.0]),aspect='auto',interpolation='bicubic')
        fig.savefig(pp,format='pdf')
        
        plt.close(fig)
        
        fig = plt.figure()
    
        for i in  xrange(6):
            canvas = fig.add_subplot(2,3,i+1)
            canvas.imshow(self.svd['V'][i,:,:],aspect='auto')
        
        fig.savefig(pp,format='pdf')
        plt.close(fig)
        
        fig = plt.figure()
    
        for i in  xrange(6):
            canvas = fig.add_subplot(2,3,i+1)
            #print self.svd['R'][0].shape
            canvas.imshow(np.roll(np.roll(np.real(self.svd['R'][i]),self.ny/2,axis=1),self.nx/2, axis=0),aspect='auto')
        
        fig.savefig(pp,format='pdf')

        fig = plt.figure()
    
        for i in  xrange(6):
            canvas = fig.add_subplot(2,3,i+1)
            print self.svd['R'][0].shape
            #canvas.imshow(np.real(np.roll(np.roll(np.fft.fft2(self.svd['R'][i]),self.ny/2,axis=1),self.nx/2, axis=0)),aspect='auto')
            canvas.imshow(np.real(np.fft.fft2(self.svd['R'][i]))[0:self.nx/8,0:self.ny/16],aspect='auto',interpolation='bicubic')
        fig.savefig(pp,format='pdf')


        fig = plt.figure()
        for i in  xrange(0,12,2):
            canvas = fig.add_subplot(3,4,i+1)

            canvas.plot(self.svd['P'][i/2][0:self.nx/4.0,0])
            
            canvas = fig.add_subplot(3,4,i+2)
            canvas.plot(self.svd['P'][i/2][0,0:self.ny/32.0],'r')
    
        fig.savefig(pp,format='pdf') 

        #### we will extract some parameters from  the autocorrection 
        def fit_2gauss(x,y):
            from scipy.optimize import curve_fit
            from scipy.signal import argrelextrema  
        #find local maxima
            max_indx = argrelextrema(y, np.greater,
                                     order = 4,mode='wrap')
            print 'ymax: ', y[max_indx],max_indx
        
            #print max_indx[0][1]
            if len(max_indx) > 1:
                p=[np.max(y),1,max_indx[0][1],0,1]
            else:
                p=[np.max(y),1,len(x)/2,0,1]

            def func(x, a, b, x1,c,d):
                return a*np.exp(-b * (x)**2) + c*np.exp(-d * (x-x1)**2)  
            
        
            popt, pcov= curve_fit(func, x[0:len(x)/4], y[0:len(x)/4],p0=p)
            return popt,pcov, popt[0]*np.exp(-popt[1] * (x)**2) + popt[3]*np.exp(-popt[4] * (x-popt[2])**2) 

        fig = plt.figure()
        for i in  xrange(0,12,2):
            canvas = fig.add_subplot(3,4,i+1)

            print 'sizes: ', self.x[:,0].shape, self.svd['R'][0][0,:].shape
            canvas.plot(self.x[:,0],self.svd['R'][i/2][:,0])
            try:
                popt,pcov,estimatex =  fit_2gauss(self.x[:,0], np.real(self.svd['R'][i/2][:,0]))
                canvas.plot(self.x[:,0],estimatex)
            except:
                print 'no fit'

            canvas = fig.add_subplot(3,4,i+2)
            canvas.plot(self.y[0,:],self.svd['R'][i/2][0,:],'r')
            try:
                popt,pcov,estimatey =  fit_2gauss(self.y[0,:], np.real(self.svd['R'][i/2][0,:]))
                canvas.plot(self.y[0,:],estimatey)
            except:
                print 'no fit'
            
            #canvas.imshow(np.roll(np.roll(np.real(self.svd['R'][i]),self.ny/2,axis=1),self.nx/2, axis=0),aspect='auto')
        
        fig.savefig(pp,format='pdf')


        plt.close(fig)
        
        fig = plt.figure()
        
        print self.svd['U'].shape
        for i in  xrange(6):
            canvas = fig.add_subplot(2,3,i+1)
            canvas.plot(self.svd['U'][:,i])
    
        fig.savefig(pp,format='pdf')
        plt.close(fig)

        fig = plt.figure()
        canvas = fig.add_subplot(121)
        canvas.plot(np.cumsum(self.svd['s']**2)/(self.svd['s']**2).sum())
        canvas.set_yscale('log')
        canvas.set_xscale('log')

        canvas = fig.add_subplot(122)
        canvas.plot(self.svd['s'] * self.svd['U'][-1,:])
        #canvas.set_yscale('log')

        fig.savefig(pp,format='pdf')
        plt.close(fig)
        
        pp.close()
        return 0
    
  

    def present(self):
        return 0
        
    def save(self):
        return 0
        
        
#blob1 = Blob2D(
