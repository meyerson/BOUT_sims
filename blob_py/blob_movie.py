import os, sys, inspect
import sqlite3 as sql
import pickle as pkl
from datetime import datetime

HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')
utc = datetime.utcnow()


sys.path.append('/usr/local/pylib')
sys.path.append(HOME+'/lib/python')
sys.path.append(HOME+'/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib/boutdata')
sys.path.append(BOUT_TOP+'/tools/pylib/boututils')
sys.path.append(BOUT_TOP+'/tools/pylib/post_bout')

from blob_info import Blob2D
from blob_draw import BlobDraw
from frame import Frame

# from pb_corral import LinRes
# from pb_nonlinear import NLinResDraw
# from pb_transport import Transport

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages

try:
    import matplotlib.artist as artist 
    import matplotlib.ticker as ticker
    import matplotlib.animation as animation  
    from matplotlib.lines import Line2D
    import mpl_toolkits.axisartist as axisartist
    import pylab
    from matplotlib import rc
    from matplotlib import rc
except:
    print "Can't import what you need"

#this calls will create figure object whose data will be update at each frame and
#the result saved

class BlobMovie(BlobDraw):
    def __init__(self,data,meta=None,fast_center=True):
        BlobDraw.__init__(self,data,meta=meta,fast_center=fast_center)
        if meta is not None:
            for k, v in meta.items():
                setattr(self, k, v)
    
    def movieframe(self,frames,moviename='output',
                   fast=True,bk=None,outline=True,
                   t_array=None,encoder='mencoder',fps=5.0):
        
        if fast:
            dpi = 100
        else:
            dpi = 240
            
        lin_formatter = ticker.ScalarFormatter()
        lin_formatter.set_powerlimits((-2, 2))
        
        font = {'family' : 'normal',
                'weight' : 'normal',
                'size'   : 4}
        
        axes = {'linewidth': .5}
        tickset ={'markeredgewidth': .25}
 
        
        rc('font', **font)
        rc('axes',**axes)
        rc('lines',**tickset)
        
        plt.tick_params(axis='both',direction='in',which='both')
    
        jet = plt.get_cmap('jet',2000) 


        
        fig = plt.figure()

        nrow = 2
        ncol = 2
        
        for i,elem in enumerate(frames):
            elem.render(fig,nrow,ncol,i)
        
        def update_img(t):
            for elem in frames:
                elem.update()
    
        ani = animation.FuncAnimation(fig,update_img,self.nt-2)    
        ani.save(moviename+'.mp4',writer=encoder,dpi=dpi,bitrate=20000,fps=5)

        plt.close(fig)
       
        return 0
   
    def movie(self,data2=None,moviename='output',norm=False,
              overcontour=True,aspect='auto',meta=None,
              cache='/tmp/',hd=False,nlevels = 9,removeZero=True,
              t_array=None,outline=True,bk=None,fps=5.0,fast=True,
              encoder='mencoder'):
        
        data = self.raw_data
        size = data.shape
        ndims = len(size)
        nt,nx,nz = size  

        if norm: 
            data_n = normalize(data)
        else:
            data_n = data
        
        if data2 != None:
            data_c = data2
        else:
            data_c = data
            
        if norm:
            data_c = normalize(data_c)

        
        if fast:
            dpi = 100
        else:
            dpi = 240

        amp = self.amp

        kx_max = self.kx_max
        kz_max = self.ky_max
        fft_data = self.fft
        power = self.power
        acorr = self.acorr
        kx = self.kx
        kz = self.ky
        k = np.sqrt(kx**2 + kz**2)

        axhandle = []

        lin_formatter = ticker.ScalarFormatter()
        lin_formatter.set_powerlimits((-2, 2))
        
        font = {'family' : 'normal',
                'weight' : 'normal',
                'size'   : 4}
        
        axes = {'linewidth': .5}
        tickset ={'markeredgewidth': .25}
 
        
        rc('font', **font)
        rc('axes',**axes)
        rc('lines',**tickset)
        
        plt.tick_params(axis='both',direction='in',which='both')
    
        jet = plt.get_cmap('jet',2000) 
        
        
        #this is where we start using the frame object to make it happen
        
    
        frm_data = Frame(self.raw_data,meta={'mask':True,'data_c':data_c})
        frm_power = Frame(self.power,meta={'zoom':True,'coords':'k','center':True})
        frm_amp = Frame(self.amp)
        frm_cm = Frame(np.array(self.xmoment[1]))
        frm_acorr = Frame(self.acorr,meta={'center':True})
        
        fig = plt.figure()

        frm_data.render(fig,231)
        frm_power.render(fig,232)
        frm_amp.render(fig,233)
        frm_cm.render(fig,234)
        frm_acorr.render(fig,235)

        def update_img(t):
            print t
            frm_data.update()
            frm_power.update()
            frm_amp.update()
            frm_cm.update()
            frm_acorr.update()
        
        ani = animation.FuncAnimation(fig,update_img,self.nt-2)    
        ani.save(moviename+'.mp4',writer=encoder,dpi=dpi,bitrate=20000,fps=5)
        
        #plt.savefig('output.pdf',dpi=dpi) 
        
        # frames[1].render(fig,221)
        # plt.savefig('output.pdf',dpi=dpi)
        
        plt.close(fig)
       
        return 0   
    
