from blob_info import Blob2D
# from pb_corral import LinRes
# from pb_nonlinear import NLinResDraw
# from pb_transport import Transport

import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
#from matplotlib.ticker import FuncFormatter
#from matplotlib.ticker import ScalarFormatter 

# from reportlab.platypus import *
# from reportlab.lib.styles import getSampleStyleSheet
# from reportlab.rl_config import defaultPageSize
# from reportlab.lib.unitas import inch
# from reportlab.graphics.charts.linecharts import HorizontalLineChart
# from reportlab.graphics.shapes import Drawing
# from reportlab.graphics.charts.lineplots import LinePlot
# from reportlab.graphics.widgets.markers import makeMarker
# from reportlab.lib import colors

# from replab_x_vs_y import RL_Plot
#for movie making
from multiprocessing import Queue,Pool
import multiprocessing
import subprocess


#the function of this class to extend the base blob class 
#and create the images to be rendered by the present class


#all the draw methods require display object to write to, they do 
#however prodive own canvas/figure




class BlobDraw(Blob2D):
    def __init__(self,data,meta=None,fast_center=True):
        Blob2D.__init__(self,data,meta=meta,fast_center=fast_center)

    def basic_2Dview(self,t_i = None):
        fig = plt.figure()
        canvas = fig.add_subplot(121)
        canvas.imshow(np.log(self.raw_data[-1,:,:]+.1),aspect='auto')
        canvas = fig.add_subplot(122)
        fftimg = np.fft.fft2(self.raw_data[-1,:,:])
        canvas.imshow(np.log(np.real(fftimg*fftimg.conj())[0:self.nx/8,0:self.ny/16.0]),
                      aspect='auto',interpolation='bicubic')
        plt.close(fig)
        
        return fig

    def basic_view(self,t_i = None):
        fig = plt.figure()
        canvas = fig.add_subplot(121)
        canvas.imshow(np.log(self.raw_data[-1,:,:]+.1),aspect='auto')
        canvas = fig.add_subplot(122)
        fftimg = np.fft.fft2(self.raw_data[-1,:,:])
        plt.close(fig)
        
        return fig
    
    def draw_svd(self,comp='V',nmodes=6,center=False,zoom = False,
                 interp='none'):
        #compute svd is not present . . later

        fig = plt.figure()
        try:
            y = self.svd[comp]
        except:
            y = self.svd['V']

        for i in  xrange(nmodes):
            canvas = fig.add_subplot(2,3,i+1)
            
            # if comp=='V':
            #     y_i = y[i,:,:]
            # else:
            #     y_i = y[:,i]
            y_i = y[i]
        
            nx,ny = y_i.shape
            if center:
                y_i = np.roll(np.roll(y_i,nx/2,axis=0),ny/2,axis=1)
            if zoom and not center:
                y_i = y_i[0:nx/8.0,0:ny/8.]
            if zoom and center:
                y_i = y_i[nx/2-nx/8:nx/2.0+nx/8,
                          ny/2-ny/16:ny/2.0+ny/16,]
    
            canvas.imshow(np.real(y_i),aspect='auto',interpolation= interp)
        
        plt.close(fig)
        return fig

            
    # def summary()
    # def plot_svdnew(self):
        
    #     pp = PdfPages('svd.pdf')
    #     fig = plt.figure()
    #     canvas = fig.add_subplot(121)
    #     canvas.imshow(np.log(self.center_data[-1,:,:]+.1),aspect='auto')
    #     # fig.savefig(pp,format='pdf')
        
    #     # fig = plt.figure()
    #     canvas = fig.add_subplot(122)
    #     fftimg = np.fft.fft2(self.center_data[-1,:,:])
        
    #     canvas.imshow(np.log(np.real(fftimg*fftimg.conj())[0:self.nx/8,0:self.ny/16.0]),aspect='auto',interpolation='bicubic')
    #     fig.savefig(pp,format='pdf')
        
    #     plt.close(fig)
        
    #     fig = plt.figure()
    
    #     for i in  xrange(6):
    #         canvas = fig.add_subplot(2,3,i+1)
    #         canvas.imshow(self.svd['V'][i,:,:],aspect='auto')
        
    #     fig.savefig(pp,format='pdf')
    #     plt.close(fig)
        
    #     fig = plt.figure()
    
    #     for i in  xrange(6):
    #         canvas = fig.add_subplot(2,3,i+1)
    #         print self.svd['R'][0].shape
    #         canvas.imshow(np.roll(np.roll(np.real(self.svd['R'][i]),self.ny/2,axis=1),self.nx/2, axis=0),aspect='auto')
        
    #     fig.savefig(pp,format='pdf')

    #     fig = plt.figure()
    
    #     for i in  xrange(6):
    #         canvas = fig.add_subplot(2,3,i+1)
    #         print self.svd['R'][0].shape
    #         #canvas.imshow(np.real(np.roll(np.roll(np.fft.fft2(self.svd['R'][i]),self.ny/2,axis=1),self.nx/2, axis=0)),aspect='auto')
    #         canvas.imshow(np.real(np.fft.fft2(self.svd['R'][i]))[0:self.nx/8,0:self.ny/16],aspect='auto',interpolation='bicubic')
    #     fig.savefig(pp,format='pdf')


    #     fig = plt.figure()
    #     for i in  xrange(0,12,2):
    #         canvas = fig.add_subplot(3,4,i+1)

    #         canvas.plot(self.svd['P'][i/2][0:self.nx/4.0,0])
            
    #         canvas = fig.add_subplot(3,4,i+2)
    #         canvas.plot(self.svd['P'][i/2][0,0:self.ny/32.0],'r')
    
    #     fig.savefig(pp,format='pdf') 

    #     #### we will extract some parameters from  the autocorrection 
    #     def fit_2gauss(x,y):
    #         from scipy.optimize import curve_fit
    #         from scipy.signal import argrelextrema  
    #     #find local maxima
    #         max_indx = argrelextrema(y, np.greater,
    #                                  order = 4,mode='wrap')
    #         print 'ymax: ', y[max_indx],max_indx
        
    #         #print max_indx[0][1]
    #         if len(max_indx) > 1:
    #             p=[np.max(y),1,max_indx[0][1],0,1]
    #         else:
    #             p=[np.max(y),1,len(x)/2,0,1]

    #         def func(x, a, b, x1,c,d):
    #             return a*np.exp(-b * (x)**2) + c*np.exp(-d * (x-x1)**2)  
            
        
    #         popt, pcov= curve_fit(func, x[0:len(x)/4], y[0:len(x)/4],p0=p)
    #         return popt,pcov, popt[0]*np.exp(-popt[1] * (x)**2) + popt[3]*np.exp(-popt[4] * (x-popt[2])**2) 

    #     fig = plt.figure()
    #     for i in  xrange(0,12,2):
    #         canvas = fig.add_subplot(3,4,i+1)

    #         print 'sizes: ', self.x[:,0].shape, self.svd['R'][0][0,:].shape
    #         canvas.plot(self.x[:,0],self.svd['R'][i/2][:,0])
    #         try:
    #             popt,pcov,estimatex =  fit_2gauss(self.x[:,0], np.real(self.svd['R'][i/2][:,0]))
    #             canvas.plot(self.x[:,0],estimatex)
    #         except:
    #             print 'no fit'

    #         canvas = fig.add_subplot(3,4,i+2)
    #         canvas.plot(self.y[0,:],self.svd['R'][i/2][0,:],'r')
    #         try:
    #             popt,pcov,estimatey =  fit_2gauss(self.y[0,:], np.real(self.svd['R'][i/2][0,:]))
    #             canvas.plot(self.y[0,:],estimatey)
    #         except:
    #             print 'no fit'
            
    #         #canvas.imshow(np.roll(np.roll(np.real(self.svd['R'][i]),self.ny/2,axis=1),self.nx/2, axis=0),aspect='auto')
        
    #     fig.savefig(pp,format='pdf')


    #     plt.close(fig)
        
    #     fig = plt.figure()
        
    #     print self.svd['U'].shape
    #     for i in  xrange(6):
    #         canvas = fig.add_subplot(2,3,i+1)
    #         canvas.plot(self.svd['U'][:,i])
    
    #     fig.savefig(pp,format='pdf')
    #     plt.close(fig)

    #     fig = plt.figure()
    #     canvas = fig.add_subplot(121)
    #     canvas.plot(np.cumsum(self.svd['s']**2)/(self.svd['s']**2).sum())
    #     canvas.set_yscale('log')
    #     canvas.set_xscale('log')

    #     canvas = fig.add_subplot(122)
    #     canvas.plot(self.svd['s'] * self.svd['U'][-1,:])
    #     #canvas.set_yscale('log')

    #     fig.savefig(pp,format='pdf')
    #     plt.close(fig)
        
    #     pp.close()
    #     return 0   
