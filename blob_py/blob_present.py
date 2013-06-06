from blob_info import Blob2D
from blob_draw import BlobDraw

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

class BlobPresent(BlobDraw):
    def __init__(self,data,meta=None,fast_center=True):
        BlobDraw.__init__(self,data,meta=meta,fast_center=fast_center)
    
    def plot_svdnew(self):
        
        pp = PdfPages('svd.pdf')
        fig = self.basic_2Dview()
        fig.savefig(pp,format='pdf')
        
 
        
        #fig = plt.figure()
    
        # for i in  xrange(6):
        #     canvas = fig.add_subplot(2,3,i+1)
        #     canvas.imshow(self.svd['V'][i,:,:],aspect='auto')
        
        fig = self.draw_svd()
        fig.savefig(pp,format='pdf')

        fig = self.draw_svd(comp='R',center=True)
        fig.savefig(pp,format='pdf')


        fig = self.draw_svd(comp='P',zoom=True,center=True,interp='bicubic')
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
