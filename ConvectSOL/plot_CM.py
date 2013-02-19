import numpy as np
import math
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker

import subprocess 

def CM_mass(data,meta=None,label=None):

    
    
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
    

    pos = np.mgrid[xmin:xmax:nx*complex(0,1),ymin:ymax:ny*complex(0,1)]

    x = []#np.zeros(nt)
    y = []

    print pos[0,::]

    for t in xrange(nt):
        print np.mean( pos[1,::])
        x.append(np.sum(data[t,:,:] * pos[0,::])/np.sum(data[t,:,:]))
        y.append(np.sum(data[t,:,:] * pos[1,::])/np.sum(data[t,:,:]))
        
    cm = {'x':x,'y':y,'t':dt*np.arange(nt),'label':label,'dt':dt} 
    #print cm
    print ymin
    print ymax
    print xmin
    print xmax
    return cm

def present(cm,pp,xcanvas=None,vcanvas=None,compare_png_x=None,
            compare_png_v=None):
    
    xownpage  = False
    vownpage = False
    
    if xcanvas is None:
        xownpage = True
        figX = plt.figure()
        xcanvas = figX.add_subplot(1,1,1) 
        figX.subplots_adjust(bottom=0.14)

    if vcanvas is None:
        vownpage = True     
        figV = plt.figure()
        vcanvas = figV.add_subplot(1,1,1)

    colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
    styles = ['s','^']

    
    #fig.subplots_adjust(top=0.80)
    #fig.subplots_adjust(right=0.83)
    #fig.subplots_adjust(left=0.17)
    #adj = fig.subplots_adjust(hspace=0.4,wspace=0.4)
    for i,elem in enumerate(cm): 
        j = np.int(np.random.rand(1))
        xcanvas.plot(elem['t'],elem['x'],colors[i]+styles[j],alpha = .2,markersize=2)
        xcanvas.plot(elem['t'],elem['x'],colors[i],alpha=.7)
        xcanvas.annotate(elem['label'], (1.02*elem['t'][-1],elem['x'][-1]), xytext=None, xycoords='data',
                        textcoords='data', arrowprops=None,fontsize = 10)
    xcanvas.set_title(' center of mass')
    #xcanvas.ylabel('CM - x')
    xcanvas.set_ylabel(r'$\frac{x}{\rho_{ci}}$',fontsize=20,rotation='horizontal')
    xcanvas.set_xlabel(r'$\frac{t}{\tau_{ci}}$',fontsize=20)
   
    xcanvas.set_xscale('linear')
    xcanvas.grid(True,linestyle='-',color='.75')
    
    if compare_png_x is not None:
        # compare_png[:,:,3] = compare_png[:,:,3]/2
        # compare_png[:,:,0] = compare_png[:,:,0]*2
        # compare_png[:,:,1] = compare_png[:,:,1]*2

        xcanvas.imshow(compare_png_x,extent=[0,25,0,15],aspect='auto')

   
    
    # fig = plt.figure()
    # fig.subplots_adjust(bottom=0.14)
    for i,elem in enumerate(cm):     
        j = np.int(np.round(np.random.rand(1)))
        vx = np.gradient(np.array(elem['x']))/elem['dt']
        vcanvas.plot(elem['t'],vx,colors[i]+styles[j],alpha = .4,markersize=4)
        vcanvas.plot(elem['t'],vx,colors[i],alpha=.2)
        vcanvas.annotate(elem['label'], (1.02*elem['t'][-1],vx[-1]), xytext=None, xycoords='data', textcoords='data', arrowprops=None,fontsize = 10)  
    vcanvas.set_title(' center of mass velocity')
    #vcanvas.set_ylabel('CM - x')
    vcanvas.set_ylabel(r'$\frac{V_x}{\omega_{ci} \rho_{ci}}$',fontsize=20,rotation='horizontal')
    vcanvas.set_xlabel(r'$\frac{t}{\tau_{ci}}$',fontsize=20)
    # #plt.xlabel('$\frac{x}{\rho_{ci}}$',fontsize=20)
    # plt.xscale('linear')
    vcanvas.grid(True,linestyle='-',color='.75')
    
    if compare_png_v is not None:
        # compare_png[:,:,3] = compare_png[:,:,3]/2
        # compare_png[:,:,0] = compare_png[:,:,0]*2
        # compare_png[:,:,1] = compare_png[:,:,1]*2
        vcanvas.imshow(compare_png_v,extent=[0,25,0,1],aspect='auto')
    
    # fig.savefig(pp, format='pdf')
    # plt.close() 
    if xownpage is True:
        figX.savefig(pp, format='pdf')
        figV.savefig(pp, format='pdf')
        plt.close(figX) 
        plt.close(figV)
