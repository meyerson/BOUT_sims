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
    max_val =[]
    net = []
    print pos[0,::]

    for t in xrange(nt):
        
        x.append(np.sum(data[t,:,:] * pos[0,::])/np.sum(data[t,:,:]))
        y.append(np.sum(data[t,:,:] * pos[1,::])/np.sum(data[t,:,:]))
        max_val.append(np.max(data[t,:,:]))
        net.append(np.sum(data[t,:,:]))

        print t,':',np.max(data[t,:,:])
    cm = {'x':x,'y':y,'t':dt*np.arange(nt),'label':label,
          'dt':dt,'max':max_val,'net':net} 
    #print cm
    print ymin
    print ymax
    print xmin
    print xmax
    return cm

def present(cm,pp,xcanvas=None,vcanvas=None,maxcanvas=None,cons_canvas=None,
            compare_png_x=None,compare_png_v=None,
            compare_png_max=None):
    
    ownpage  = False
    #vownpage = False
    
    # for elem in canvas_stack:
    #     ownpage = True
    #     figX = plt.figure()
    #     xcanvas = figX.add_subplot(1,1,1) 
    #     figX.subplots_adjust(bottom=0.14)

    # def set_canvas(canvas):
    #     fig = plt.figure()
    #     xcanvas = figX.add_subplot(1,1,1) 
    #     figX.subplots_adjust(bottom=0.14)
        
    if xcanvas is None:
        ownpage = True
        figX = plt.figure()
        xcanvas = figX.add_subplot(1,1,1) 
        figX.subplots_adjust(bottom=0.14)

    if vcanvas is None:
        ownpage = True     
        figV = plt.figure()
        vcanvas = figV.add_subplot(1,1,1)

    if maxcanvas is None:
        ownpage = True
        figM = plt.figure()
        maxcanvas = figM.add_subplot(1,1,1) 
        figM.subplots_adjust(bottom=0.14)
        
    if cons_canvas is None:
        ownpage = True
        figC = plt.figure()
        cons_canvas = figC.add_subplot(1,1,1) 
        figC.subplots_adjust(bottom=0.14)

    colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
    styles = ['s','^']

    
 
    handles=[]

    for i,elem in enumerate(cm): 
        j = np.int(np.random.rand(1))
        label = str(elem['label'])
        xcanvas.plot(elem['t'],elem['x'],colors[i]+styles[j],
                     label=label,alpha = .3,markersize=3)
        xcanvas.plot(elem['t'],elem['x'],colors[i],alpha=1)
        #xcanvas.annotate(elem['label'], (1.02*elem['t'][-1],elem['x'][-1]), xytext=None, xycoords='data',
         #               textcoords='data', arrowprops=None,fontsize = 10)
        

    handles, labels = xcanvas.get_legend_handles_labels()  
    leg = xcanvas.legend(handles,labels,ncol=2,loc='best',prop={'size':8},fancybox=True) 
    leg.get_frame().set_alpha(0.5)
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
###########
    handles=[]

    for i,elem in enumerate(cm):     
        j = np.int(np.round(np.random.rand(1)))
        label = str(elem['label'])
        vx = np.gradient(np.array(elem['x']))/elem['dt']
        vcanvas.plot(elem['t'],vx,colors[i]+styles[j],alpha = .3,
                     label=label,markersize=4)
        vcanvas.plot(elem['t'],vx,colors[i],alpha=1)
        vcanvas.annotate(elem['label'], (1.02*elem['t'][-1],vx[-1]), xytext=None, xycoords='data', textcoords='data', arrowprops=None,fontsize = 10)  

    handles, labels = xcanvas.get_legend_handles_labels()  
    leg = vcanvas.legend(handles,labels,ncol=2,loc='best',prop={'size':8},fancybox=True) 
    leg.get_frame().set_alpha(0.3)

    vcanvas.set_title(' center of mass velocity')
    #vcanvas.set_ylabel('CM - x')
    vcanvas.set_ylabel(r'$\frac{V_x}{\omega_{ci} \rho_{ci}}$',fontsize=20,rotation='horizontal')
    vcanvas.set_xlabel(r'$\frac{t}{\tau_{ci}}$',fontsize=20)
    vcanvas.grid(True,linestyle='-',color='.75')
    
    if compare_png_v is not None:
        vcanvas.imshow(compare_png_v,extent=[0,25,0,1],aspect='auto')
   


########
# Plot amplitude history
#########
    handles=[]
    
    for i,elem in enumerate(cm):     
        j = np.int(np.round(np.random.rand(1)))
        label = str(elem['label'])
        if 'max' in elem.keys():
            val_max = np.array(elem['max'])
        else:
            val_max = np.array(elem['x'])*0

        maxcanvas.plot(elem['t'],val_max,colors[i]+styles[j],alpha = .4,
                     label=label,markersize=4)
        maxcanvas.plot(elem['t'],val_max,colors[i],alpha=.8)
        maxcanvas.annotate(elem['label'], (1.02*elem['t'][-1],val_max[-1]), xytext=None, xycoords='data', textcoords='data', arrowprops=None,fontsize = 10)  

    handles, labels = maxcanvas.get_legend_handles_labels()  
    leg = maxcanvas.legend(handles,labels,ncol=2,loc='best',prop={'size':8},fancybox=True) 
    leg.get_frame().set_alpha(0.3)

    maxcanvas.set_title('maximum amp')
    #vcanvas.set_ylabel('CM - x')
    maxcanvas.set_ylabel(r'$\frac{n}{n(t=0)}$',fontsize=20,rotation='horizontal')
    maxcanvas.set_xlabel(r'$\frac{t}{\tau_{ci}}$',fontsize=20)
    maxcanvas.grid(True,linestyle='-',color='.75')
    
    if compare_png_max is not None:
        maxcanvas.imshow(compare_png_max,extent=[0,25,0,1],aspect='auto')

######
# Conservation
#####
    handles=[]
    
    for i,elem in enumerate(cm):     
        j = np.int(np.round(np.random.rand(1)))
        label = str(elem['label'])
        if 'net' in elem.keys():
            val_net = np.gradient(np.array(elem['net']))/np.array(elem['net'])
        else:
            val_net = np.array(elem['x'])*0

        cons_canvas.plot(elem['t'],val_net,colors[i]+styles[j],alpha = .4,
                     label=label,markersize=4)
        cons_canvas.plot(elem['t'],val_net,colors[i],alpha=.8)
        cons_canvas.annotate(elem['label'], (1.02*elem['t'][-1],val_net[-1]), xytext=None, xycoords='data', textcoords='data', arrowprops=None,fontsize = 10)  

    handles, labels = cons_canvas.get_legend_handles_labels()  
    leg = cons_canvas.legend(handles,labels,ncol=2,loc='best',prop={'size':8},fancybox=True) 
    leg.get_frame().set_alpha(0.3)

    cons_canvas.set_title('total density')
    #vcanvas.set_ylabel('CM - x')
    cons_canvas.set_ylabel(r'$\frac{n}{n(t=0)}$',fontsize=20,rotation='horizontal')
    cons_canvas.set_xlabel(r'$\frac{t}{\tau_{ci}}$',fontsize=20)
    cons_canvas.grid(True,linestyle='-',color='.75')
    
    if compare_png_max is not None:
        cons_canvas.imshow(compare_png_max,extent=[0,25,0,1],aspect='auto')


    if ownpage is True:
        figX.savefig(pp, format='pdf')
        figV.savefig(pp, format='pdf')
        figM.savefig(pp, format='pdf')
        figC.savefig(pp, format='pdf')
        plt.close(figX) 
        plt.close(figV)
        plt.close(figM)
        plt.close(figC)
