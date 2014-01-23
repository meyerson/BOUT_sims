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

    colors = ['b','g','r','k','c','y','k','b','g','r','c','m','y','k']
    styles = ['s','^']

    
 
    handles=[]
    
    labeler={'264_mu=1e-3_HD':'BOUT++, Ra='+r'$10^6$','264_mu=1e-2_HD':'BOUT++, Ra='+r'$10^4$','264_mu=1e-1_HD':'Ra = '+r'$10^1$','Ra1e4':'Garcia ref., Ra='+r'$10^4$','Ra1e6':'Garcia ref., Ra='+r'$10^6$'}
    markers = {'264_mu=1e-3_HD':'-','264_mu=1e-2_HD':'--','264_mu=1e-1_HD':'Ra = '+r'$10^1$','Ra1e4':'s','Ra1e6':'^'}

    
    for i,elem in enumerate(cm): 
        j = np.int(np.random.rand(1))
        label = labeler[str(elem['label'])]
        #label = 'adsfasf'
        # xcanvas.plot(elem['t'],elem['x'],colors[i]+styles[j],
        #              label=label,alpha = .5,markersize=6)
        # xcanvas.plot(elem['t'],elem['x'],colors[i],alpha=1)
        #xcanvas.annotate(elem['label'], (1.02*elem['t'][-1],elem['x'][-1]), xytext=None, xycoords='data',
         #               textcoords='data', arrowprops=None,fontsize = 10)
        
        if 'ref' in elem.keys():
            mstyle =  markers[str(elem['label'])]
            #print colors[i]+markers[i]
            xcanvas.plot(elem['t'],elem['x'],colors[i]+mstyle,alpha = .4,
                         label=label,linewidth=4,markersize=8)
            xcanvas.plot(elem['t'],elem['x'],colors[i],alpha = .5,
                         linewidth=1)
        else:
            
            xcanvas.plot(elem['t'],elem['x'],colors[i]+markers[str(elem['label'])],
                         alpha=.8,label=label,linewidth=3)


    handles, labels = xcanvas.get_legend_handles_labels()  
    leg = xcanvas.legend(handles,labels,ncol=2,loc='best',prop={'size':14},fancybox=True) 
    leg.get_frame().set_alpha(0.9)
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
        label = labeler[str(elem['label'])]
        #label = 'adsfasf'
        print elem.keys()
        try:
            vx = elem['V']
            x = elem['t_v']
            print vx.shape
        except:
            vx = np.gradient(np.array(elem['x']))/elem['dt']
            x = elem['t']
        print x.shape,vx.shape,elem['label']
        
        if 'ref' in elem.keys():
            
            mstyle =  markers[str(elem['label'])]

            vcanvas.plot(x,vx,colors[i]+mstyle,alpha = .4,
                         label=label,linewidth=4,markersize=8)
            vcanvas.plot(x,vx,colors[i],alpha = .5,
                         linewidth=1)
        else:
        
            vcanvas.plot(x,vx,colors[i]+markers[str(elem['label'])],alpha=1,label=label,linewidth=3)

        # if 'ref' in elem.keys():
        #     vcanvas.plot(elem['t'],elem['x'],colors[i]+styles[j],alpha = .2,
        #                  label=label,linewidth=4,markersize=8)
        #     xcanvas.plot(elem['t'],elem['x'],colors[i],alpha = .5,
        #                  linewidth=1)
        # else:
            
        #     xcanvas.plot(elem['t'],elem['x'],colors[i]+'--',alpha=.8,label=label,
        #                  linewidth=3)
       # vcanvas.annotate(elem['label'], (1.02*x[-1],vx[-1]), xytext=None, xycoords='data', textcoords='data', arrowprops=None,fontsize = 10)  

    handles, labels = xcanvas.get_legend_handles_labels()  
    leg = vcanvas.legend(handles,labels,ncol=2,loc='best',prop={'size':14},fancybox=True) 
    leg.get_frame().set_alpha(0.8)

    vcanvas.set_title('center of mass velocity')
    #vcanvas.set_ylabel('CM - x')
    vcanvas.set_ylabel(r'$\frac{V_x}{\omega_{ci} \rho_{ci}}$',fontsize=20,rotation='horizontal')
    vcanvas.set_xlabel(r'$\frac{t}{\tau_{ci}}$',fontsize=20)
    vcanvas.grid(True,linestyle='-',color='.75')
    
    if compare_png_v is not None:
        vcanvas.imshow(compare_png_v,extent=[0,25,0,1],aspect='auto')
   
    #figX.savefig(pp, format='pdf')
    #figV.savefig(pp, format='pdf')
    # figM.savefig(pp, format='pdf')
    # figC.savefig(pp, format='pdf')


########
# Plot amplitude history
#########
    handles=[]
    
    for i,elem in enumerate(cm):     
        j = np.int(np.round(np.random.rand(1)))
        #label = str(elem['label'])
        label = labeler[str(elem['label'])]
        if 'max' in elem.keys():
            val_max = np.array(elem['max'])
            x = np.array(elem['t'])
        elif 'Namp' in elem.keys():
            val_max = np.array(elem['Namp'])
            x = np.array(elem['t_a'])
        else:
            val_max = np.array(elem['t'])*0
            x = np.array(elem['t'])
        if 'ref' in elem.keys():
            print x.shape,val_max.shape
            mstyle =  markers[str(elem['label'])]
            maxcanvas.plot(x,val_max,colors[i]+mstyle,alpha = .4,
                           label=label,linewidth=4,markersize=8)
            maxcanvas.plot(x,val_max,colors[i],alpha = .5,
                         linewidth=1)
        else:
            maxcanvas.plot(x,val_max,colors[i]+markers[str(elem['label'])],alpha=.8,label=label, linewidth=3)
            # maxcanvas.annotate(elem['label'], (1.02*elem['t'][-1],val_max[-1]), xytext=None, xycoords='data', textcoords='data', arrowprops=None,fontsize = 10)  
        print elem.keys()
    handles, labels = maxcanvas.get_legend_handles_labels()  
    print labels
    
    #exit()
    leg = maxcanvas.legend(handles,labels,ncol=2,loc='best',prop={'size':14},fancybox=True) 
    leg.get_frame().set_alpha(0.8)

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
    # handles=[]
    
    # for i,elem in enumerate(cm):     
    #     j = np.int(np.round(np.random.rand(1)))
    #     #label = labeler[str(elem['label'])]
    #     label = str(elem['label'])
    #     if 'net' in elem.keys():
    #         val_net = np.gradient(np.array(elem['net']))/np.array(elem['net'])
    #     else:
    #         val_net = np.array(elem['x'])*0

    #     cons_canvas.plot(elem['t'],val_net,colors[i]+styles[j],alpha = .5,
    #                  label=label,markersize=4)
    #     cons_canvas.plot(elem['t'],val_net,colors[i],alpha=.8)
    #     #cons_canvas.annotate(elem['label'], (1.02*elem['t'][-1],val_net[-1]), xytext=None, xycoords='data', textcoords='data', arrowprops=None,fontsize = 10)  

    # handles, labels = cons_canvas.get_legend_handles_labels()  
    # leg = cons_canvas.legend(handles,labels,ncol=2,loc='best',prop={'size':14},fancybox=True) 
    # leg.get_frame().set_alpha(0.8)

    # cons_canvas.set_title('total density')
    # #vcanvas.set_ylabel('CM - x')
    # cons_canvas.set_ylabel(r'$\frac{n}{n(t=0)}$',fontsize=20,rotation='horizontal')
    # cons_canvas.set_xlabel(r'$\frac{t}{\tau_{ci}}$',fontsize=20)
    # cons_canvas.grid(True,linestyle='-',color='.75')
    
    # if compare_png_max is not None:
    #     cons_canvas.imshow(compare_png_max,extent=[0,25,0,1],aspect='auto')


    if ownpage is True:
        figX.savefig(pp, format='pdf')
        figV.savefig(pp, format='pdf',bbox_inches='tight',pad_inches=.5)
        figM.savefig(pp, format='pdf',bbox_inches='tight',pad_inches=.5)
        figC.savefig(pp, format='pdf')
        plt.close(figX) 
        plt.close(figV)
        plt.close(figM)
        plt.close(figC)
