import os,sys
boutpath = os.environ['BOUT_TOP']
pylibpath = boutpath+'/tools/pylib'
pbpath = pylibpath+'/post_bout'
boutdatapath = pylibpath+'/boutdata'
boututilpath = pylibpath+'/boututils'
allpath = [boutpath,pylibpath,pbpath,boutdatapath,boututilpath]

[sys.path.append(elem) for elem in allpath]
#print sys.path
#from ordereddict import OrderedDict
from scipy.interpolate import interp2d,interp1d
from boutdata import collect
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
#from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
#from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
from boutdata import collect
from boututils import showdata
import numpy as np
from scipy import ndimage
from frame import Frame
from matplotlib import patches
try:
    path = sys.argv[1]
except:
    path = '/tmp/SOLblob/data_blob_nlog_min'
def fast2Dplot(pp,data,title=None,xlabel=None,ylabel=None,addcurve=None):
    
    fig = plt.figure()
    
    sm = fig.add_subplot(1,1,1)
    
    #sm.imshow(np.flipud(np.rot90(data)),aspect='auto',interpolation='none',origin='lower')
    sm.imshow(data,aspect='auto',interpolation='none',origin='lower')
    sm.contour(data,colors='k',linestyles='solid',alpha=.5)
    #sm.imshow(data, interpolation='none', aspect='auto',origin='lower')
    if addcurve != None:
        boundary = sm.plot(addcurve['x'],addcurve['y'])
        plt.setp(boundary, linewidth=4,alpha = 1.0,color=addcurve['color'])
  
    if title != None:
        sm.set_title(title)
    if xlabel != None:
        sm.set_xlabel('x index')
    if ylabel != None:
        sm.set_ylabel('y index')

    fig.savefig(pp, format='pdf')


def aveplot(pp,data,title='chaotic '+r"$\alpha \pm \sigma$",x_label=r'$\rho_s$',
            label=r"$\alpha$",format='pdf',dx=1.0):
    fig, sm = plt.subplots(1)
    
    alpha = data.mean(axis=1)
    sigma = data.std(axis=1)

    r = np.arange(data.shape[0])*dx
#line, = sm.plot(alpha, color='blue', lw=2)
    sm.plot(r,alpha, lw=2, label=label, color='blue')
    sm.fill_between(r,alpha+sigma, alpha-sigma, facecolor='yellow', alpha=0.5)
    sm.set_title(title)
    #sm.set_yscale('symlog')
    sm.legend(loc='upper left')
    sm.set_xlabel(x_label)
    formatter = ticker.ScalarFormatter()
    formatter.set_powerlimits((-2, 2))  #force scientific notation
    sm.yaxis.set_major_formatter(formatter)
#sm.set_ylabel('$\alpha$')
    sm.grid()
    fig.savefig(pp, format=format)
    plt.close(fig)


a = np.squeeze(collect("alpha",path=path))
mask = np.squeeze(collect("alpha_mask",path=path))
n = np.squeeze(collect("n",path=path))

a_smooth = np.squeeze(collect("alpha_smooth",path=path))
a_smooth = a_smooth[2:-2,:]
a_smoothpy = ndimage.gaussian_filter(a, 15)
nx,ny = a.shape

q0= 3.
aa = 45.
b = 55
nu = 2.0

x = b*(.8+.2*np.arange(nx)/nx)
dx = np.squeeze(collect("dx",path=path,xind=[0,0]))
dy = np.squeeze(collect("dz",path=path,xind=[0,0]))

q = q0*np.power(x/aa,2.0)/(1.0- np.power(1.0-x/aa,nu+1.0)*(x<aa))
#pri
#q = q0*np.power(x/aa,2.0)

#pp = PdfPages('Ullmann.pdf')

fig = plt.figure()
a_mean = Frame(a.mean(axis=1),meta={'dx':dx,'sigma':a.std(axis=1),'stationary':True,
                                    'ylabel':r'$\frac{2 \rho_s}{L_{\parallel}}$','linewidth':3,
                                    'xlabel':r'$x/\rho_s$','title':'average ' +r'$\alpha$',
                                    'fontsz':30})
a_mean.render(fig,111)

a_mean.ax.set_xlim(0, dx*nx)
plt.tick_params(axis='both',direction='in',which='both',labelsize=20)
formatter = ticker.ScalarFormatter()
formatter.set_powerlimits((-2, 2))
a_mean.ax.yaxis.set_major_formatter(formatter)
a_mean.ax.xaxis.set_label_coords(.5, -0.07)
plt.tight_layout()
fig.savefig('AverageAlpha.eps')#,pad_inches=.25,bbox_inches='tight')
fig.savefig('AverageAlpha.pdf')
fig = plt.figure()
a_mean = Frame(a.mean(axis=1),
                meta={'dx':dx,'sigma':[np.max(a[:,2:-2],axis=1)-a.mean(axis=1),-np.min((a[:,2:-2]),axis=1)+a.mean(axis=1)],'stationary':True,
                'ylabel':r'$\frac{2 \rho_s}{L_{\parallel}}$','linewidth':3,
                'xlabel':r'$x/\rho_s$','title':'average ' +r'$\alpha$',
                'fontsz':30,'yscale':'linear'})

# a_mean = Frame(a.mean(axis=1),
#                meta={'dx':dx,'sigma':[.0001,.00001],'stationary':True,
#                 'ylabel':r'$\frac{2 \rho_s}{L_{\parallel}}$','linewidth':3,
#                 'xlabel':r'$x/\rho_s$','title':'average ' +r'$\alpha$',
#                 'fontsz':30,'yscale':'linear'})

a_mean.render(fig,111)
plt.tick_params(axis='both',direction='in',which='both',labelsize=20)
formatter = ticker.ScalarFormatter()
formatter.set_powerlimits((-2, 2))
a_mean.ax.yaxis.set_major_formatter(formatter)
a_mean.ax.xaxis.set_label_coords(.5, -0.05)
fig.savefig('AverageAlphaAlt.eps',pad_inches=.25,bbox_inches='tight')
exit()
fig = plt.figure()

# fig_junk = plt.figure()
# img = fig_junk.add_subplot(111)
# cbar = fig.colorbar(img)

a_frm = Frame(a,meta={'stationary':True,'dx':dx,'dy':dy,
                      'xlabel':r'$\rho_s$','fontsz':30})
a_frm.render(fig,111)
cbar = fig.colorbar(a_frm.img,format='%.1g')

# a_frm = Frame(a,meta={'stationary':True,'dx':dx,'dy':dy,
#                                'xlabel':r'$\rho_s$','fontsz':20})

#print dir(fig.colorbar)#.formatter.set_scientific(True)
a_frm.render(fig,111)
#fig.savefig('UllmannContour.eps')


#fig.savefig(pp, format='pdf')
#fig.savefig('UllmannContour.eps')


#fig = plt.figure()
a_frm = Frame(np.log(a+.00026),
              meta={'stationary':True,'dx':dx,'dy':dy,'ylabel':r'$\frac{y}{\rho_s}$',
                    'xlabel':r'$x/\rho_s$','fontsz':30,'contour_only':True,})

#cbar = fig.colorbar(a_frm.img,format='%.1g')   
#print dir(fig.colorbar)#.formatter.set_scientific(True)
a_frm.render(fig,111)

a_frm = Frame(np.log2(a+.000021),
              meta={'stationary':True,'dx':dx,'dy':dy,'ylabel':r'$\frac{y}{\rho_s}$',
                    'xlabel':r'$x/\rho_s$','fontsz':20,'contour_only':True})
a_frm.render(fig,111)

a_frm = Frame(np.log2(a+.000007),
              meta={'stationary':True,'dx':dx,'dy':dy,'ylabel':r'$\frac{y}{\rho_s}$',
                    'xlabel':r'$x/\rho_s$','fontsz':20,'contour_only':True})
a_frm.render(fig,111)
#cbar = fig.colorbar(a_frm.img,format='%.1g')
#fig.savefig(pp, format='pdf')
plt.tick_params(axis='both',direction='in',which='both',labelsize=20)

a_frm.ax.xaxis.set_label_coords(.5, -0.070)

fig.savefig('UllmannMask.pdf',pad_inches=.25,bbox_inches='tight')
fig.savefig('UllmannMask.eps',pad_inches=.25,bbox_inches='tight')


A1 = a_frm.ax.add_patch(patches.Rectangle((50,60),\
70,20,color='g',alpha=.5))
handles, labels = a_frm.ax.get_legend_handles_labels()
handles.append(A1)
labels.append("Area Detail")

leg = a_frm.ax.legend(handles,labels,ncol=2,loc='best',prop={'size':15},fancybox=True)
plt.tick_params(axis='both',direction='in',which='both',labelsize=20)

fig.savefig('UllmannMaskLabel.pdf',pad_inches=.25,bbox_inches='tight')
fig.savefig('UllmannMaskLabel.eps',pad_inches=.25,bbox_inches='tight')
a_detail = a[50/dx:120/dx,60/dy:80/dy]
#levels = np.arange(0, np.max(a_detail),np.max(a_detail)/20.)
a_frm = Frame(a_detail,
              meta={'stationary':True,'dx':dx,'dy':dy,
                    'ylabel':r'$\frac{y}{\rho_s}$','xlabel':r'$x/\rho_s$',
                    'fontsz':20,'contour_only':False,'x0':50,'y0':60})
fig2 = plt.figure()
a_frm.render(fig2,111)

#cbar = fig2.colorbar(a_frm.img,format='%.1g')

a_frm = Frame(np.log2(a_detail),
              meta={'stationary':True,'dx':dx,'dy':dy,
                    'ylabel':r'$\frac{y}{\rho_s}$',
                    'xlabel':r'$x/\rho_s$','fontsz':20,
                    'x0':50,'y0':60,'colors':'g',
                    'contour_only':True,'nlevels':50,
                    'linewidth':.5})
a_frm.render(fig2,111)
a_frm.ax.xaxis.set_label_coords(.5, -0.070)


labels = dict(zip(a_frm.cset.levels,["{:0.1e}".format(np.exp(val)) for val in a_frm.cset.levels]))


plt.clabel(a_frm.cset, a_frm.cset.levels[1::2],  # label every second level
           inline=True,
           fmt=labels,
           fontsize=16)
plt.tick_params(axis='both',direction='in',which='both',labelsize=20)
fig2.savefig('UllmannDetail.eps',pad_inches=.25,bbox_inches='tight')

#pp.close()



#, ticks=[-1, 0, 1])
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])

# pp.close()

# A



# pp = PdfPages('sm.pdf')

# fast2Dplot(pp,a,title='chaotic ' +r"$\alpha$"+' in BOUT++/C++')

# fast2Dplot(pp,np.log(a),title='chaotic ' +r"$\alpha$"+' in BOUT++/C++')
# aveplot(pp,a)


# fast2Dplot(pp,np.log(a_smooth),title='chaotic ' +r"$\alpha$"+'smoothed in BOUt++/C++')
# aveplot(pp,a_smooth)

# fast2Dplot(pp,np.log(a_smoothpy),title='chaotic ' +r"$\alpha$"+' smoothed outside of BOUT++/C++ - no MPI to deal with')
# aveplot(pp,a_smoothpy)

# fig, sm = plt.subplots(1)
# sm.plot(x,q)
# fig.savefig(pp, format='pdf')
# pp.close()
