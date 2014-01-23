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
from boutdata import collect2 as collect
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
#from boutdata import collect
from boututils import showdata
import numpy as np
from scipy import ndimage
from frame import Frame

try:
    path = sys.argv[1]
except:
    path = '/tmp/SOLblob/data_blob_jpar'
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
#n = np.squeeze(collect("n",path=path))

a_smooth = np.squeeze(collect("alpha_smooth",path=path))
print a_smooth.shape
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

pp = PdfPages('Ullmann.pdf')
fig = plt.figure()

# fig_junk = plt.figure()
# img = fig_junk.add_subplot(111)
# cbar = fig.colorbar(img)

a_frm = Frame(a,meta={'stationary':True,'dx':dx,'dy':dy,
                      'xlabel':r'$\rho_s$','fontsz':20})


#print dir(fig.colorbar)#.formatter.set_scientific(True)
a_frm.render(fig,111)

cbar = fig.colorbar(a_frm.img,format='%.1g')
fig.savefig(pp, format='pdf')

a_frm = Frame(mask,meta={'stationary':True,'dx':dx,'dy':dy,
                      'xlabel':r'$\rho_s$','fontsz':20})


#print dir(fig.colorbar)#.formatter.set_scientific(True)
a_frm.render(fig,111)

#cbar = fig.colorbar(a_frm.img,format='%.1g')
fig.savefig(pp, format='pdf')


# a_frm = Frame(n[-1,:,:],meta={'stationary':True,'dx':dx,'dy':dy,
#                       'xlabel':r'$\rho_s$','fontsz':20})


#print dir(fig.colorbar)#.formatter.set_scientific(True)
a_frm.render(fig,111)

cbar = fig.colorbar(a_frm.img,format='%.1g')
fig.savefig(pp, format='pdf')



a_frm = Frame(np.log(a),meta={'stationary':True,'dx':dx,'dy':dy,
                              'xlabel':r'$\rho_s$','fontsz':20,
                              'contour_only':True,'alpha':.4,
                              'title':r'$Ullmann$ '+ r'$\alpha$'})


#print dir(fig.colorbar)#.formatter.set_scientific(True)
a_frm.render(fig,111,rasterized=False)

#cbar = fig.colorbar(a_frm.img,format='%.1g')
#cbar.formatter.set_scientific(True)
#cbar.formatter.set_powerlimits((-2, 2))
#a_frm.render(fig,111)
#cbar = fig.colorbar(a_frm.img,ax = a_frm.ax)


fig.savefig(pp, format='pdf')
aveplot(pp,a,dx =dx)
fig = plt.figure()
a_mean = Frame(a.mean(axis=1),meta={'dx':dx,'sigma':a.std(axis=1),'stationary':True,
                                    'xlabel':r'$\rho_s$','title':'average ' +r'$\alpha$',
                                    'fontsz':20})
a_mean.render(fig,111)
formatter = ticker.ScalarFormatter()
formatter.set_powerlimits((-2, 2))
a_mean.ax.yaxis.set_major_formatter(formatter)
a_mean.ax.xaxis.set_label_coords(.5, -0.05)
fig.savefig(pp,format='pdf')
#, ticks=[-1, 0, 1])
#cbar.ax.set_yticklabels(['< -1', '0', '> 1'])

pp.close()





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
