import os,sys
boutpath = os.environ['BOUT_TOP']
pylibpath = boutpath+'/tools/pylib'
pbpath = pylibpath+'/post_bout'
boutdatapath = pylibpath+'/boutdata'
boututilpath = pylibpath+'/boututils'
allpath = [boutpath,pylibpath,pbpath,boutdatapath,boututilpath]
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/boutdata')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/boututils')
    # sys.path.append('/home/cryosphere/BOUT/tools/pylib/post_bout')
    #sys.path.append(allpath)
[sys.path.append(elem) for elem in allpath]
print sys.path
from ordereddict import OrderedDict
from scipy.interpolate import interp2d,interp1d
from boutdata import collect
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker
from boutdata import collect
from boututils import showdata
import numpy as np
from scipy import ndimage

try:
    path = sys.argv[1]
except:
    path = '/tmp/SOLblob/data_blob'
def fast2Dplot(pp,data,title=None,xlabel=None,ylabel=None,addcurve=None):
    
    fig = plt.figure()
    
    sm = fig.add_subplot(1,1,1)
    
    sm.imshow(np.flipud(np.rot90(data)),aspect='auto',interpolation='none',origin='lower')
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


def aveplot(pp,data,title='chaotic '+r"$\alpha \pm \sigma$",x_label='x index',
            label=r"$\alpha$",format='pdf'):
    fig, sm = plt.subplots(1)
    
    alpha = data.mean(axis=1)
    sigma = data.std(axis=1)

    r = np.arange(data.shape[0])
#line, = sm.plot(alpha, color='blue', lw=2)
    sm.plot(r,alpha, lw=2, label=label, color='blue')
    sm.fill_between(r,alpha+sigma, alpha-sigma, facecolor='yellow', alpha=0.5)
    sm.set_title(title)
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
edge = np.squeeze(collect("edge",path=path))
a_smooth = np.squeeze(collect("alpha_smooth",path=path))
a_smooth = a_smooth[2:-2,:]
a_smoothpy = ndimage.gaussian_filter(a, 15)

pp = PdfPages('sm.pdf')

fast2Dplot(pp,np.log(a),title='chaotic ' +r"$\alpha$"+' in BOUT++/C++')
aveplot(pp,a)

fast2Dplot(pp,np.log(a_smooth),title='chaotic ' +r"$\alpha$"+'smoothed in BOUt++/C++')
aveplot(pp,a_smooth)

fast2Dplot(pp,np.log(a_smoothpy),title='chaotic ' +r"$\alpha$"+' smoothed outside of BOUT++/C++ - no MPI to deal with')
aveplot(pp,a_smoothpy)

pp.close()