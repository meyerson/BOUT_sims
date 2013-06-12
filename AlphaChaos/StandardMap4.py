
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


#import gobject
import numpy as np
print 'in post_bout/post_bout.py'
#from ordereddict import OrderedDict
#from scipy.interpolate import interp2d,interp1d
from scipy import ndimage

from read_cxx import read_cxx, findlowpass
from boutdata import collect
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker

import subprocess 
    


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


def go_forward(x,y,k):
    x_new = (x + k*np.sin(y))
    y_new = np.mod(y+x_new,2*np.pi)
    return x_new,y_new

def go_back(x,y,k):
    y_old = np.mod(y-x,2*np.pi)
    x_old = x-k*np.sin(y_old)
    return x_old,y_old

def to_index_coord(x,y,nx,ny):
    x_i = x*(nx/(2*np.pi))+(nx/2.)
    y_i = y*(ny/(2*np.pi))  
    
    return x_i,y_i
    
def StandardMap(x,y,L,k,q0):
    
    hit_divert = (x>=np.pi)
    inSOL = x>0
    #toolongL = L>1000.0
    #stopevolve = (toolongL + hit_divert)>0
    stopevolve = hit_divert

    #print 'x: ',x,x+k*np.sin(y),y,np.sin(y)
    
    #one can argue that x = 2*M_PI*fmod(q(R),1)
    #so given x and some q0 s.t. fmod(q0,1) = 0;  q(x) = q0 + x/(2*np.pi)

    #Chirikov-Taylor
    x_new = (stopevolve == 0)*(x + k*np.sin(y))+\
         (stopevolve)*x
    y_new = (stopevolve == 0)*(np.mod(y+x_new,2*np.pi))+\
        (stopevolve)*y
        
    new_inSOL = (x_new >0) 
    #full_orbit = (new_inSOL & inSOL) == False #can visit SOL,but can't stay
    full_orbit = (new_inSOL & inSOL) == False 
    half_orbit = (new_inSOL & inSOL) == True #ok, you hit the divertor

    #print new_inSOL + inSOL
    q = q0 + x/(2*np.pi)
    L = L + (stopevolve ==0)*((full_orbit)*q *100* 2*np.pi + \
                               half_orbit *100* 2*q* np.pi)

    #make sure that fieldlines that hit the divertor with 0<x<pi don't get remapped
    #print 'x: ',half_orbit
    x_new = x_new + half_orbit*np.pi
    #will satisfy hit_divert at new iteration, so stopevole = True at next call
   # print 'x: ',x_new

    return x_new,y_new,L


    

#let's keep this vectorizedd

def edge_finder(nx,nz,k=1.0):
    z_b = np.arange(0,nz,1)
    z_b = 2*np.pi*(1.0*z_b/nz)
    x_b = k*np.sin(z_b)
    #print x_b
    return x_b,z_b

# def edge_detect(data,axis =1,sampe):
#     np.gradient(data)
    

 
def setup_xz(nx=128,nz=128,edge=None):
    if edge==None:
        dx = (2*np.pi)/(nx)
        x = np.arange(-1*np.pi,np.pi,dx)
        #print x
        x = np.repeat(x,nz)
        x = x.reshape(nx,nz)
    else:
        #x_b = edge #in rads, nz long
        x =[]   
        #edge  = -edge
        #edge[1:10] = np.pi/2
        for i,xmax in enumerate(edge):
            #xmax = -2.5
            #xmax = np.min(edge)
            x.append(np.linspace(xmax-np.pi/10.0,xmax+np.pi/20.0,nx))
            #x.append(np.linspace(-np.pi,xmax-.4,nx))
        x = np.array(x)
        print x.shape
        x = np.transpose(x)
        
    dz = (2*np.pi)/(nz)
    z = np.arange(0,2*np.pi,dz)
    z = np.repeat(z,nx)
    z = np.transpose(z.reshape(nz,nx))


    return x,z

def StandardLength(x,z,k=1,max_pol_orbits=200,q=5.0):
    
    print  max_pol_orbits
    #print x.shape
    keep_i= list(np.where(x < np.pi))
    
    
    
    L = 0.0*z
    count = 0

    while count<max_pol_orbits and len(x[keep_i]) !=0:
        x[keep_i], z[keep_i],L[keep_i] = StandardMap(x[keep_i],z[keep_i],L[keep_i],k,q)
        print count,' : ',100.0*len(x[keep_i])/np.size(x),'% of the field-lines jumping'
        
    
        keep_i= list(np.where((x < np.pi)))
        count+=1

    return L


def saveAlphaMap(ncells =32,k=1.5,q=5):
    x,z = setup_xz(nx=ncells,nz=ncells)
    
    L = StandardLength(x,z)

    try:
        from boututils import write_grid
    except:
        'no write_grid in $BOUT_TOP/tools/pylib/boututils'
    #we will need to create a 3D field with y extents =1 for BOUT
    
    write_grid(gridfile='alpha_map.nc',nx=ncells+4,ny=1,dx=1,dy=1)
    #write_grid(

def showLmap(ncells=32,k=1.5,q=5):

    #from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('sm.pdf')

    x,z = setup_xz(nx=ncells,nz=ncells)
    x_b,z_b = edge_finder(ncells,ncells,k)
    L = StandardLength(x,z,k=1,max_pol_orbits = 40)
    Lsmooth = ndimage.gaussian_filter(L, 3)
    fast2Dplot(pp,Lsmooth)

    x_b,z_b = to_index_coord(x_b,z_b,ncells,ncells)
    
    # x,z = setup_xz(nx=ncells,nz=ncells)
    # L = StandardLength(x,z,k=k,max_pol_orbits = 20)
    a = .2/L   
   
    
    #sloppy . .
    da = (np.gradient(Lsmooth,ncells/10.0))
    detect_edge = da[0]**2 + da[1]**2
    x_b = np.argmax(detect_edge,axis=0)
    
    # x_b_fft = np.fft.rfft(x_b)
    # x_b_fft[ncells/2::]=0.0
    # x_b = np.real(np.fft.irfft(x_b_fft))
    x_b = ndimage.gaussian_filter(x_b,10)

    
    edge = (x_b - ncells/2.0)*(2*np.pi/ncells)
    print edge
    altx,altz = setup_xz(nx=ncells,nz=ncells,edge =edge)
    
    altL = StandardLength(altx,altz,k=k,max_pol_orbits=1000)
    a_new = .2/altL
    #fast2Dplot(pp,np.log(detect_edge),title='find the edge')
    fast2Dplot(pp,detect_edge,title='find the edge',
               addcurve={'x':x_b,'y':z_b,'color':'b'})
   
    fast2Dplot(pp,np.log(Lsmooth),title='original chaotic Log('+r"$\alpha$"+")",
               addcurve={'x':x_b,'y':z_b,'color':'r'})

    
    fast2Dplot(pp,np.log(a_new),title='new chaotic Log('+r"$\alpha$"+")")
    fast2Dplot(pp,np.log(ndimage.gaussian_filter(a_new, 10)),title='new smoothed chaotic Log('+r"$\alpha$"+")")

    
    # fig = plt.figure()
    # sm = fig.add_subplot(1,1,1)

    fig, sm = plt.subplots(1)
    
    alpha = a_new.mean(axis=1)
    sigma = a_new.std(axis=1)
#print alpha
#print sigma/alpha


    r = np.arange(a.shape[0])
#line, = sm.plot(alpha, color='blue', lw=2)
    sm.plot(r,alpha, lw=2, label=r"$\alpha$", color='blue')
    sm.fill_between(r,alpha+sigma, alpha-sigma, facecolor='yellow', alpha=0.5)
    sm.set_title('chaotic '+r"$\alpha \pm \sigma$")
    sm.legend(loc='upper left')
    sm.set_xlabel('x index')
    formatter = ticker.ScalarFormatter()
    formatter.set_powerlimits((-2, 2))  #force scientific notation
    sm.yaxis.set_major_formatter(formatter)
#sm.set_ylabel('$\alpha$')
    sm.grid()

    # plt.contourf(np.flipud(np.rot90(np.log(L))),80)
    fig.savefig(pp, format='pdf')
    plt.close(fig)

    
    pp.close()
    
    
showLmap(ncells=100,k=.95)
#saveAlphaMap()
