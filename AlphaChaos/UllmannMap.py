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

from scipy.optimize import curve_fit, root,newton_krylov
from scipy.signal import argrelextrema  
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
    


def fast2Dplot(pp,data,title=None,xlabel=None,ylabel=None,addcurve=None,extent=[0,1,0,1]):
    
    fig = plt.figure()
    
    sm = fig.add_subplot(1,1,1)
    
    im = sm.imshow(np.flipud(np.rot90(data)),aspect='auto',interpolation='none',origin='lower')
    im.set_extent(extent)
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
    
def StandardMap(x,y,L,k,q0,b=30.0):
    
 

    #print 'x: ',x,x+k*np.sin(y),y,np.sin(y)
    
    #one can argue that x = 2*M_PI*fmod(q(R),1)
    #so given x and some q0 s.t. fmod(q0,1) = 0;  q(x) = q0 + x/(2*np.pi)

    #Ullmann Map 
    # x_new = (stopevolve == 0)*(x + k*np.sin(y))+\
    #      (stopevolve)*x
    # y_new = (stopevolve == 0)*(np.mod(y+x_new,2*np.pi))+\
    #     (stopevolve)*y
   
    print 'b: ', b
    aa = -.00
    B_0 = 1.0 #in tesla
    #b = 30 #minor rad
    R_0 = 90 #major rad
    m = 3 #external mode
    l = 10 #coil width
    a= 40
    
    hit_divert = (x>b)
    inCORE = x<b
    stopevolve = hit_divert

    x_new = (stopevolve == 0)*x/(1-aa*np.sin(y))
    q = q0*(x/a)**2
    y_new =  (stopevolve == 0)*(y+ 2*np.pi/q + aa*np.cos(y))
    y_new = np.mod(y_new,2*np.pi)

    #
    

    #see  "DIFFUSIVE TRANSPORT THROUGH A NONTWIST BARRIER IN TOKAMAKS"
    eps = .2
    print m,l,a,R_0,q0,b
    C = ((2*m*l*a**2)/(R_0*q0*b**2))*eps
    print 'C: ', C/eps
    #eps is the ration between limited and plasma currents
 
    #need to find roots of this thing 
    def func(x_out):
        return (-x_new + x_out +(m*b*C)/(m-1)*(x_out/b)**(m-1) *np.sin(m*y_new))**2
    
    #print root(func,x_new)
    x_new2 =  (stopevolve == 0)* (newton_krylov(func,x_new)) + (stopevolve)*x

   # print (-x_new + x_new2 +(m*b*C)/(m-1)*(x_new2/b)**(m-1) *np.sin(m*y_new))**2

    # x_new2 = (stopevolve == 0)*(x_new +(m*b*C)/(m-1)*(x_new/b)**(m-1) *np.sin(m*y_new)) + (stopevolve)*x
    
    
    y_new2 = (stopevolve == 0)*(y_new - C*(x_new2/b)**(m-2) * np.cos(m*y_new))+ (stopevolve)*y
                                
    #print 'xchange:', x_new2/x

    x_new = x_new2
    y_new = np.mod(y_new2,2*np.pi)

    new_inSOL = (x_new >b) 
    #full_orbit = (new_inSOL & inSOL) == False #can visit SOL,but can't stay
    full_orbit = (new_inSOL & inCORE) == False 
    half_orbit = (new_inSOL & inCORE) == True #ok, you hit the divertor

    #print new_inSOL + inSOL
    #q = q0 + x/(2*np.pi)
   # L = L + (stopevolve ==0)*((full_orbit)*q *100* 2*np.pi + \
   #                            half_orbit *100* 2*q* np.pi)
    L = L + (full_orbit)# + .5*half_orbit
    #make sure that fieldlines that hit the divertor with 0<x<pi don't get remapped
    #print 'x: ',half_orbit
    #x_new = x_new + half_orbit*np.pi
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
    

 
def setup_xz(nx=128,nz=128,b = .3,edge=None,rmin=0.0,rmax=1.0):
    # if edge==None:
    #     #dx = (2*np.pi)/(nx)
    #     dx = np.float(b*(rmax-rmin))/nx
    #     print b,nx,dx
    #     #x = np.arange(-1*np.pi,np.pi,dx)
    #     x = np.arange(rmin*b,rmax*b,dx)
    #     print x.shape,nx,nz
    #     x = np.repeat(x,nz)
    #     x = x.reshape(nx,nz)
    # else:
    #     #x_b = edge #in rads, nz long
    #     x =[]   
    #     #edge  = -edge
    #     #edge[1:10] = np.pi/2
    #     for i,xmax in enumerate(edge):
    #         #xmax = -2.5
    #         #xmax = np.min(edge)
    #         x.append(np.linspace(xmax-np.pi/10.0,xmax+np.pi/20.0,nx))
    #         #x.append(np.linspace(-np.pi,xmax-.4,nx))
    #     x = np.array(x)
    #     print x.shape
    #     x = np.transpose(x)
        
    # dz = (2*np.pi)/(nz)
    # z = np.arange(0,2*np.pi,dz)
    # z = np.repeat(z,nx)
    # z = np.transpose(z.reshape(nz,nx))

    
    x, z = np.mgrid[b*rmin:b*rmax:complex(0,nx),0:2*np.pi:complex(0,nz)]

    #print x.shape
    return x,z

def StandardLength(x,z,k=1,max_pol_orbits=100,q=5.0,b=.3):
    
    print  max_pol_orbits
    print 'xhape ',  x.shape
    keep_i= list(np.where(x < b))
    
    
    
    L = 0.0*z
    count = 0

    while count<max_pol_orbits and len(x[keep_i]) !=0:
        #print count
        x[keep_i], z[keep_i],L[keep_i] = StandardMap(x[keep_i],z[keep_i],L[keep_i],k,q,b=b)
        print count,' : ',100.0*len(x[keep_i])/np.size(x),'% of the field-lines jumping'
        
    
        keep_i= list(np.where((x < b)))
        count+=1

    return L


def saveAlphaMap(ncells =32,k=1.5,q=5):
    x,z = setup_xz(nx=ncells,nz=ncells)
    
    print 'x.shape: ',x.shape

    L = StandardLength(x,z)

    try:
        from boututils import write_grid
    except:
        'no write_grid in $BOUT_TOP/tools/pylib/boututils'
    #we will need to create a 3D field with y extents =1 for BOUT
    
    write_grid(gridfile='alpha_map.nc',nx=ncells+4,ny=1,dx=1,dy=1)
    #write_grid(

def showLmap(ncells=32,k=1.5,q=5,b=45,rmin =0.0,rmax = 1.0):

    #from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('sm.pdf')

    x,z = setup_xz(nx=ncells,nz=ncells,b=b,rmin=rmin,rmax=rmax)
    #print x
    x_b,z_b = edge_finder(ncells,ncells,k)
    L = StandardLength(x,z,k=1,max_pol_orbits = 100,b=b)
    Lsmooth = ndimage.gaussian_filter(L, 3)
    #fast2Dplot(pp,np.log(L),extent=[rmin,rmax,0,2*np.pi])
    fast2Dplot(pp,L,extent=[rmin,rmax,0,2*np.pi])
    x_b,z_b = to_index_coord(x_b,z_b,ncells,ncells)
    
    # x,z = setup_xz(nx=ncells,nz=ncells)
    # L = StandardLength(x,z,k=k,max_pol_orbits = 20)
    #a = .2/L   
    a = L
 
    fig, sm = plt.subplots(1)
    
    alpha = a.mean(axis=1)
    sigma = a.std(axis=1)
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
    #sm.set_yscale('log')

    # plt.contourf(np.flipud(np.rot90(np.log(L))),80)
    fig.savefig(pp, format='pdf')
    plt.close(fig)

    
    pp.close()
    
    
showLmap(ncells=100,q=3.0,rmin = .92,rmax = 1.0,b=45)
#saveAlphaMap()
