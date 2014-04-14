import os,sys
boutpath = os.environ['BOUT_TOP']
pylibpath = boutpath+'/tools/pylib'
pbpath = pylibpath+'/post_bout'
boutdatapath = pylibpath+'/boutdata'
boututilpath = pylibpath+'/boututils'

allpath = [boutpath,pylibpath,pbpath,boutdatapath,boututilpath]

[sys.path.append(elem) for elem in allpath]
print(sys.path)

from scipy.optimize import newton_krylov
from scipy.signal import argrelextrema  
#import gobject
import numpy as np
print ('in post_bout/post_bout.py')
#from ordereddict import OrderedDict
#from scipy.interpolate import interp2d,interp1d
from scipy import ndimage
from copy import copy


from read_cxx import read_cxx, findlowpass
from boutdata import collect
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.artist as artist 
import matplotlib.ticker as ticker


    

from multiprocessing import Process, Manager
import platform,multiprocessing

import os
import time


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


def go_forward(x,y,a=40,b=50,R_0 = 90,l=10,m=3,aa=0.0,q0=3.0):

    hit_divert = (x>b)
    inCORE = x<b
    stopevolve = hit_divert
    
    eps = .2
    C = ((2*m*l*a**2)/(R_0*q0*b**2))*eps

    x_new = x/(1-aa*np.sin(y))
    q = q0*(x_new/a)**2
    y_new =  (y+ 2*np.pi/q + aa*np.cos(y))
    y_new = np.mod(y_new,2*np.pi)

    def func(x_out):
        return (-x_new + x_out +(m*b*C)/(m-1)*(x_out/b)**(m-1) *np.sin(m*y_new))**2
    
    x_new2 = (newton_krylov(func,x_new,method='gmres',maxiter=50))
    y_new2 = (y_new - C*(x_new2/b)**(m-2) * np.cos(m*y_new))
                                
    #print 'xchange:', x_new2/x

    x_new = x_new2
    y_new = np.mod(y_new2,2*np.pi)

    
    return x_new,y_new

def go_back(x,y,a=40,b=50,R_0 = 90,l=10,m=3,aa=0.0,q0=3.0,eps=.07):

    hit_divert = (x>b)
    inCORE = x<b
    stopevolve = hit_divert
    
    #eps = .2
    C = ((2*m*l*a**2)/(R_0*q0*b**2))*eps
    
    def func(y_out):
        return (-y + y_out - C*(x/b)**(m-2) *np.cos(m*y_out))**2
    
    def func2(y_out):
        return (-y_old + y_out + (2.0*np.pi/q) + aa*np.cos(y_out))**2
    
    y_old = copy(y)
    y_old = (newton_krylov(func,y))
    y_old = np.mod(y_old,2.0*np.pi)

    x_old = x + (m*b*C)/(m-1)*(x/b)**(m-1) *np.sin(m*y_old)

    q = q0*(x_old/a)**2

    y_old2 = copy(y_old)
    y_old2 = (newton_krylov(func2,y_old))

    #y_old2 = y_old - 2*np.pi/q #- aa*np.cos(
    y_old2 = np.mod(y_old2,2.0*np.pi)
    x_old2 = x_old*(1.0 -aa*np.sin(y_old2))



    return x_old2,y_old2


def to_index_coord(x,y,nx,ny):
    x_i = x*(nx/(2*np.pi))+(nx/2.)
    y_i = y*(ny/(2*np.pi))  
    
    return x_i,y_i
    
def StandardMap(x,y,L,k,q0,b=30.0,aa=0.0,eps=.3,m =3):

    print('b: ', b)
    #aa = -.00
    B_0 = 1.0 #in tesla
    #b = 30 #minor rad
    R_0 = 90 #major rad
    #m = 3. #external mode
    l = 10 #coil width
    a= 40
    beta = 2.0
    mu = 1.0
    beta_p = beta*(mu+1)/(beta+mu+1)
    
 
    
    hit_divert = (x>b)
    inCORE = x<b
    stopevolve = hit_divert
    
    x_new = (stopevolve == 0)*x/(1-aa*np.sin(y))
    
    # q = (q0*(x_new/a)**2) #*(1 -(1 + beta_p * (x_new/a)**2)*
    #  (1.0- (x_new/a)**2)*(a>x_new))**-1
    
    xx = x_new/a
    nu = 2.

    q = q0*(xx**2)/(1-(1-xx)**(nu+1))
   
    print('q: ', q.min(),q.max())
    y_new =  (stopevolve == 0)*(y+ 2*np.pi/q + aa*np.cos(y))
    y_new = np.mod(y_new,2*np.pi)

    #
    

    #see  "DIFFUSIVE TRANSPORT THROUGH A NONTWIST BARRIER IN TOKAMAKS"
    #eps = .3
    print(m,l,a,R_0,q0,b,eps)
    C = ((2*m*l*a**2)/(1.0*R_0*q0*b**2))*eps
    print('C: ', C/eps)
    #eps is the ration between limited and plasma currents
 
    #need to find roots of this thin2g 
    def func(x_out):
        return (-x_new + x_out +((m*b*C)/(m-1))*((x_out/b)**(m-1) )*np.sin(m*y_new))**2
    
    x_new2 = copy(x_new)

    #x_new2 =  (stopevolve == 0)* (newton_krylov(func,x_new2)) + (stopevolve)*x

    x_new2 =  (stopevolve == 0)* (newton_krylov(func,x_new2,method='bicgstab')) + (stopevolve)*x
    #print (-x_new + x_new2 +(m*b*C)/(m-1)*((x_new2/b)**(m-1) )*np.sin(m*y_new))**2
    
    y_new2 = (stopevolve == 0)*(y_new - C*(x_new2/b)**(m-2) * np.cos(m*y_new))+ (stopevolve)*y
                                
    #print 'xchange:', x_new2/x
    #print(x_new2,x)
    x_new = x_new2
    y_new = np.mod(y_new2,2*np.pi)

    stay_inCORE = (x_new <b) 
    #full_orbit = (new_inSOL & inSOL) == False #can visit SOL,but can't stay
    full_orbit = (stay_inCORE & inCORE) == True 
    half_orbit = (stay_inCORE  & inCORE) == True #ok, you hit the divertor

    #print new_inSOL + inSOL
    #q = q0 + x/(2*np.pi)
    L = L + (stopevolve ==0)*((full_orbit)*q *R_0* 2*np.pi)# + \
                               #half_orbit *100* 2*q* np.pi)
    #L = L +(full_orbit)# + .5*half_orbit
 

    return x_new,y_new,L


    
def showXhist(a=40,b=50,R_0 = 90,l=10,m=3,aa=0.0):
    
    rmin = 1.0
    rmax = 1.2

    ncells  = 8
    x0,y0 = setup_xz(nx=ncells,nz=ncells,b=b,rmin=rmin,rmax=rmax)

    #x0 = b/2.0
    #y0 = np.pi
    xhist = []
    yhist= []
    x = x0
    y = y0
  
    
    for i,xmax in enumerate(xrange(10)):
        xhist.append(x)
        yhist.append(y)
        x, y = go_forward(x,y,a=a,b=b,R_0 = R_0,l=l,m=3,aa=aa,eps=.07)
        print(i,x.shape)

    pp = PdfPages('xhist.pdf')  
 
    xhist = np.array(xhist)
    yhist = np.array(yhist)

    
    fig, sm = plt.subplots(1)
    sm.plot(xhist.flatten(),yhist.flatten(), lw=2,linestyle='None',marker='.')
    #sm.plot(xhist,yhist, lw=2)
    sm.grid()
    fig.savefig(pp, format='pdf')

    fig, sm = plt.subplots(1)
    for i,ic in enumerate(xrange(ncells)):
        sm.plot(xhist[:,i,:].flatten(),yhist[:,i,:].flatten(), lw=2,linestyle='None',marker='.',markersize=3,alpha = .5)
    sm.grid()
    fig.savefig(pp, format='pdf')

    # fig, sm = plt.subplots(1)
    # sm.plot(xhist,lw=2,linestyle='None',marker='.')
    # #sm.plot(xhist,yhist, lw=2)
    # sm.grid()
    # fig.savefig(pp, format='pdf')

    plt.close(fig)  
    pp.close()


def showXrev(a=40,b=50,R_0 = 90,l=10,m=3,aa=0.0,throw_away = True,
             name='Ullmann2',cached=False,compare=False):
    ncells = 50
    xcells = 20
    x0,y0 = np.mgrid[b:1.01*b:complex(0,xcells),0:2*np.pi:complex(0,ncells)]
  
    xhist = []
    yhist= []
    x = x0
    y = y0


 

    if cached:
        Hist = (np.load('lastX.npy')).item()
        xhist = Hist['x']
        yhist = Hist['y']
    else:
        for i,xmax in enumerate(xrange(1000)):
            xhist.append(x)
            yhist.append(y)
            x, y = go_back(x,y,a=a,b=b,R_0 = R_0,l=l,m=m,aa=aa,eps=.3) # jump back
            keep_i= list(np.where(x < b)) #see which ones return to the CORE
            
            print(i ," " ,len(keep_i),100.0*len(x[keep_i])/np.size(x))

    
        xhist = np.squeeze(np.array(xhist))
        yhist = np.squeeze(np.array(yhist)) 
        Hist={'x':xhist,'y':yhist}
        np.save('lastX',Hist)


    pp = PdfPages(name+'.pdf')  

    print(xhist.shape)

    fig, sm = plt.subplots(1)
    sm.plot(xhist.flatten(),yhist.flatten(), lw=2,linestyle='None',marker='.',rasterized=True)
    #sm.plot(xhist,yhist, lw=2)
    sm.grid()
    fig.savefig(pp, format='pdf')

    fig, sm = plt.subplots(1)
    sm.plot(xhist[0:20,12,5]/b, lw=2,linestyle='None',marker='.',rasterized=True)
    #sm.plot(xhist,yhist, lw=2)
    sm.grid()
    fig.savefig(pp, format='pdf')

    #for all x<b find the last point oging back that where x>b and get the connection length
    #for x in xhist[:,15,6]:
        


    fig, sm = plt.subplots(1)
 
    #cached = True
    if compare:
        Ldict = (np.load('lastL.npy')).item()
        im = sm.imshow((np.flipud(np.rot90(Ldict['data']))),aspect='auto',interpolation='none',origin='lower')
        im.set_extent([Ldict['rmin'],Ldict['rmax'],0,2*np.pi])
        #fast2Dplot(pp,L,extent=[rmin,rmax,0,2*np.pi])

    for i,ic in enumerate(xrange(xcells)):
        sm.plot(xhist[:,i,:].flatten()/b,yhist[:,i,:].flatten(), lw=2,linestyle='None',
                marker='.',markersize=1,alpha = .05,rasterized=True)
    sm.set_ylim([0,2*np.pi])
    sm.set_xlim([.8,1])
    sm.grid()
    fig.savefig(pp, format='pdf')

    
    plt.close(fig)  
    pp.close()


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

def StandardLength(x,z,procnum,return_dict,k=1,max_pol_orbits=100,q=5.0,b=.3,
                   aa=0.0,eps = .3,m=3):
    
    print (max_pol_orbits)
    print ('xhape ',  x.shape)
    print('max orbits: ', max_pol_orbits)
    keep_i= list(np.where(x < b))
    
    
    
    L = 0.0*z
    count = 0

    while count<max_pol_orbits and len(x[keep_i]) !=0:
        #print count, x and y get replaced
        x[keep_i], z[keep_i],L[keep_i] = StandardMap(x[keep_i],z[keep_i],L[keep_i],k,q,b=b,aa=aa,eps=eps,m = m)
        # if lyapunov:
        #     x[keep_i], z[keep_i],L[keep_i] = StandardMap(x[keep_i],z[keep_i],L[keep_i],k,q,b=b,aa=aa,eps=eps)  
        
        
        print (count,' : ',100.0*len(x[keep_i])/np.size(x),'% of the field-lines jumping')
        
    
        keep_i= list(np.where((x < b)))
        count+=1

    R_0 = 90
    L[keep_i] = count*10*R_0*2*np.pi
    #print(dir(return_dict))
    
    return_dict[procnum] = L

def simple(i,return_dict):
    return_dict[i] = 5

def saveAlphaMap(ncells =32,k=1.5,q=5):
    x,z = setup_xz(nx=ncells,nz=ncells)
    
    print ('x.shape: ',x.shape)

    L = StandardLength(x,z)

    try:
        from boututils import write_grid
    except:
        'no write_grid in $BOUT_TOP/tools/pylib/boututils'
    #we will need to create a 3D field with y extents =1 for BOUT
    
    write_grid(gridfile='alpha_map.nc',nx=ncells+4,ny=1,dx=1,dy=1)
    #write_grid(

def showLmap(ncells=32,k=1.5,q=5,b=45,rmin =0.0,rmax = 1.0,
             aa=0.0,max_orbits = 100,cached = True,eps = .3,
             m=3):

    #from matplotlib.backends.backend_pdf import PdfPages
    pp = PdfPages('sm.pdf')

    #x,z = setup_xz(nx=ncells,nz=ncells,b=b,rmin=rmin,rmax=rmax)
    x,z = np.mgrid[b*rmin:b*rmax:complex(0,ncells),0:2*np.pi:complex(0,ncells)]
    
    dx = 1e-11
    dz = 1e-11
    
    x1,z1 = x+dx,z+dz
    
    
  
    #print x
    #x_b,z_b = edge_finder(ncells,ncells,k)

    #f = open('lastL', 'w')

    if cached:
        L = np.load('lastL.npy')
        L = L.item()
        rmin = L['rmin']
        rmax = L['rmax']
        L = L['data']
        
    else:
        L = StandardLength(x,z,q=q,k=k,max_pol_orbits = max_orbits,
                           b=b,aa=aa,eps = eps,m =m)
        Ldict={'data':L,'rmin':rmin,'rmax':rmax}
        np.save('lastL',Ldict)
    
    #print 'L.shape: '  
    #f.close()

    fast2Dplot(pp,np.log(L),extent=[rmin,rmax,0,2*np.pi])
    fast2Dplot(pp,L,extent=[rmin,rmax,0,2*np.pi])
   
    #fast2Dplot(pp,np.log(L),extent=[rmin,rmax,0,2*np.pi])
    #x_b,z_b = to_index_coord(x_b,z_b,ncells,ncells)
    
    # x,z = setup_xz(nx=ncells,nz=ncells)
    # L = StandardLength(x,z,k=k,max_pol_orbits = 20)
    a = .2/L   
    #a = L
    fast2Dplot(pp,a,extent=[rmin,rmax,0,2*np.pi])
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
    
    
# showLmap(ncells=20,q=3.0,rmin = .82,rmax = 1,b=50,aa=-0.00,
#          cached=False,max_orbits = 100,eps = .2)


#print 'Main Process start'
starttime = time.time()
processlist = []
Ncpu = multiprocessing.cpu_count()
# p1 = Process(target=showLmap(ncells=300,q=3,rmin = .90,
#                              rmax = 1.05,b=50,aa=-0.02,
#                              cached=False,max_orbits = 10,
#                              eps = 0.07,m=7))
chunk = 50 
nx = Ncpu*chunk
nz = 200

b = 50
rmax = .98
rmin = .92

manager = Manager()
return_L = manager.dict()

x,z = np.mgrid[b*rmin:b*rmax:complex(0,nx),0:2*np.pi:complex(0,nz)]

#break up points amongs cpus

# return_L[0] = 5
for i in range(Ncpu):
    print(i)
    #p = Process(target=simple, args =(i,return_L))
    p = Process(target=StandardLength, 
                args =(x[chunk*i:chunk*(i+1),:],z[chunk*i:chunk*(i+1),:],i,return_L),
                kwargs = {'q':3,'k':1.5,'max_pol_orbits':100,
                          'b':b,'eps':.07,'m':7,'aa':-.02})
    processlist.append(p)
    p.start()
    
for p in processlist:
    p.join()

L = np.concatenate(return_L)

print(return_L[3].shape)
#print(x[chunk*1:(chunk)*2,:].shape)


pp = PdfPages('sm.pdf')
fast2Dplot(pp,np.log(L),extent=[rmin,rmax,0,2*np.pi])
pp.close()
#fast2Dplot(pp,L,extent=[rmin,rmax,0,2*np.pi])
#p1 = Process(target=StandardLength, args =(x,z,0,return_L))
#print return_L.values()
             # kwargs=(q=3,k=1.5,
             #         max_pol_orbits = 10,
             #         b=50,aa=-.01,eps = .07,m =7))





#saveAlphaMap()
#showXrev(aa=-.04)
#showXrev(aa=-.001,compare=True,cached=True,eps=1.0)
#showLmap(ncells=20,q=3.0,rmin = .5,rmax = .9,b=50,aa=-.05)
