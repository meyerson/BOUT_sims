import os,sys
boutpath = os.environ['BOUT_TOP']
pylibpath = boutpath+'/tools/pylib'
pbpath = pylibpath+'/post_bout'
boutdatapath = pylibpath+'/boutdata'
boututilpath = pylibpath+'/boututils'

allpath = [boutpath,pylibpath,pbpath,boutdatapath,boututilpath]

[sys.path.append(elem) for elem in allpath]
print(sys.path)

from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

from scipy.optimize import newton_krylov
from scipy.signal import argrelextrema  
#import gobjectW
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
from sympy import *
dr, dth = symbols("dr dth")
from sympy import log, sin, cos, tan, Wild, Mul, Add
from sympy import Function, Symbol
import sympy as sp
#from mpmath import *
#create some vars with global scope, these will be changed 
#elsewhere in the script with user input

q0,a,b,m,mu = 3.,60.,68.,3.,1.
r0, r1,r2,th0,th1,th2 = symbols("r0 r1 r2 theta0 th1 th2")
dr0, dr1,dr2,dth0,dth1,dth2 = symbols("dr0 dr1 dr2 dth0 dth1 dth2")
eps = symbols("epsilon")
d = symbols("delta")
perturb = [(r0,r0+d*dr0),(th0,th0+d*dth0),
               (r1,r1+d*dr1),(th1,th1+d*dth1),
               (r2,r2+d*dr2),(th2,th2+d*dth2)]
#r1n = r1/a
bprime = 2
l= 10
R_0 = 85.0
aa = -.01
eps = .2
init_printing()

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


# def go_forward(x,y,a=40,b=50,R_0 = 90,l=10,m=3,aa=0.0,q0=3.0):

#     hit_divert = (x>b)
#     inCORE = x<b
#     stopevolve = hit_divert
    
#     eps = .2
#     C = ((2*m*l*a**2)/(R_0*q0*b**2))*eps

#     x_new = x/(1-aa*np.sin(y))
#     q = q0*(x_new/a)**2
#     y_new =  (y+ 2*np.pi/q + aa*np.cos(y))
#     y_new = np.mod(y_new,2*np.pi)

#     def func(x_out):
#         return (-x_new + x_out +(m*b*C)/(m-1)*(x_out/b)**(m-1) *np.sin(m*y_new))**2
    
#     x_new2 = (newton_krylov(func,x_new,method='gmres',maxiter=50))
#     y_new2 = (y_new - C*(x_new2/b)**(m-2) * np.cos(m*y_new))
                                
#     #print 'xchange:', x_new2/x

#     x_new = x_new2
#     y_new = np.mod(y_new2,2*np.pi)

    
#     return x_new,y_new

# def go_back(x,y,a=40,b=50,R_0 = 90,l=10,m=3,aa=0.0,q0=3.0,eps=.07):

#     hit_divert = (x>b)
#     inCORE = x<b
#     stopevolve = hit_divert
    
#     #eps = .2
#     C = ((2*m*l*a**2)/(R_0*q0*b**2))*eps
    
#     def func(y_out):
#         return (-y + y_out - C*(x/b)**(m-2) *np.cos(m*y_out))**2
    
#     def func2(y_out):
#         return (-y_old + y_out + (2.0*np.pi/q) + aa*np.cos(y_out))**2
    
#     y_old = copy(y)
#     y_old = (newton_krylov(func,y))
#     y_old = np.mod(y_old,2.0*np.pi)

#     x_old = x + (m*b*C)/(m-1)*(x/b)**(m-1) *np.sin(m*y_old)

#     q = q0*(x_old/a)**2

#     y_old2 = copy(y_old)
#     y_old2 = (newton_krylov(func2,y_old))

#     #y_old2 = y_old - 2*np.pi/q #- aa*np.cos(
#     y_old2 = np.mod(y_old2,2.0*np.pi)
#     x_old2 = x_old*(1.0 -aa*np.sin(y_old2))



#     return x_old2,y_old2


def to_index_coord(x,y,nx,ny):
    x_i = x*(nx/(2*np.pi))+(nx/2.)
    y_i = y*(ny/(2*np.pi))  
    
    return x_i,y_i
    
def StandardMap(x,y,M,L,K,q0,b=30.0,aa=0.0,eps=.3,
                m =3,a=40,R_0 = 90):
   
    print('b: ', b)
    #aa = -.00
    B_0 = 1.0 #in tesla0
    #b = 30 #minor rad
    #R_0 = 90 #major rad
    #m = 3. #external mode
    l = 10 #coil width
    #a= 40
    beta = 2.0
    mu = 1.0
    beta_p = beta*(mu+1)/(beta+mu+1)
    
 
    #with PyCallGraph(output=GraphvizOutput()):
    hit_divert = (x>b)
    inCORE = x<b
    stopevolve = hit_divert

    x_new = (stopevolve == 0)*x/(1-aa*np.sin(y))

    # q = (q0*(x_new/a)**2) #*(1 -(1 + beta_p * (x_new/a)**2)*
    #  (1.0- (x_new/a)**2)*(a>x_new))**-1

    xx = x_new/a
    #nu = 1.

    q = q0*(xx**2)/(1-(1+2*xx**2)*(1-xx)**(mu+1))

    print('q: ', q.min(),q.max())
    y_new =  (stopevolve == 0)*(y+ 2*np.pi/q + aa*np.cos(y))
    y_new = np.mod(y_new,2*np.pi)


    # r2rule = -r1+r2 + (m*C*eps*b/(m-1))*(r2/b)**(m-1) * sin(m*th1)
    # th2rule = -th2 + th1 - C *eps *(r2/b)**(m-2) * cos(m*th1)th1rule = -th1 + th0 + aa*cos(th0)+ 2*pi/q


    #


    #see  "DIFFUSIVE TRANSPORT THROUGH A NONTWIST BARRIER IN TOKAMAKS"
    #eps = .3
    print(m,l,a,R_0,q0,b,eps)
    C = ((2*m*l*a**2)/(1.0*R_0*q0*b**2))*eps
    print('C: ', C/eps)


    def func(x_out):
        return (-x_new + x_out +((m*b*C)/(m-1))*((x_out/b)**(m-1) )*np.sin(m*y_new))**2

    x_new2 = copy(x_new)

    #x_new2 =  (stopevolve == 0)* (newton_krylov(func,x_new2)) + (stopevolve)*x

    x_new2 =  (stopevolve == 0)* (newton_krylov(func,x_new2,method='bicgstab')) + (stopevolve)*x
    #print (-x_new + x_new2 +(m*b*C)/(m-1)*((x_new2/b)**(m-1) )*np.sin(m*y_new))**2

    y_new2 = (stopevolve == 0)*(y_new - C*(x_new2/b)**(m-2) * np.cos(m*y_new))+ (stopevolve)*y

    #print 'xchange:', x_new2/x
    #print(x_new2,x)
    #x_new = x_new2
    y_new2 = np.mod(y_new2,2*np.pi)

    stay_inCORE = (x_new2 <b) 
    #full_orbit = (new_inSOL & inSOL) == False #can visit SOL,but can't stay
    full_orbit = (stay_inCORE & inCORE) == True 
    half_orbit = (stay_inCORE  & inCORE) == True #ok, you hit the divertor

    L = L + (stopevolve ==0)*((full_orbit)*q *R_0* 2*np.pi)# + \
                               #half_orbit *100* 2*q* np.pi)
 
    
    Kf = lambdify((r0,th0,r1,th1,r2,th2),K)
    numK = np.array([np.matrix(Kf(x[i],y[i],x_new[i],y_new[i],x_new2[i],y_new2[i])) 
                     for i,elem in enumerate(x_new)])

    # print numK[0].shape
    # exit()



    M_new = np.array([numK[i] * eleM for i,eleM in enumerate(M)])

    return x_new2,y_new2,M_new,L


 
def setup_xz(nx=128,nz=128,b = .3,edge=None,rmin=0.0,rmax=1.0):

    x, z = np.mgrid[b*rmin:b*rmax:complex(0,nx),0:2*np.pi:complex(0,nz)]

    return x,z

def StandardLength(x,z,M,lmbda,K,procnum,return_dict,return_lam,k=1,
                   max_pol_orbits=100,q=5.0,b=50,
                   aa=0.0,eps = .3,m=3,a=40,R_0 = 90):
   
    keep_i= list(np.where(x < b))
    
    
    
    L = 0.0*z
    count = 0

 
    while count<max_pol_orbits and len(x[keep_i]) !=0:
        #print count, x and y get replaced
        x[keep_i], z[keep_i],M[keep_i],L[keep_i] = StandardMap(x[keep_i],z[keep_i],M[keep_i],L[keep_i],K,q,b=b,aa=aa,eps=eps,m = m)
        #with PyCallGraph(output=GraphvizOutput()): #
        #    x[keep_i], z[keep_i],M[keep_i],L[keep_i] = StandardMap(x[keep_i],z[keep_i],M[keep_i],L[keep_i],K,q,b=b,aa=aa,eps=eps,m = m)
        #exit()                  # 
            
        lmbda[keep_i] = np.max(np.real(np.log([np.linalg.eig(np.inner(eleM.T,eleM)**(1./(2.*(count+1))))[0] for eleM in M[keep_i]])),axis=1)

       
       
        # if lyapunov:
        #     x[keep_i], z[keep_i],L[keep_i] = StandardMap(x[keep_i],z[keep_i],L[keep_i],k,q,b=b,aa=aa,eps=eps)  
        
        if procnum==0:
            print (count,' : ',100.0*len(x[keep_i])/np.size(x),'% of the field-lines jumping')
        
    
        keep_i= list(np.where((x < b)))
        count+=1

    R_0 = 85.
    L[keep_i] = count*10*R_0*2*np.pi #make the trapped ones MMUCH longer

    return_dict[procnum] = L
    return_lam[procnum] = lmbda
  
    

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

b = 68.
rmax = 1.1
rmin = .92

#R_0 = 90

manager = Manager()
return_L = manager.dict()
return_lam = manager.dict()

x,z = np.mgrid[b*rmin:b*rmax:complex(0,nx),0:2*np.pi:complex(0,nz)]

lmbda = np.zeros((nx,nz))

I = np.identity(2)
M = np.array([[np.matrix(I) for zz in range(nz)] for xx in range(nx)])



r1n = r1/a
C = ((2*m*l*a**2)/(1.0*R_0*q0*b**2))

r1rule= r1 - r0/(1.0 - aa*sin(th0))
qq = q0*(r1n**2)/( 1- (1+bprime*r1n**2) * (1-r1n**2)**(mu+1) )
th1rule = 2*np.pi/qq -th1 + th0 + aa*cos(th0) #+ 2*np.pi/qq
r2rule = -r1+r2 + (m*C*eps*b/(m-1))*(r2/b)**(m-1) * sin(m*th1)
th2rule = -th2 + th1 - C *eps *(r2/b)**(m-2) * cos(m*th1)
allrules = [r1rule,th1rule,r2rule,th2rule]

    
alleqns = [diff(elem.subs(perturb).series(d,n=2),d).subs(Order(d),0) for elem in allrules]
alleqns2 = [elem.subs(dr1,solve(alleqns[0],dr1)[0])for elem in alleqns[1:4]]
alleqns3 = [alleqns2[k].subs(dth1,solve(alleqns2[0],dth1)[0])for k in [1,2]]
    
r2soln = solve(alleqns3[0],dr2)[0]
soln = [r2soln,solve(alleqns3[1].subs(dr2,r2soln),dth2)[0]]
K = sp.Matrix(soln)
K = K.jacobian([dr0,dth0])


for i in range(Ncpu):
    print(i)
    #p = Process(target=simple, args =(i,return_L))
    #with PyCallGraph(output=GraphvizOutput()):
    p = Process(target=StandardLength, 
                args =(x[chunk*i:chunk*(i+1),:],
                       z[chunk*i:chunk*(i+1),:],
                       M[chunk*i:chunk*(i+1),:],
                       lmbda[chunk*i:chunk*(i+1),:],K,
                       i,return_L,return_lam),
                kwargs = {'q':3,'k':1.5,'max_pol_orbits':10.,
                          'b':b,'eps':.2,'m':2,'aa':-.01,'R_0':85.,
                          'a':a})
    # p = Process(target=StandardLength, 
    #             args =(x[chunk*i:chunk*(i+1),:],z[chunk*i:chunk*(i+1),:],i,return_L),
    #             kwargs = {'q':3,'k':1.5,'max_pol_orbits':10,
    #                       'b':b,'eps':.07,'m':7,'aa':-.02})

    processlist.append(p)
    p.start()
    
for p in processlist:
    p.join()


L = np.concatenate(return_L)
lmbda = np.concatenate(return_lam)

#print(return_L[3].shape)
#print(x[chunk*1:(chunk)*2,:].shape)


pp = PdfPages('sm.pdf')
fast2Dplot(pp,(L),extent=[rmin,rmax,0,2*np.pi])
pp.close()

pp = PdfPages('lam.pdf')
fast2Dplot(pp,lmbda,extent=[rmin,rmax,0,2*np.pi])
pp.close()


                #kwargs = {'q':3,'k':1.5,'max_pol_orbits':2, # 
                #          'b':b,'eps':.07,'m':7,'aa':-.02}) # 
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
