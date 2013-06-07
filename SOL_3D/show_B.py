#!/opt/apps/python/epd/7.2.2/bin/python
import sys, os
HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')

sys.path.append('/usr/local/pylib')
sys.path.append(HOME+'/lib/python')
sys.path.append(HOME+'/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib/boutdata')
sys.path.append(BOUT_TOP+'/tools/pylib/boututils')
sys.path.append(BOUT_TOP+'/tools/pylib/post_bout')

from read_inp import metadata
import sys
import os
from boutdata import collect
from boututils import *
from post_bout import read_grid
import numpy as np
#from plot_CM import CM_mass, present
import pickle
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
from post_bout import read_inp,parse_inp,read_cxx
from read_cxx import *
import matplotlib 
matplotlib.use('pdf')

from reportlab.platypus import *
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch
from reportlab.graphics.charts.linecharts import HorizontalLineChart
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.lineplots import LinePlot
from reportlab.graphics.widgets.markers import makeMarker
from reportlab.lib import colors

from matplotlib.ticker import ScalarFormatter, FormatStrFormatter, MultipleLocator

#data_dir='/tmp/SOLblob/data_3D'
data_dir=sys.argv[1]

#data_dir='/scratch/01523/meyerson/BOUT_sims/SOL_3D/data_parallel_2e-3'
#output_grid.nc
def get_IC(file='output_grid.nc'):
    IC = read_grid(gridfile=file)
    

    meta={}
    
    for elem in IC.variables:
        print elem
        meta[elem] = np.array(IC.variables[elem][:])
  
    return meta
################
#get metadata from  inp, grid file, output nc as well as the source code

#grid
IC_rmp = get_IC() #grid file, refered to as IC here, contains dx,dy and 
# and unnormalized B field component

#BOUT.inp
inp = read_inp(path='./',boutinp='BOUT_3D.inp')
inp = parse_inp(inp)

#some data from a single .dmp.xx.nc
outinfo = file_import(data_dir+'/BOUT.dmp.0.nc') #output data
availkeys = np.array([str(x) for x in outinfo])

#grab some info that is computed at runtime info
metric = ['g11','g22','g33','b0xcvx','b0xcvy','b0xcvz']
for elem in metric:
    IC_rmp[elem] = collect(elem,path=data_dir)
try:
    print path
    cxxinfo = no_comment_cxx(path=data_dir,boutcxx='physics_code.cxx.ref')
       #evolved = get_evolved_cxx(cxxinfo)
    fieldkeys = get_evolved_cxx(cxxinfo)
    fieldkeys = ['['+elem+']' for elem  in fieldkeys]
except:
    print 'cant find the cxx file'

#appent metadata dict to contain some scalar values
for elem in availkeys:
    data =  np.array(outinfo[elem][:])
    if data.shape == (1,):
        IC_rmp[elem]=data[:]
        print elem,data[:] 

#################
#COMPUTE SOME DERIVED QUANTITIES
ny = IC_rmp['ny'][0]
nx = IC_rmp['nx'][0]


#as of now inp, IC.nc and some OUT.nc data is here
#normalize a few things
rho_s= IC_rmp['rho_s']
IC_rmp['Rxy'] = IC_rmp['Rxy']/rho_s
IC_rmp['hthe'] = IC_rmp['hthe']/rho_s
IC_rmp['Zxy'] = IC_rmp['Zxy']/rho_s
IC_rmp['dx'] = IC_rmp['dx']/(rho_s**2 * IC_rmp['bmag'])#dy is already unitless/rescaled
dx = IC_rmp['dx']

#calculate some derived quantities
nz = np.int(inp['[main]']['MZ'])-1

zmax = np.double(inp['[main]']['ZMAX'])
L_z = IC_rmp['Rxy']*2*np.pi*zmax



IC_rmp['L_z'] = L_z
IC_rmp['dz'] = L_z/nz
dz =  L_z/nz
IC_rmp['intdx'] =np.cumsum(dx,0)

r = IC_rmp['Rxy']
hthe = IC_rmp['hthe']
#z = IC_rmp['Zxy'][x_i,:]
dy = IC_rmp['dy']
bt = IC_rmp['Btxy']
bp = IC_rmp['Bpxy']
dlthe = dy*hthe
Lpar = np.cumsum((bt/bp)*dlthe,1)
#print np.repeat(np.min(r,1),)
#print (np.repeat(np.min(r,1),ny)).reshape(nx,ny)
L_perp = np.sqrt((np.abs(r-(np.repeat(np.min(r,1),ny)).reshape(nx,ny)))**2)
L_perp = np.max(r,0)-np.min(r,0) #ny long
L_perp = np.sqrt(np.max(L_z,0)**2 + L_perp**2) #ny
#L_perp = r-np.transpose(np.repeat(np.min(r,0),nx).reshape(ny,nx))
IC_rmp['Lpar'] = Lpar 
IC_rmp['L_perp'] = L_perp  



#IC_rmp['b0xcvz_n'] = IC_rmp['b0xcvz']/(2*np.pi*zmax/nz)
IC_rmp['b0xcvz_n'] = IC_rmp['b0xcvz']*(r*bp/bt)
IC_rmp['b0xcvy_n'] = IC_rmp['b0xcvy']*hthe
IC_rmp['b0xcvx_n'] = IC_rmp['b0xcvx']/(r*bp)


#let's average along the field line
#beta = (np.repeat(np.mean(beta,1),ny)).reshape(nx,ny)
#IC_rmp['beta'] = beta
###################
#SHOW AND PRINTx

for elem in IC_rmp:
    data =  np.array(IC_rmp[elem][:])
    if data.shape == (1,):
        print elem,data[:] 
print IC_rmp.keys()
    




colors = ['b','g','r','c','m','y','k','b','g','r','c','m','y','k']
styles = ['s','^']

pp = PdfPages('grid.pdf')
#list of field to plot and label


field_list =['Bpxy','Btxy','Bxy','dy','dx','dz','intdx','Rxy','hthe','L_z','Lpar',
             'b0xcvx','b0xcvy','b0xcvz','Ni0','psixy','g11','g22','g33','b0xcvz_n']

psixy = IC_rmp['psixy'][:]

#plot some stuff
hi_low=['outer','inner']
for elem in field_list:
    fig = plt.figure()
    canvas = fig.add_subplot(1,1,1)
    lines = canvas.plot(psixy[:,ny/2],IC_rmp[elem][:,ny/2],
                        psixy[:,0],IC_rmp[elem][:,0],
                        label=elem,alpha = .8,markersize=2)
    plt.setp(lines, linewidth=2)

    handles, labels = canvas.get_legend_handles_labels()
    leg = canvas.legend(handles,hi_low,ncol=2,loc='best',prop={'size':12},fancybox=True) 
    leg.get_frame().set_alpha(0.3)

    formatter = ticker.ScalarFormatter()
    formatter.set_powerlimits((-2, 2))  #force scientific notation
    canvas.yaxis.set_major_formatter(formatter)

    canvas.grid(True,linestyle='-',color='.75')
    canvas.set_title(elem)
    
    fig.savefig(pp, format='pdf')
    plt.close(fig)

#lets see how curvature varies along the field line
alongB = ['bxcvx','bxcvy','bxcvz']
fig = plt.figure()
canvas = fig.add_subplot(1,1,1)
for elem in alongB:
    line, = canvas.plot(IC_rmp[elem][nx/2,:],label=elem,alpha = .8,markersize=2)
    canvas.grid(True,linestyle='-',color='.75')

handles, labels = canvas.get_legend_handles_labels()
leg = canvas.legend(handles,labels,ncol=2,loc='best',prop={'size':8},fancybox=True) 
leg.get_frame().set_alpha(0.3)

canvas.set_title('b0xcv along the field line')
plt.setp(line, linewidth=2, color='r')
fig.savefig(pp, format='pdf')
plt.close(fig)

alongB = ['b0xcvx_n','b0xcvy_n','b0xcvz_n']
fig = plt.figure()
canvas = fig.add_subplot(1,1,1)
for elem in alongB:
    line, = canvas.plot(IC_rmp[elem][nx/2,:],label=elem,alpha = .8,markersize=2)
    canvas.grid(True,linestyle='-',color='.75')

handles, labels = canvas.get_legend_handles_labels()
leg = canvas.legend(handles,labels,ncol=2,loc='best',prop={'size':8},fancybox=True) 
leg.get_frame().set_alpha(0.3)

canvas.set_title('b0xcv/dx along the field line')
plt.setp(line, linewidth=2, color='r')
fig.savefig(pp, format='pdf')
plt.close(fig)


#let's compare to curvature computed at run-time - carefull with normalization

#plot the perp scale length as a funtion of parallel displacement
fig = plt.figure()
canvas = fig.add_subplot(1,1,1)
line, = canvas.plot(L_perp,label=elem,alpha = .8,markersize=2)
canvas.grid(True,linestyle='-',color='.75')

handles, labels = canvas.get_legend_handles_labels()
leg = canvas.legend(handles,labels,ncol=2,loc='best',prop={'size':8},fancybox=True) 
leg.get_frame().set_alpha(0.3)

canvas.set_title('L_perp along the field line')
plt.setp(line, linewidth=2, color='r')
fig.savefig(pp, format='pdf')
plt.close(fig)


# field line with the polarization vector as we follow it
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')

x_locs=[nx/2]

for x_i in x_locs:
    r = IC_rmp['Rxy'][x_i,:]
    hthe = IC_rmp['hthe'][x_i,:]
    z = IC_rmp['Zxy'][x_i,:]
    dy = IC_rmp['dy'][x_i,:]
    bt = IC_rmp['Btxy'][x_i,:]
    bp = IC_rmp['Bpxy'][x_i,:]
    dlthe = dy*hthe
    
    phi = np.cumsum((bt/bp)*(dlthe/r))
    x = r * np.sin(phi)
    y = r * np.cos(phi)
    ax.plot(x, y, z, label='a field line')
    #ax.plot(x, y, 0*z+np.min(z))
    ax.plot(0*x+np.min(x), y, z,alpha=.2)
    ax.plot(x, 0*y+np.min(y), z,alpha=.2)
 

ax.legend()


fig.savefig(pp, format='pdf')
plt.close(fig)




pp.close()

             
