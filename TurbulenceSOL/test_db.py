import os, sys, inspect,types,hashlib
import sqlite3 as sql
import pickle as pkl
from datetime import datetime
from pymongo import Connection, MongoClient
from pymongo import errors as mongoErr
from bson.binary import Binary
import cPickle as pickle



HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')
utc = datetime.utcnow()

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import matplotlib.artist as artist 
import matplotlib.ticker as ticker

sys.path.append('/usr/local/pylib')
sys.path.append(HOME+'/lib/python')
sys.path.append(HOME+'/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib')
sys.path.append(BOUT_TOP+'/tools/pylib/boutdata')
sys.path.append(BOUT_TOP+'/tools/pylib/boututils')
sys.path.append(BOUT_TOP+'/tools/pylib/post_bout')

cmd_folder = HOME+'/BOUT_sims/blob_py'

print(cmd_folder)

if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)



from turb_info import field_info as sim 

import numpy as np
from boutdata import collect
import argparse
from frame import Frame, FrameMovie

def is_numeric_paranoid(obj):
    try:
        obj+obj, obj-obj, obj*obj, obj**obj, obj/obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        return True
   
def serialize_obj(obj_in):
     obj_out = {}

     for key,value in obj_in.iteritems():
          #if type(value) not in copy_types:
          if is_numeric_paranoid(value):# or type(value) in copy_types:
               #print key, 'not serializing',value,type(value)
               if 'all' in dir(value):
                    obj_out[key] = np.asscalar(value)
               else:
                    obj_out[key] = value
               
          elif type(value) in copy_types:
                #print key, type(value)
                obj_out[key] = value

          elif type(value) is types.ListType:
                #print key, type(value)
                obj_out[key] = value

          elif type(value) is types.DictType:
               #print key, type(value),value.keys()
               obj_out[key] = value
          
          else:
               print key, 'serializing',value,type(value)
               obj_out[key] = Binary(pickle.dumps(value,protocol=2))

     return obj_out

parser = argparse.ArgumentParser()


parser.add_argument("-s","--server", type=str,
                    help="server url",default='128.83.61.211')

parser.add_argument("--path", type=str,
                    help="data location",nargs='?')
                    #default=SCRATCH+'/BOUT_sims/AlphaChaos/data_chaos_a1.0e-2_eps1.0e-1')

parser.add_argument("--tchunk", type=int,
                    help="time chunk",nargs='?', default=10)




args = parser.parse_args()

path = args.path
tchunk = args.tchunk
server = args.server

#print path

#CONNECT TO DATABASE
local_c = MongoClient(host='localhost')
c = MongoClient(host=server)
db = c.test_database

#copy_database(from_name, to_name[, from_host=None

#local_c.copy_database('alpha_runs','beavis_copy',from_host=server)
#local_c.copy_database('test_database','beavis_copy',from_host=server)


db_local = local_c['beavis_copy']

print type(local_c)
print type(db_local),dir(db_local)

print db_local.beavis_copy.find({"author":"Dmitry"}).count()

#c.drop_database(db)
alpha_runs = db.alpha_runs
# for run in alpha_runs.find({"author": "Dmitry"}):
#      print 'how many results: ', len(run)
# print run

if path:
     n = np.squeeze(collect("n",path=path,tind =[0,0]))
     # u = np.squeeze(collect("u",path=path,tind =[1,50]))
     # phi = np.squeeze(collect("phi",path=path,tind =[1,50]))
     nx,ny = n.shape

     meta=[]
     #blob_db = blob_info(n)

     time = np.squeeze(collect("t_array",path=path,xind=[0,0]))


     sim_blob = {"author": "Dmitry",
                 "text": "My 2nd blog post!",
                 "tags": ["fusion", "in", "50 years"],
                 "date": datetime.utcnow()}



     copy_types = [types.StringType,types.LongType,types.FloatType,types.IntType,type( datetime.utcnow())]
     ban_types = [type(np.array([123]))]#types.DictType]


     meta={'y0':0.0,'x0':0.0,'dx':1.0,'dy':1.0,'location':path}
     blob = sim(path,meta=meta,debug=False)

     for key,value in blob.__dict__.iteritems():
          print key, type(value)
          #print 'take data from sim to a big dictionary'
          if type(value) not in ban_types:
               #print key, type(value)
               sim_blob[key] = value
          else:
               if type(value) == type(np.array([12])):
                    #print value.size
                    if value.ndim == 1:
                         sim_blob[key] = value.tolist()
                    elif value.size == 1:
                         sim_blob[key] = value
                    else:          
                         print 'not adding to the db obj',key, type(value),value.shape
     
     ser_dict = serialize_obj(sim_blob)
     alpha_runs.ensure_index('md5',unique=True,dropDups=True)#,{'unique': True, 'dropDups': True})
     try:
          alpha_runs.insert(ser_dict,db)
     except mongoErr.DuplicateKeyError:
          print 'the hash: ',ser_dict['md5'],sim_blob['md5'],ser_dict.keys()

          ser_dict.pop("_id",None) #rip off the id 

          alpha_runs.update({'md5':ser_dict['md5']}, {"$set":ser_dict})#, upsert=False)
          print 'Duplicate run not adding to db, but updating'
     #print blob.__dict__.keys()

     
     f=open(path+'/BOUT.log.0', 'r')
     pathprint = hashlib.md5(f.read()).hexdigest()
     f.close()


     for run in alpha_runs.find({"md5": pathprint}):
          print 'how many results: ', len(run)



 #print alpha_runs.find({"md5": pathprint}).count()
print alpha_runs.find({"author":"Dmitry"}).count()




#print foo 
#


for run in alpha_runs.find({"author": "Dmitry"}):
     print 'how many results: ', len(run)

      #for key in run:
      #     print key, run[key].__class__

print db.collection_names()
pp = PdfPages('sm.pdf')

compare = plt.figure()

i = 1
lamdas = []


def magic(numList):
        s = ''.join(map(str, numList))
        return int(s)

lam_chaos = []
lam_smooth = []
eps=[]
eps_smooth = []
jump = []

for run in alpha_runs.find():
     if ("location" in run) and ("eps" not in run):
          print 'updating eps field'
          eps_i = run['location'].find('eps')
          chaos = run['location'].find('chaos')
          jump = run['location'].find('jump')
          smooth = run['location'].find('smooth')
          if eps_i > 1: 
               run['eps'] = np.float((run['location'][eps_i::]).lstrip('eps').rstrip('/'))
          else:
               run['eps'] = 0
          if jump >1:
               run['tags'].append('jump')
          if chaos >1:
               run['tags'].append('chaos')
          #later figure out how to save
          


#############EXTRACT FROM DATABASE AND PLOT #############
fig = plt.figure()
fig_pdf = plt.figure()
for run in alpha_runs.find():
     if "location" in run:
          print run['location'], len(run['lam'])

     if "xmoment" in run:
          print "xmoment", len(run['xmoment']['1'])
          
     

     print run.keys()
     if "pdf_y" in run:
          print 'histograms . . .'
          
          y = pickle.loads(run['pdf_y'])
          # x = pickle.loads(run['pdf_x'])
          x = run['pdf_x']

          print y[0].shape,len(y),len(x)

          pdf_frm = Frame(y[-1][600,:],meta={'ylabel':r'$\rho_s$',
                                  'xlabel':r'$t$',
                                  'title':r'$\lambda$'+' for '+ run['path'][-20:-1],
                                  'ticksize':30,'fontsz':10,
                                  'stationary':True})
          pdf_frm.x = x[-1][600,:]
          pdf_frm.render(fig_pdf,magic([1,1,1]))
          fig_pdf.savefig(pp, format='pdf')

     if "lam" in run and len(run['lam']) > 10 :
     
          temp = (np.array(run['lam'])).transpose()
          #lam_frame = Frame
          #canvas.plot(temp[0],temp[1])
          if len(temp) == 2:
               y = temp[1]
               dx = temp[0][1]- temp[0][0]
          else:
               y = run['lam']
               dt = run['t']
               dt = dt[1]- dt[0]


          lam_frm = Frame(y,meta={'dx':dt,
                                        'ylabel':r'$\rho_s$',
                                        'xlabel':r'$t$',
                                        'title':r'$\lambda$'+' for '+ run['path'][-20:-1],
                                        'ticksize':30,'fontsz':10,
                                        'stationary':True})

          
          if i < 10:
               lam_frm.render(fig,magic([3,3,i]))
               lam_frm.ax.set_ylim(0,1.3*np.max(temp[1]))

          if len(y) > 10:
               eps_i = run['location'].find('eps')
               chaos = run['location'].find('chaos')
               jump = run['location'].find('jump')
               smooth = run['location'].find('smooth')

               if eps_i > 1: 
                    if jump >1:
                         eps_smooth.append(0)
                         lam_smooth.append(np.mean(y[-40:-1]))
                    elif smooth > 1:
                         eps_smooth.append(np.float((run['location'][eps_i::]).lstrip('eps').rstrip('/')))
                         lam_smooth.append(np.mean(y[-40:-1])) 
                    elif chaos >1:
                         eps.append(np.float((run['location'][eps_i::]).lstrip('eps').rstrip('/')))
                         lam_chaos.append(np.mean(y[-40:-1]))
                         
          #lam_frm.ax = None
          #lam_frm.render(compare,111)

          
          i = i+1
          #print run['location'].find('eps')
          #print temp[1]
         
          #canvas.plot(run['time'][0:len(run['lam'])],run['lam'])
     #sm.plot(run['xmoment'])
     #sm.plot(run['time'],run['xmoment']['1'])
     #print len(run['time']), len(run['xmoment']['1'])


lam_meta = {'style':'.','markersize':20,
            'stationary':True,'ylabel':r'$\rho_s$',
            'xlabel':r'$\epsilon$','fontsz':20, 
            'title':r'$\lambda_n$'+r', $\alpha=1e-2, \beta= 1e-2$',
            'alpha':.5}

lam_compare  = Frame(lam_smooth,meta=lam_meta)
lam_compare.x = eps_smooth

lam_chaos =  Frame(lam_chaos,meta=lam_meta)
lam_chaos.x = eps
lam_chaos.style = 'r.'
lam_chaos.alpha = .8
#lam_chaos =  Frame(lam_chaos,meta={'x':eps,'style':'r.','alpha':.8})



lam_compare.render(compare,111)
lam_chaos.ax = lam_compare.ax
lam_chaos.stationary=True
lam_chaos.render(compare,111)
#lam_chaos.ax.legend(loc='upper left')
leg = lam_chaos.ax.legend([lam_compare.img,lam_chaos.img],['w/out chaos','w/ chaos'],
                    loc='best', fancybox=True)   
leg.get_frame().set_alpha(0.5)        

compare.savefig(pp,format='pdf')
fig.savefig(pp, format='pdf')

pp.close()

# for elem,val in enumerate(post):
          #  print elem,post,post[val]
