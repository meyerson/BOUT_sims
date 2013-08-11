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

parser.add_argument("--debug", type=bool,
                    help="debug",nargs='?', default=False)


args = parser.parse_args()

path = args.path
tchunk = args.tchunk
server = args.server
debug = args.debug

#print path

#CONNECT TO DATABASE
c = MongoClient(host=server)
db = c.test_database
#c.drop_database(db)
alpha_runs = db.alpha_runs


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
     blob = sim(path,meta=meta,debug=debug)

     print dir(blob)
     blob.to_db(server=server)
    
print alpha_runs.find({"author":"Dmitry"}).count()




#print foo 
#


for run in alpha_runs.find({"author": "Dmitry"}):
     print 'how many results: ', len(run)

      #for key in run:
      #     print key, run[key].__class__
for run in alpha_runs.find():
     for key in run.keys():
          print key, run[key].__class__

print db.collection_names()
pp = PdfPages('sm.pdf')
#fig, canvas = plt.subplots(1)
fig = plt.figure()
    

i=0
for run in alpha_runs.find():

     if "lam" in run:
          sm = fig.add_subplot(2,2,i+1)
          temp = (np.array(run['lam'])).transpose()
          sm.plot(temp[0],temp[1])
          i = i+1
          #canvas.plot(run['time'][0:len(run['lam'])],run['lam'])
     #sm.plot(run['xmoment'])
     #sm.plot(run['time'],run['xmoment']['1'])
     try:
          print len(run['time']), len(run['xmoment']['1'])
     except:
          print 'nevermind'
fig.savefig(pp, format='pdf')
pp.close()

# for elem,val in enumerate(post):
          #  print elem,post,post[val]
