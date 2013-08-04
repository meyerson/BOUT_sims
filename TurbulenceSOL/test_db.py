import os, sys, inspect,types
import sqlite3 as sql
import pickle as pkl
from datetime import datetime
from pymongo import Connection, MongoClient
from bson.binary import Binary
import cPickle as pickle


HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')
utc = datetime.utcnow()

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

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

#prepare a list of directories
sim_key='Ra1e4_turb'

#path="/tmp/SOLblob/data_"+sim_key
path=SCRATCH+'/BOUT_sims/AlphaChaos/data_chaos_a1.0e-2_eps1.0e-1'


#read the data and process

n = np.squeeze(collect("n",path=path,tind =[0,0]))
# u = np.squeeze(collect("u",path=path,tind =[1,50]))
# phi = np.squeeze(collect("phi",path=path,tind =[1,50]))
nx,ny = n.shape

# pp = PdfPages('svd0.pdf')
#         #self.canvas.imshow(self.raw_data.reshape(self.nt,self.nx*self.ny)=)
# fig = plt.figure()
# canvas = fig.add_subplot(111)
# canvas.imshow(n[:,:,:].reshape(nt,nx*ny),aspect='auto')
# fig.savefig(pp,format='pdf')
# plt.close(fig)
# pp.close()

# print n.shape
meta=[]
#blob_db = blob_info(n)

meta={'y0':0.0,'x0':0.0,'dx':1.0,'dy':1.0,'location':path}
blob = sim(path,meta=meta)



time = np.squeeze(collect("t_array",path=path,xind=[0,0]))


sim_run  = {"author": "Dmitry",
             "text": "My first blog post!",
             "tags": ["fusion", "in", "40 years"],
             "date": datetime.utcnow()}

sim_run2 = {"author": "Dmitry",
             "text": "My 2nd blog post!",
             "tags": ["fusion", "in", "50 years"],
             "date": datetime.utcnow()}

sim_blob = {"author": "Dmitry",
             "text": "My 2nd blog post!",
             "tags": ["fusion", "in", "50 years"],
             "date": datetime.utcnow()}

# for key,value in sim_blob.iteritems():
#      print key,value
#      blob.__dict__[key] = value

save_to_db = ['md5','location',]

copy_types = [types.StringType,types.LongType,types.FloatType,types.IntType,type( datetime.utcnow())]
ban_types = [type(np.array([123])),types.ListType]#types.DictType]

for key,value in blob.__dict__.iteritems():
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
     

#print blob.__dict__.keys()



def create_db():
     c = MongoClient(host='tselitel.no-ip.org')
     return  c.test_database

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

def push_to_db(obj_dict,db):
     posts = db.posts
     posts.insert(obj_dict)
     #print obj_dict
    # print posts.find_one({"author": "Dmitry"})


#c = MongoClient()
c = MongoClient(host='tselitel.no-ip.org')
db = c.test_database
c.drop_database(db)
print 'after'
#sys.exit("Connected to TSELITEL")
#print blob.fft.shape

# for key,value in blob.iteritems():
#      print key,value, value.__class__

posts = db.posts

ser_dict = serialize_obj(sim_blob)



posts.insert(ser_dict,db)
#push_to_db(sim_run2,db)

for post in posts.find({"author": "Dmitry"}):
      print post.keys()
# for elem,val in enumerate(post):
          #  print elem,post,post[val]
