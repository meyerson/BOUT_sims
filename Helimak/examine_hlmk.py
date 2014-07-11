#!/opt/apps/python/epd/7.2.2/bin/python
import sys, os, gc
from datetime import datetime
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

from pymongo import Connection, MongoClient
from pymongo import errors as mongoErr
from bson.binary import Binary
import cPickle as pickle
import zlib
import numpy as np
server ='dataplay.no-ip.org'


c = MongoClient(host=server)
db = c.hlmk
db.authenticate("lonestar","britoebito")
n_coll = db.n_profiles
phi_coll = db.phi_profiles
Te_coll = db.Te_profiles
v_coll = db.V_profiles
u_coll = db.u_profiles
meta_coll = db.meta



#use the meta collection to isolate a subset or records
#specfically grab the last set of  documents the corresponds to each run
#temp = meta_coll.find()
def to_dict(record):
    return dict_obj

meta = meta_coll.find({"nx": {"$exists": True}})
#meta =  meta_coll.find()
ids = list(record['_id'].encode('ascii') for record in meta)

#pull and group the corresponding profiles
n_profiles = n_coll.find({"_id":{"$in":ids}})
phi_profiles = phi_coll.find({"_id":{"$in":ids}})

for elem in ids:
    meta =  meta_coll.find({"_id":elem})
    n_profile = n_coll.find({"_id":elem})
    phi_profile = phi_coll.find({"_id":elem})
    nx = meta[0]['nx']
    nt = meta[0]['nt']
    n = np.array(n_profile[0]['xt'] ).reshape(nt,nx)
    dn =  np.array(n_profile[0]['std'] ).reshape(nt,nx)
    print np.max(dn/n)
# n_xt = np.array(n_profiles[0]['xt']).reshape(nt,nx)
