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
import matplotlib.pyplot as plt
from pymongo import Connection, MongoClient
from pymongo import errors as mongoErr
from bson.binary import Binary
import cPickle as pickle
import zlib
import numpy as np
import mdp 

from frame import Frame

server ='dataplay.no-ip.org'
import matplotlib
matplotlib.use('Agg')
from pylab import legend

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

meta = meta_coll.find({"[main]": {"$exists": True}})
n_profile = n_coll.find({"std": {"$size": 364}},{'$_id': 1})

good_ids = map(lambda row:row['_id'],n_profile)
print len(good_ids)
print meta.count()
#good_ids = print n_profile
# for elem in meta:
#     print elem['_id'] in good_ids
# exit()
ids={}
for elem in meta:
    
    bias = elem["[hlmk]"]['bias_phi_0']
    if elem['_id'] in good_ids:
        if bias in ids:
            if elem['nt'] > ids[bias]['nt'] and elem['t1'] != 0 and elem['_id'] in good_ids:
                ids[bias]['nt'] = elem['nt']
                ids[bias]['_id'] = elem['_id']
                ids[bias]['t1'] = elem['t1']
        else:
            ids[bias]={}
            ids[bias]['nt'] = elem['nt']
            ids[bias]['_id'] = elem['_id']
            ids[bias]['t1'] = elem['t1']

# for elem in ids:
#     print elem,ids[elem]['nt']

# exit()

def plot_profiles(key,prof_coll,prefix=None):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for elem in ids:
        _id = ids[elem]['_id']
        print _id, _id in good_ids
        # try:
        meta =  meta_coll.find({"[main]": {'$exists': True},"_id":_id})
        n_profile = prof_coll.find({"_id":_id})
        print 'how many:', n_profile.count(),meta.count()
        if key == 'std':
            try:
                print '###################################'
                print 'bias phi:', meta[0]["[hlmk]"]["bias_phi_0"]
                print meta[0]["t1"], meta[0]["t2"]
                print 'std: ',np.array(n_profile[0][key]).shape
                
                print  np.array(n_profile[0]['xt']).shape
            except:
                pass
            #exit()

        try:
            nx = meta[0]['nx']
            nt = meta[0]['nt']
            
            dx = meta[0]["dx"]
            ## print meta[0]["t1"]
            ## meta[0]["t2"],meta[0].keys()

            n = np.array(n_profile[0]['xt'] )
            value =  np.array(n_profile[0][key] )

            reshape =  value.shape != n.shape
            n = n.reshape(nt,nx)
            if key is 'xt':
                v_mean_norm = np.mean(n,axis=0)
            elif reshape:
                print 'ploting std'
                v_mean_norm = value/np.mean(n,axis=0)
            else:
                value =  value.reshape(nt,nx)
                v_mean_norm = np.mean(value/n,axis=0)

            ax.plot(dx*np.arange(nx),v_mean_norm,label= meta[0]["[hlmk]"]["bias_phi_0"])
            fig.savefig(prefix+key+'.eps')
        except:
            pass
    ax.legend()
    plt.axvspan(nx*dx*.25,nx*dx*.5, facecolor='g', alpha=0.2,fill=False,hatch='/')
    fig.savefig(prefix+key+'.eps')
    plt.axvspan(nx*dx*.25,nx*dx*.5, facecolor='g', alpha=0.2)
    fig.savefig(prefix+key+'.pdf')
    plt.close(fig)   


def V_profiles(key,x_coll,y_coll,prefix=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for elem in ids:
        _id = ids[elem]['_id']
        meta =  meta_coll.find({"[main]": {'$exists': True},"_id":_id})
        y_prof = y_coll.find({"_id":_id})
        x_prof = x_coll.find({"_id":_id})
        #try:
        nx = meta[0]['nx']
        nt = meta[0]['nt']
        
        dx = meta[0]["dx"]
        n = np.array(y_prof[0]['xt'] )
        y = np.array(y_prof[0][key] )/np.mean(n,axis=0)

        x = np.array(x_prof[0]['V_y_ph'] ).reshape(nt,nx)
        x = np.mean(x,axis=0)
        #x = np.gradient(np.mean(x,axis=0))/dx
        print x.shape

        ax.plot(x[80:180],y[80:180],'.',label= meta[0]["[hlmk]"]["bias_phi_0"],marker='.') 
        fig.savefig(prefix+key+'.eps')
        #ax.scatter(x,y,label= meta[0]["[hlmk]"]["bias_phi_0"])    
        

    ax.legend()
    fig.savefig(prefix+key+'.eps')
    plt.close(fig)   


V_profiles('std',v_coll,n_coll,prefix='Image/n_')
     

exit()    
for key in ['xt','std','min',"max"]:
    plot_profiles(key,n_coll,prefix='n_')
# for key in ['xt', 'std','min',"max"]:
#     plot_profiles(key,phi_coll,prefix='phi_')
# for key in ['xt', 'std','min',"max"]:
#     plot_profiles(key,Te_coll,prefix='Te_')
# for key in ['xt', 'std','min',"max"]:
#     plot_profiles(key,u_coll,prefix='u_')






# n_xt = np.array(n_profiles[0]['xt']).reshape(nt,nx)
