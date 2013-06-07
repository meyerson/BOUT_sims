import os, sys, inspect
import sqlite3 as sql
import pickle as pkl
from datetime import datetime

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

print cmd_folder

if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

from blob_info import blob_info, Blob2D
from blob_draw import BlobDraw
from blob_present import BlobPresent
from blob_movie import BlobMovie
from blob_db import BlobDB

import numpy as np
from boutdata import collect

#prepare a list of directories
sim_key='Ra1e4_turb'
path="/tmp/SOLblob/data_"+sim_key



#read the data and process
try:
     n = np.squeeze(collect("n",path=path,tind =[1,55]))
     u = np.squeeze(collect("u",path=path,tind =[1,55]))
     phi = np.squeeze(collect("phi",path=path,tind =[1,55]))
     nt,nx,ny = n.shape
except:
     print "fail"
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

meta={'y0':0.0,'x0':0.0,'dx':1.0,'dy':1.0}
#blob1 = Blob2D(n,meta=meta)
#blob2 = BlobDraw(n,meta=meta)

#create a means of saving a created blob obj
#meta['data2']=phi

#blob = BlobMovie(n,meta=meta)
#vort_movie = BlobMovie(u,meta=meta)
print 'movie'


#blob = BlobDB(n,meta=meta)
#blob = BlobMovie(n,meta=meta)
blob = Blob2D(n,meta=meta)

time = np.squeeze(collect("t_array",path=path,xind=[0,0]))
#how to quickly make a movie
from frame import Frame, FrameMovie
data_c = phi
frm_data = Frame(blob.raw_data,meta={'mask':True,'data_c':data_c})
frm_amp = Frame(blob.amp)

frm_u = Frame(u,meta={'mask':True,'data_c':data_c})

sigma = blob.raw_data.std(axis=2)
frm_data1D = Frame(np.average(blob.raw_data,axis=2),meta={'sigma':sigma,'t_array':time})

frames= [frm_data,frm_u,frm_data1D,frm_amp]
#FrameMovie(frames,fast=True)

#let's examine the boundary
bc_left = Frame(blob.raw_data[:,0:20,:],meta={'mask':True})
bc_right = Frame(blob.raw_data[:,-5:-1,:],meta={'mask':True})

sig_lft = blob.raw_data[:,0:20,:].std(axis=2)
bc1d_left = Frame(np.average(blob.raw_data[:,0:20,:],axis=2),meta={'sigma':sig_lft,'t_array':time})
FrameMovie([bc_left,bc1d_left],fast=False,moviename='boundary',fps=2)

#attempt to use sqlite - fail for now


def to_db(blob,db_name='blob.db'):
     conn = sql.connect(db_name)
     conn.row_factory = sql.Row

     cursor = conn.cursor()

     cursor.execute("""CREATE TABLE albums
                   (title text, artist text, release_date text, 
                    publisher text, media_type text) 
                """)

     # insert some data
     #cursor.execute("INSERT INTO albums VALUES ('Glow', 'Andy Hunter', '7/24/2012', 'Xplore Records', 'MP3')")
     
# save data to database
     conn.commit()
     cursor.execute("SELECT * FROM albums")
          #cur.execute('PRAGMA table_info(blob)')
          
     rows = cursor.fetchall()
     
     print rows
   
     
# insert multiple records using the more secure "?" method
     albums = [('Exodus', 'Andy Hunter', '7/9/2002', 'Sparrow Records', 'CD'),
               ('Until We Have Faces', 'Red', '2/1/2011', 'Essential Records', 'CD'),
               ('The End is Where We Begin', 'Thousand Foot Krutch', '4/17/2012', 'TFKmusic', 'CD'),
               ('The Good Life', 'Trip Lee', '4/10/2012', 'Reach Records', 'CD')]
     cursor.executemany("INSERT INTO albums VALUES (?,?,?,?,?)", albums)
     conn.commit()
     cursor.execute("SELECT * FROM albums")
          #cur.execute('PRAGMA table_info(blob)')
          
     rows = cursor.fetchall()
     
     print rows

def update_db(blob,db_name='blob.db'):
     conn = sql.connect(db_name)
     cursor = conn.cursor()
     
     update_str = """UPDATE albums SET artist = 'John Doe' WHERE artist = 'Andy Hunter'"""
     cursor.execute(update_str)
     conn.commit()     


def add_to_db(blob_info,db_name='blob.db'):


     con = sql.connect(db_name)
     #con.row_factory = sql.Row
     
     with con:
          
          cur = con.cursor()    
          cur.execute('SELECT SQLITE_VERSION()')
          
          data = cur.fetchone()
          
          print "SQLite version: %s" % data                
          
    #some fields/columns/enties in the database
          datelist = {'year':['INTEGER',utc.year],'month':['INTEGER',utc.year],
                      'day':['INTEGER',utc.year],'hour':['INTEGER',utc.year],
                      'minute':['INTEGER',utc.minute]}
          
          
          meta_db = {}
          grid_list =['dx','dy','dz','dt','nt','nx','ny','nz','Lz','Lx','Ly']
          grid_db = {}
          for elem in grid_list:
               grid_db[elem] = ['blob',None]
               
          meta_db.update(grid_db)    
               
          moment_db = {'xmoment':['blob',None],'ymoment':['blob',None],'Vxmax':['REAL',None],
                            'Vymax':['REAL',None],'Vmax':['REAL',None]}

          meta_db.update(moment_db)

          # wavelet_db = {}
          # meta_db.update(wavelet_db)
    
          BC_db = {'xplus':['TEXT',None],'xminus':['TEXT',None],'yup':['TEXT',None],'ydown':['TEXT',None]}
    
          meta_db.update(BC_db)
          
          IC_db  = {'description':['TEXT',None]}
    
          meta_db.update(IC_db)
          
          meta_row = meta_db
 
          cur.execute("CREATE TABLE blob_tbl (Id INTEGER PRIMARY KEY);")
          
         

          #create the databasefrom the meta_row dictionary
          for colm in meta_row:
               print meta_row[colm],colm
               cur.execute("alter table blob_tbl ADD "+colm+" "+meta_row[colm][0]+";")
               con.commit()

               
          cur.execute("INSERT INTO blob_tbl(xplus) VALUES ('Rebecca');")
          cur.execute("INSERT INTO blob_tbl(xplus) VALUES ('Butts');")
           # albums = [('Exodus', 'Andy Hunter', '7/9/2002', 'Sparrow Records', 'CD'),
           #           ('Until We Have Faces', 'Red', '2/1/2011', 'Essential Records', 'CD'),
           #           ('The End is Where We Begin', 'Thousand Foot Krutch', '4/17/2012', 'TFKmusic', 'CD'),
           #           ('The Good Life', 'Trip Lee', '4/10/2012', 'Reach Records', 'CD')]
           # cursor.executemany("INSERT INTO albums VALUES (?,?,?,?,?)", albums)
          con.commit()
                #lid = cur.lastrowid
                #print "The last Id of the inserted row is %d" % lid

     # cur.execute("UPDATE blob SET amp = ? WHERE Id = 1", [sql.Binary(pkl.dumps(blob_db['x']))])


           #cur.execute("INSERT INTO blob (Id, xminus) VALUES(1,'3000AD')")

     #amp = pkl.dumps(np.random.rand(10,20)) #serialize
     #print pkl.loads(amp) #get it back
          cur.execute("SELECT * FROM blob_tbl")
           #cur.execute('PRAGMA table_info(blob)')

          rows = cur.fetchall()

           #print rows
          for row in rows:
               print row
          # print cur.rowcount
          
          # for row in cur.execute("SELECT rowid, * FROM blob"):
          #      print row

#add_to_db(blob)
#update_db(blob)
