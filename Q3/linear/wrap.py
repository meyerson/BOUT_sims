#! /share/home/01523/meyerson/local/bin/python
import post_bout as pb
from pb_present import LinResPresent
import sys
s = pb.corral(cached=True,logname=sys.argv[1],debug=False)
y = LinResPresent(s.db)
y.show(pdfname=sys.argv[2])
