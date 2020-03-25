#!/usr/bin/env python
#coding: utf-8

import os
import sys
from StarRW import *

if len(sys.argv) < 3: 
        print "usage: python joinCtfAndAllStack.py TMV/ CtfFind/job002/"
        sys.exit()

os.system("mkdir %s/StarFiles"%sys.argv[1])
l = [i for i in os.listdir(sys.argv[1]) if os.path.isdir("%s/%s"%(sys.argv[1], i)) and not i.startswith("StarFiles")]
ctfdir = sys.argv[2]
ctfstar = Star(ctfdir+"/micrographs_ctf.star")

primdir = os.path.abspath(sys.argv[1])

for f in l:
        os.chdir(primdir)
        wdir = f
        os.chdir(wdir)
        ctfparams = ctfstar.get_parameter_names()
        ctfnameindex = ctfstar.find_parameter_index("_rlnMicrographName")
        indexes = range(len(ctfparams))
        stars = [Star(f) for f in os.listdir(os.getcwd()) if ((not (f.endswith("preprocessed.star") or f.endswith("particles.star"))) and ("stack" in f) and f.endswith(".star"))]
        if len(stars) == 0: continue
        newheader = stars[0].header+ctfparams
        newcontent = []
        for star in stars:
                nameindex = star.find_parameter_index("_rlnMicrographName")
                for i in range(len(star.content)):
                        for ctfcont in ctfstar.content:
                                if os.path.basename(ctfcont[ctfnameindex]) == os.path.basename(star.content[i][nameindex]):
                                        additional = [ctfcont[index] for index in indexes]
                                        star.content[i] += additional
                newcontent += star.content
        write_star(wdir.rstrip("/")+".star", newheader, newcontent)
        os.system("mv %s.star ../StarFiles/"%wdir.rstrip("/"))

os.chdir(primdir)
dirlist =  os.listdir("StarFiles/")
dirlist.sort()

stars = []
for filename in dirlist:
        stars.append(Star("StarFiles/"+filename))

content = []
for star in stars:
        micnameindex = star.find_parameter_index("_rlnMicrographName")
        star.remove_this_parameter(micnameindex)
        star.add_parameter("_rlnOriginX")
        star.add_parameter("_rlnOriginY")
        content += star.content

write_star("joined.star", star.header, content)

