#!/usr/bin/python
#coding: utf-8

import os
import sys
from StarRW import *

if len(sys.argv) != 3:
        print "python star2par.py input.star output.par"
        sys.exit()

os.system("python regroup.py %s %s"%(sys.argv[1], "tmp.star"))
print "regrouped."

star = Star("tmp.star")
par = open(sys.argv[2], "w")

apixindex = star.find_parameter_index("_rlnDetectorPixelSize")
psiindex = star.find_parameter_index("_rlnAnglePsiPrior")
thetaindex = star.find_parameter_index("_rlnAngleTiltPrior")
phiindex = star.find_parameter_index("_rlnAngleRot")
xindex = star.find_parameter_index("_rlnOriginX")
yindex = star.find_parameter_index("_rlnOriginY")
magindex = star.find_parameter_index("_rlnMagnification")
gnindex = star.find_parameter_index("_rlnGroupNumber")
df1index = star.find_parameter_index("_rlnDefocusU")
df2index = star.find_parameter_index("_rlnDefocusV")
dangindex = star.find_parameter_index("_rlnDefocusAngle")

cont = []
APIX = star.content[0][apixindex]/star.content[0][magindex]*10000
print "APIX=%f"%APIX

l1 = "C	  PSI	  THETA	 PHI	   SHX	    SHY	      MAG    FILM   DF1	     DF2     ANGAST  OCC 	 LogP	  SIGMA	   SCORE    CHANGE"
cont.append(l1)

filmnum = 0
prevgn = 0
for i in range(len(star.content)):
        num = i+1
        psi = star.content[i][psiindex]
        theta = star.content[i][thetaindex]
        phi = star.content[i][phiindex]
        shx = -star.content[i][xindex]*APIX
        shy = -star.content[i][yindex]*APIX
        #mag = star.content[i][magindex]
        mag = 10000
        if prevgn != star.content[i][gnindex]:
                prevgn = star.content[i][gnindex]
                filmnum += 1
        film = filmnum
        df1 = star.content[i][df1index]
        df2 = star.content[i][df2index]
        dang = star.content[i][dangindex]
        line = "%7d%8.2f%8.2f%8.2f%10.2f%10.2f%8d%6d%9.1f%9.1f%8.2f%7.2f%11d%11.4f%8.2f%8.2f"%(num, psi, theta, phi, shx, shy, mag, film, df1, df2, dang, 100., -500, 1., 20., 0.)
        cont.append(line)        

par.write("\n".join(cont))
par.close()
print "par file created."
os.system("python shorten.py %s %s"%(sys.argv[2], sys.argv[2]))
print "shortened."
os.remove("tmp.star")
