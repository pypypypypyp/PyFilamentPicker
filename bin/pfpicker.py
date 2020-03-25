#!/usr/bin/python
#coding: utf-8

import time
import sys
from EMAN2 import *
from math import *
from progressbar import ProgressBar
from StarRW import *
import optparse
import shutil

def setupParserOptions():
        parser = optparse.OptionParser(usage="Usage: %prog <.in file> [options]")
        parser.add_option("--continue", dest="continue", action="store_true", default=False, help="Continue from continue.in")
        parser.add_option("--safemode", dest="safemode", action="store_true", default=False, help="Continue from safemode.in")
        parser.add_option("--loadnum", dest="loadnum", action="store", type="int", default=50, help="# of micrographs to load")
        parser.add_option("--sigma", dest="sigma", action="store", type="int", default=20, help="sigma value for AGK")
        parser.add_option("--rho", dest="rho", action="store", type="float", default=1.2, help="rho value for AGK")
        parser.add_option("--debug", dest="debug", action="store_true", default=False, help="Save debug files")
        options, args = parser.parse_args()
        if len(args) > 1: parser.error("Unknown command-line options: %s"%str(args))
        if len(sys.argv) < 2:
                parser.print_help()
                sys.exit()
        params = {}
        for i in parser.option_list:
                if isinstance(i.dest, str): params[i.dest] = getattr(options, i.dest)
        return params

def openInFile(infile):
       cont = [cont.rstrip("\n") for cont in open(infile, "r").readlines()]
       if cont[0] == "--all":
                return sorted([cont[1]+"/"+i for i in os.listdir(cont[1])])
       if cont[0] == "--re":
                list = os.listdir(cont[1])
                matches = []
                for name in list:
                        m = re.match(cont[2], name)
                        if m is not None: matches.append(cont[1]+"/"+name)
                return sorted(matches)
       if cont[0] == "--list":
                nofmics = len(cont)-2
                for i in range(nofmics):
                        return sorted([cont[1]+"/"+cont[2+i] for i in range(nofmics)])

def initialize(params):
        # read .in file
        head = sys.argv[1][:-3]
        if params["continue"]: firstin = head+"/continue.in"
        elif params["safemode"]: firstin = head+"/safemode.in"
        else:
                if os.path.exists(head+"/continue.in"):
                        shutil.copyfile(head+"/continue.in", head+"/continue.in.backup")
                firstin = sys.argv[1]
        print "Micrographs in %s will be loaded ..."%firstin
        mrcfiles_initial = openInFile(firstin)
        mrcfiles_initial = ["../"+os.path.normpath(i) for i in mrcfiles_initial]
        # set up directory & file structure
        if not os.path.exists(head): os.system("mkdir %s"%head)
        if not os.path.exists("%s/__init__.py"%head): open("%s/__init__.py"%head, "w").write("# !/usr/bin/python\n# coding: utf-8")
        if not os.path.exists("%s/config.py"%head): os.system("cp config.py %s/"%head)
        return mrcfiles_initial, head

def load_images(mrcfiles_initial, head):
        # read first LOADNUM  images
        LOADNUM = params["loadnum"] # Number of micrographs loaded in one time
        if len(mrcfiles_initial) <= LOADNUM:
                args = readMrcFiles(mrcfiles_initial)
                del mrcfiles_initial[:]
        else:
                print "Loading first %d files to save memory ..."%LOADNUM
                args = readMrcFiles(mrcfiles_initial[:LOADNUM])
                del mrcfiles_initial[:LOADNUM]
        return mrcfiles_initial, args

if __name__ == "__main__":
        params = setupParserOptions()
        mrcfiles, head = initialize(params)
        os.chdir(head)
        from pfpicker_interactive import *
        from pfpicker_utils import *
        mrcfiles_initial, args = load_images(mrcfiles, head)
        window = pfpicker_interactive(params, mrcfiles_initial, *args)
        gtk.main()

