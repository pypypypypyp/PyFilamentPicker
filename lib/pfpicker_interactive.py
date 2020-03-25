#!/usr/bin/env python
#coding: utf-8

VERSION = "3.0.20.1"
POST = "Beta"

DESCRIPTION = """
                PyFilamentPicker ver.%s (%s)

PyFilamentPicker is a semi-automated filament tracer for 
cryo-EM particle picking. If you find this software useful,
please cite xxx.

(C) 2019 Yuta Komori


"""%(VERSION, POST)


import warnings

with warnings.catch_warnings():
        warnings.simplefilter("ignore")
warnings.simplefilter("ignore")

from pfpicker_utils import *
import time
import sys
import numpy as np
from EMAN2 import *
from math import *
from scipy.interpolate import spline, splrep, splprep, splev
from scipy.ndimage import map_coordinates as interpolate2d
from scipy.ndimage import zoom, rotate, interpolation
from scipy.ndimage.filters import gaussian_laplace as log
from scipy.ndimage.filters import gaussian_filter as gauss
from scipy.signal import *
from progressbar import ProgressBar
from StarRW import *
from matplotlib import pyplot as plt
import pygtk
pygtk.require("2.0")
import gtk
import cairo
import multiprocessing

global STARPARAMS
STARPARAMS =["_rlnMicrographName", "_rlnCoordinateX", "_rlnCoordinateY", "_rlnHelicalTubeID", "_rlnAngleTiltPrior", "_rlnAnglePsiPrior", "_rlnHelicalTrackLength", "_rlnAnglePsiFlipRatio"] # Header for output star file

def writeConfig():
        content = "APIX = %s\nHEIGHT = %s\nANGSTEP = %s\nDIST = %s\nSEGSTEP = %s\nANGRANGE = %s\nSTACK_BOXSIZE = %s\nFILAMENT_DIAMETER = %s\nINITIAL_ZOOM = %s\nINITIAL_THRES = %s\nSHRINK = %s\nLOCAL_FILAMENT_LENGTH = %s\nBIN = %s\nPERCENT_CONTRAST = %s\nSAVE = %s"
        open("config.py", "w").write(content%(str(APIX), str(HEIGHT), str(ANGSTEP), str(DIST), str(SEGSTEP), str(ANGRANGE), str(STACK_BOXSIZE), str(FILAMENT_DIAMETER), str(INITIAL_ZOOM), str(INITIAL_THRES), str(SHRINK), str(LOCAL_FILAMENT_LENGTH), str(BIN), str(PERCENT_CONTRAST), str(SAVE)))

class pfpicker_interactive:
        def __init__(self, params, mrcfiles):
                # micrographs and file names
                self.mrcfiles = mrcfiles
                
                # FDAGK
                rho = params["rho"]
                sigma = params["sigma"] / APIX / BIN
                self.FDAGK0 = FDAGK(sigma, 0, rho)
                
                # instance variables
                self.params = params
                if self.params["safemode"]:
                        self.imgindex = int(open("safemode.dat", "r").read())
                elif self.params["continue"]:
                        self.imgindex = int(open("continue.dat", "r").read())
                elif self.params["continuefrom"] != 1:
                        self.imgindex = self.params["continuefrom"]-1
                else:
                        self.imgindex = 0
                self.origimg, self.binimg, self.dirname, self.emname = readMrcFile(self.mrcfiles, self.imgindex)
                self.img = zoom(self.binimg, INITIAL_ZOOM/float(MAX_ZOOM), order=1) # img now showing
                self.imgbu = np.copy(self.img) # back up of img
                self.binheight, self.binwidth = self.binimg.shape
                self.currentheight, self.currentwidth = h, w = self.img.shape
                self.black = 0
                self.white = 255
                self.scale = float(INITIAL_ZOOM)
                self.pickthreshold = INITIAL_THRES
                self.all_positions = []
                self.all_angles = []
                self.precise_positions = []
                self.shrinks = []
                self.filament_numbers = []
                self.filament_types = []
                self.start_stop = []
                self.displayables = []
                self.savethreads = []
                self.needtomovescale = False
                self.active_filament = None
                self.toggled = False
                self.accept_keyevent = True
                self.shrink_clicks = 0
                self.active_filament_type = 0

                #load widgets from glade file
                builder = gtk.Builder()
                builder.add_from_file("%s/pfpicker.glade"%os.path.dirname(__file__))
                self.window = builder.get_object("window1")
                self.scrolledwindow = builder.get_object("scrolledwindow1")
                self.drarea = builder.get_object("drawingarea1")
                self.shrink = builder.get_object("image1")
                self.active_filamentcolor = builder.get_object("image2")
                self.label1 = builder.get_object("label1")
                self.label2 = builder.get_object("label2")
                self.label3 = builder.get_object("label3")
                self.button1 = builder.get_object("button1")
                self.button2 = builder.get_object("button2")
                self.button3 = builder.get_object("button3")
                self.button4 = builder.get_object("button4")
                self.button5 = builder.get_object("button5")
                self.prevbutton = builder.get_object("prevbutton")
                self.nextbutton = builder.get_object("nextbutton")
                self.togglebutton1 = builder.get_object("togglebutton1")
                self.combobox = builder.get_object("combobox1")
                self.comboboxentry = builder.get_object("comboboxentry1")
                self.eventbox = builder.get_object("eventbox1")
                self.vscale1 = builder.get_object("vscale1")
                self.vscale2 = builder.get_object("vscale2")
                self.vscale3 = builder.get_object("vscale3")
                self.vscale4 = builder.get_object("vscale4")
                self.saveprev = builder.get_object("save1")
                self.savenext = builder.get_object("save2")
                self.textview1 = builder.get_object("textview1")
                self.scrolledwindow2 = builder.get_object("scrolledwindow2")
                self.filemenu = builder.get_object("menu3")
                self.configmenu = builder.get_object("menu1")
                self.helpmenu = builder.get_object("menu2")
                self.liststore1 = gtk.ListStore(str)
                self.liststore2 = gtk.ListStore(str)
                self.embuf = gtk.gdk.pixbuf_new_from_array(self.img.repeat(3).reshape(h, w, 3), colorspace=gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                self.shblank = np.zeros((400, 200, 3)).astype(np.uint8)
                self.shrinkbuf = gtk.gdk.pixbuf_new_from_array(self.shblank, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample = 8)
                self.window.set_title("%s (%d/%d) -- PyFilamentPicker ver.%s (%s) --"%(self.emname, self.imgindex+1, len(self.mrcfiles), VERSION, POST))
                self.shrink.set_from_pixbuf(self.shrinkbuf)
                self.drarea.set_size_request(self.img.shape[1], self.img.shape[0])
                self.label1.set_text("Contrast")
                self.label2.set_text("Scale")
                self.label3.set_text("Threshold\n(R vallue)")
                self.button1.set_label("SAVE")
                self.button2.set_label("REFRESH")
                self.button3.set_label("REMOVE")
                self.button4.set_label("ADD")
                self.button5.set_label("DEL")
                self.togglebutton1.set_label("SET RANGE")
                self.vscale1.set_digits(0)
                self.vscale1.set_range(0, 255)
                self.vscale1.set_value(0)
                self.vscale2.set_digits(0)
                self.vscale2.set_range(0, 255)
                self.vscale2.set_value(255)
                self.vscale3.set_digits(0)
                self.vscale3.set_range(10, MAX_ZOOM)
                self.vscale3.set_value(self.scale)
                self.vscale4.set_range(0, 1.0)
                self.vscale4.set_value(self.pickthreshold)
                self.vscale1.set_update_policy(gtk.UPDATE_DISCONTINUOUS)
                self.vscale2.set_update_policy(gtk.UPDATE_DISCONTINUOUS)
                self.vscale3.set_update_policy(gtk.UPDATE_DISCONTINUOUS)
                self.vscale4.set_update_policy(gtk.UPDATE_DISCONTINUOUS)
                self.scrolledwindow2.set_policy(gtk.POLICY_NEVER, gtk.POLICY_ALWAYS)
                self.textview1.set_wrap_mode(gtk.WRAP_CHAR)
                self.textbuffer = self.textview1.get_buffer()
                self.combobox.set_model(self.liststore1)
                self.comboboxentry.set_model(self.liststore2)
                cell = gtk.CellRendererText()
                self.combobox.pack_start(cell, True)
                self.combobox.add_attribute(cell, "text", 0)
                #cell = gtk.CellRendererPixbuf()
                #self.combobox.pack_start(cell, True)
                #self.combobox.add_attribute(cell, "pixbuf", 1)
                self.liststore2.append(["type1"])
                self.comboboxentry.set_text_column(0)
                self.comboboxentry.set_active(0)

                # menu
                self.loaditem = gtk.ImageMenuItem("Load PFP file")
                self.saveitem = gtk.ImageMenuItem("Save PFP file")
                self.configitem = gtk.ImageMenuItem("General Config")
                self.saveconfitem = gtk.ImageMenuItem("Save Config")
                self.dispconfitem = gtk.ImageMenuItem("Display Config")
                self.advopitem = gtk.ImageMenuItem("Advanced Option")
                self.descriptionitem = gtk.ImageMenuItem("Description")
                self.helpitem = gtk.ImageMenuItem("Help")
                self.filemenu.add(self.loaditem)
                self.filemenu.add(self.saveitem)
                self.configmenu.add(self.configitem)
                self.configmenu.add(self.saveconfitem)
                self.configmenu.add(self.dispconfitem)
                self.configmenu.add(self.advopitem)
                self.helpmenu.add(self.descriptionitem)
                self.helpmenu.add(self.helpitem)
                
                # connect signals to methods
                self.window.connect("delete-event", self.ExitApp)
                self.drarea.connect("expose-event", self.OnDraw)
                self.vscale1.connect("value-changed", self.OnValue1Changed)
                self.vscale2.connect("value-changed", self.OnValue2Changed)
                self.vscale3.connect("value-changed", self.OnValue3Changed)
                self.vscale4.connect("value-changed", self.OnValue4Changed)
                self.button1.connect("button-press-event", self.OnButton1Pressed)
                self.button2.connect("button-press-event", self.OnButton2Pressed)
                self.button3.connect("button-press-event", self.OnButton3Pressed)
                self.button4.connect("button-press-event", self.OnButton4Pressed)
                self.button5.connect("button-press-event", self.OnButton5Pressed)
                self.prevbutton.connect("button-press-event", self.OnPrevbuttonPressed)
                self.nextbutton.connect("button-press-event", self.OnNextbuttonPressed)
                self.togglebutton1.connect("toggled", self.OnToggleButton1Toggled)
                self.scrolledwindow.connect("button-press-event", self.OnImageClicked)
                self.combobox.connect("changed", self.OnComboboxChanged)
                self.comboboxentry.child.connect("changed", self.OnComboboxEntryChanged)
                self.eventbox.connect("button-press-event", self.OnEventBoxClicked)
                self.loaditem.connect("activate", self.OnLoaditemActivated)
                self.saveitem.connect("activate", self.OnSaveitemActivated)
                self.configitem.connect("activate", self.OnConfigitemActivated)
                self.saveconfitem.connect("activate", self.OnSaveconfitemActivated)
                self.dispconfitem.connect("activate", self.OnDispconfitemActivated)
                self.advopitem.connect("activate", self.OnAdvopitemActivated)
                self.descriptionitem.connect("activate", self.OnDescriptionitemActivated)
                self.helpitem.connect("activate", self.OnHelpitemActivated)
                
                # connect key events to methods
                self.accel = gtk.AccelGroup()
                self.window.add_accel_group(self.accel)
                self.AddAccelerator(self.accel, "<Control>w", self.OnCtrlWPressed)
                self.AddAccelerator(self.accel, "<Control>q", self.OnCtrlQPressed)
                self.AddAccelerator(self.accel, "<Control>s", self.OnCtrlSPressed)
                self.AddAccelerator(self.accel, "<Control>a", self.OnCtrlAPressed)
                self.AddAccelerator(self.accel, "<Control>r", self.OnButton3Pressed)
                
                # display message
                self.DisplayMessage("Opened %s"%self.emname)
                self.DisplayMessage("Pixel size: %f"%APIX)
                
                # show all the widgets
                self.window.show_all()

        def AddAccelerator(self, accel, accelerator, callback):
                if accelerator is not None:
                        key, mod = gtk.accelerator_parse(accelerator)
                accel.connect_group(key, mod, 0, callback)

        def ChangeMicrograph(self):
                self.LoadMicrograph()
                self.ClearAllInformationOfTracedFilaments() 
                self.window.set_title("%s (%d/%d) -- PyFilamentPicker ver.%s (%s) --"%(self.emname, self.imgindex+1, len(self.mrcfiles), VERSION, POST))
                self.img = zoom(self.binimg, self.scale/float(MAX_ZOOM), order=1) # img now showing
                self.imgbu = np.copy(self.img) # backup of img
                self.currentheight, self.currentwidth = h, w = self.img.shape
                self.embuf = gtk.gdk.pixbuf_new_from_array(self.img.repeat(3).reshape(h, w, 3), colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample = 8)
                self.shrinkbuf = gtk.gdk.pixbuf_new_from_array(self.shblank, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample = 8)
                self.MaskAndChangeScale()

        def LoadMicrograph(self):
                origimg, binimg, dirname, emname = readMrcFile(self.mrcfiles, self.imgindex)
                self.origimg = origimg
                self.binimg = binimg
                self.emname = emname
                self.dirname = dirname
        
        def OnPrevbuttonPressed(self, window, event):
                if self.saveprev.get_active():
                        self.SaveFilamentsAndChangeMicrograph(True, -1)
                else:
                        self.SaveFilamentsAndChangeMicrograph(False, -1)

        def OnNextbuttonPressed(self, window, event):
                if self.saveprev.get_active():
                        self.SaveFilamentsAndChangeMicrograph(True, 1)
                else:
                        self.SaveFilamentsAndChangeMicrograph(False, 1)
        
        def OnCtrlSPressed(self, *args):
                self.SaveFilamentsAndChangeMicrograph(False, 1)

        def OnCtrlWPressed(self, *args):
                self.SaveFilamentsAndChangeMicrograph(True, 1)

        def OnCtrlAPressed(self, *args):
                self.SaveFilamentsAndChangeMicrograph(False, -1)

        def OnCtrlQPressed(self, *args):
                self.SaveFilamentsAndChangeMicrograph(True, -1)

        def SaveFilamentsAndChangeMicrograph(self, save, move):
                if not self.accept_keyevent: return
                if move == 1:
                        if self.imgindex == len(self.mrcfiles)-1: return
                elif move == -1:
                        if self.imgindex == 0: return
                if save: self.SaveAllFilaments()
                self.SavePFPFile()
                self.imgindex += move
                self.ChangeMicrograph()

        def SaveAllFilaments(self):
                self.DisplayMessage("Saving. Please don't close the window...")
                for i in range(len(self.all_positions)):
                        filament_type = self.liststore2[self.filament_types[i]][0]
                        self.savethreads.append(multiprocessing.Process(target=self.SaveThisFilament, args=(filament_type, self.precise_positions[i], self.displayables[i], self.imgindex, i+1)))
                        self.savethreads[-1].start()
                        self.DisplayMessage(" -> Working on filament %d (PID: %d)..."%(self.filament_numbers[i]+1, self.savethreads[-1].pid))

        def DisplayMessage(self, txt):
                if txt.startswith(" "):
                        txt = txt+"\n"
                elif not txt.startswith("\t"):
                        txt = "> "+txt+"\n"
                else:
                        txt = txt+"\n"
                self.textbuffer.insert_at_cursor(txt)
                while gtk.events_pending():
                        gtk.main_iteration()
                it = self.textbuffer.get_iter_at_line(self.textbuffer.get_line_count())
                self.textview1.scroll_to_iter(it, 0.)

        def ExitApp(self, window, event):
                # open("")
                self.DisplayMessage("Closing the window ...")
                self.SavePFPFile()
                self.SaveNewInFile()
                sys.exit(0)
        
        def SaveNewInFile(self):
                open("continue.dat", "w").write(str(self.imgindex))

        def SavePFPFile(self):
                cont = []
                for i in range(len(self.all_positions)):
                        cont.append(str(self.filament_numbers[i]))
                        cont.append(str(list(self.all_positions[i][0])))
                        cont.append(str(list(self.all_positions[i][1])))
                open(self.dirname+"/"+os.path.basename(self.emname[:-4])+".pfp", "w").write("\n".join(cont))
                # save safemode.in
                open("safemode.dat", "w").write(str(self.imgindex+1))

        def SaveThisFilament(self, filament_type, precise_positions, displayables, imgindex, tubeid):
                savedir = self.dirname
                emname = self.emname
                basename = os.path.basename(emname)[:-4]
                num = tubeid
                positions, angles, displayables = prepareOutput(precise_positions, displayables)
                positions_angles = np.array([np.array([positions[j],angles[j]]) for j in range(len(positions)) if displayables[j] and positions[j][1]-STACK_BOXSIZE/2 > 0 and positions[j][1]+STACK_BOXSIZE/2 < self.binheight * BIN and positions[j][0]-STACK_BOXSIZE/2 > 0 and positions[j][0]+STACK_BOXSIZE/2 < self.binwidth * BIN])
                positions = positions_angles[:, 0]
                angles = positions_angles[:, 1]
                filename = basename+"_stack%03d_%s"%(num, filament_type)
                filename_dir = savedir+"/"+filename
                remove = [savedir+"/"+i for i in os.listdir(savedir) if i.startswith(basename+"_stack%03d"%num)]
                for i in remove: os.remove(i)
                if SAVE[0] or self.params["savestack"]:
                        [EMNumPy.numpy2em(self.origimg[int(pos[1])-STACK_BOXSIZE/2:int(pos[1])+STACK_BOXSIZE/2, int(pos[0])-STACK_BOXSIZE/2:int(pos[0])+STACK_BOXSIZE/2].astype(np.float32)).write_image(filename_dir+".mrcs", -1) for pos in positions]
                if SAVE[1]:
                        header = create_header(STARPARAMS)
                        content = [[emname[3:], positions[j][0], positions[j][1], tubeid, 90, 90-angles[j], j*SEGSTEP/APIX, 0.5] for j in range(len(positions))]
                        write_star(filename_dir+".star", header, content)
                if SAVE[2]:
                        content = ["%d\t%d\t%d\t%d"%(positions[j][0]-STACK_BOXSIZE/2, positions[j][1]-STACK_BOXSIZE/2, STACK_BOXSIZE, STACK_BOXSIZE) for j in range(len(positions))]
                        open(filename_dir+".box", "w").write("\n".join(content))

        def OnLoaditemActivated(self, window):
                pfpname = self.dirname+"/"+os.path.basename(self.emname[:-4])+".pfp"
                if os.path.exists(pfpname):
                        self.DisplayMessage("Loading PFP file ...")
                        self.ClearAllInformationOfTracedFilaments()
                        pfp = open(pfpname, "r").readlines()
                        for i in range(len(pfp)/3):
                                filnum = int(pfp[i*3+0])
                                exec("positions = (np.array(%s), np.array(%s))"%(pfp[i*3+1], pfp[i*3+2]))
                                self.all_positions.append(positions)
                                n = positions[0].size
                                self.start_stop.append([0, None])
                                self.displayables.append([True]*n)
                                if self.active_filament is None: prev_active = None
                                else: prev_active = float(self.active_filament)
                                if len(self.filament_numbers) > 0: 
                                        self.active_filament = len(self.filament_numbers)
                                        self.filament_numbers.append(max(self.filament_numbers)+1)
                                else:
                                        self.active_filament = 0
                                        self.filament_numbers.append(0)
                                # draw on a screen
                                self.RedrawAllFilaments()
                                res = self.DrawShrink(positions, -1, True)
                                if not res:
                                        self.all_positions.pop()
                                        self.start_stop.pop()
                                        self.displayables.pop()
                                        self.filament_numbers.pop()
                                        if prev_active is None: self.active_filament = None
                                        else: self.active_filament = int(prev_active)
                                        self.RedrawAllFilaments()
                                        return
                                # add a row to combobox
                                #r, g, b = COLORS[self.filament_numbers[self.active_filament]%len(COLORS)]
                                #colorarray = np.zeros(10*10*3).reshape((10, 10, 3))
                                #colorarray[:, :, 0] = int(r*255); colorarray[:, :, 1] = int(g*255); colorarray[:, :, 2] = int(b*255);
                                #colorarray = colorarray.astype(np.uint8)
                                #colorbuf = gtk.gdk.pixbuf_new_from_array(colorarray, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                                self.liststore1.append(["Filament %d"%(max(self.filament_numbers)+1)])
                                self.combobox.set_model(self.liststore1)
                                self.combobox.set_active(self.active_filament)
                        self.DisplayMessage("Done.")
                else:
                        self.DisplayMessage("PFP file doesn't exist!")

        def OnSaveitemActivated(self, window):
                self.SavePFPFile()

        def OnConfigitemActivated(self, window):
                configwindow = ConfigWindow()

        def OnSaveconfitemActivated(self, window):
                saveconfwindow = SaveconfWindow()

        def OnDispconfitemActivated(self, window):
                dispconfwindow = DispconfWindow()

        def OnAdvopitemActivated(self, window):
                advoptwindow = AdvoptWindow()

        def OnDescriptionitemActivated(self, window):
                descriptionwindow = DescriptionWindow()

        def OnHelpitemActivated(self, window):
                pass

        def OnToggleButton1Toggled(self, window):
                self.shrink_clicks = 0

        def OnEventBoxClicked(self, window, event):
                if self.active_filament >= 0 and self.togglebutton1.get_active():
                        y = event.y
                        height = window.allocation.height
                        shrinkheight = len(self.shrinks[self.active_filament])-1
                        start_y = height/2.-shrinkheight/2.
                        newy = y-start_y
                        self.shrink_clicks += 1
                        if self.shrink_clicks%2 == 1:
                                if newy < 0: self.start_stop[self.active_filament][0] = 0
                                elif newy > shrinkheight: self.start_stop[self.active_filament][0] = None
                                else: self.start_stop[self.active_filament][0] = int(newy/shrinkheight*len(self.displayables[self.active_filament]))
                                self.displayables[self.active_filament] = [True]*len(self.displayables[self.active_filament])
                                if self.start_stop[self.active_filament][0] != None: self.displayables[self.active_filament][0:self.start_stop[self.active_filament][0]] = [False]*self.start_stop[self.active_filament][0]
                                else: self.displayables[self.active_filament] = [False]*len(self.displayables[self.active_filament])
                        elif self.shrink_clicks%2 == 0 and self.start_stop[self.active_filament][0] < int(newy/shrinkheight*len(self.displayables[self.active_filament])):
                                if newy < 0: self.start_stop[self.active_filament][1] = 0
                                elif newy > shrinkheight: self.start_stop[self.active_filament][1] = None
                                else: self.start_stop[self.active_filament][1] = int(newy/shrinkheight*len(self.displayables[self.active_filament]))
                                if self.start_stop[self.active_filament][1] != None: self.displayables[self.active_filament][self.start_stop[self.active_filament][1]:] = [False]*(len(self.displayables[self.active_filament])-self.start_stop[self.active_filament][1])
                                else: pass
                        if self.shrink_clicks%2 == 0: self.togglebutton1.set_active(False)
                        # on shrink
                        shrink = self.shrinks[self.active_filament]
                        shrink = self.SetStartStopLineOnShrinkImage(shrink.copy())
                        self.shrinkbuf = gtk.gdk.pixbuf_new_from_array(shrink, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                        self.shrink.set_from_pixbuf(self.shrinkbuf)
                        # on micrograph
                        self.RefreshDisplay()
                        self.RedrawAllFilaments()

        def ClearAllInformationOfTracedFilaments(self):
                self.togglebutton1.set_active(False)
                self.RefreshDisplay()
                self.all_positions = []
                self.shrinks = []
                self.filament_numbers = []
                self.filament_types = []
                self.start_stop = []
                self.displayables = []
                self.precise_positions = []
                self.shrinkbuf = gtk.gdk.pixbuf_new_from_array(self.shblank, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample = 8)
                self.shrink.set_from_pixbuf(self.shrinkbuf)
                self.combobox.get_model().clear()

        def OnButton1Pressed(self, window, event):
                self.DisplayMessage("Saving. Please don't close the window...")
                for i in range(len(self.all_positions)):
                        filament_type = self.liststore2[self.filament_types[i]][0]
                        self.savethreads.append(multiprocessing.Process(target=self.SaveThisFilament, args=(filament_type, self.precise_positions[i], self.displayables[i], self.imgindex, i+1)))
                        self.savethreads[-1].start()
                        self.DisplayMessage(" -> Working on filament %d (PID: %d)..."%(self.filament_numbers[i]+1, self.savethreads[-1].pid))

        def OnButton2Pressed(self, window, event):
                self.ClearAllInformationOfTracedFilaments()

        def OnButton3Pressed(self, window, *args):
                if len(self.all_positions) > 0:
                        self.togglebutton1.set_active(False)
                        nowactive = int(self.active_filament)
                        del self.all_positions[self.active_filament]
                        del self.shrinks[self.active_filament]
                        del self.filament_numbers[self.active_filament]
                        del self.filament_types[self.active_filament]
                        del self.start_stop[self.active_filament]
                        del self.displayables[self.active_filament]
                        del self.precise_positions[self.active_filament]
                        self.liststore1.clear()
                        for i in range(len(self.filament_numbers)):
                                #r, g, b = COLORS[self.filament_numbers[i]%len(COLORS)]
                                #colorarray = np.zeros(10*10*3).reshape((10, 10, 3))
                                #colorarray[:, :, 0] = int(r*255); colorarray[:, :, 1] = int(g*255); colorarray[:, :, 2] = int(b*255);
                                #colorarray = colorarray.astype(np.uint8)
                                #colorbuf = gtk.gdk.pixbuf_new_from_array(colorarray, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                                self.liststore1.append(["Filament %d"%(self.filament_numbers[i]+1)])
                        self.combobox.set_model(self.liststore1)
                        self.active_filament = nowactive
                        if self.active_filament == len(self.all_positions): self.active_filament -= 1
                        if self.active_filament >= 0:
                                self.combobox.set_active(self.active_filament)
                                shrink = self.shrinks[self.active_filament]
                                shrink = self.SetStartStopLineOnShrinkImage(shrink.copy())
                                self.shrinkbuf = gtk.gdk.pixbuf_new_from_array(shrink, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                        else:
                                self.shrinkbuf = gtk.gdk.pixbuf_new_from_array(self.shblank, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                        if self.active_filament >= 0: self.combobox.set_active(self.active_filament)
                        self.shrink.set_from_pixbuf(self.shrinkbuf)
                        self.RefreshDisplay()
                        self.RedrawAllFilaments()

        def OnButton4Pressed(self, window, *args):
                col = len(self.liststore2)+1
                self.liststore2.append(["type%d"%col])
                self.comboboxentry.set_model(self.liststore2)
                self.comboboxentry.set_active(col-1)

        def OnButton5Pressed(self, window, *args):
                if len(self.liststore2) == 1: return
                d = 0
                if self.active_filament_type == len(self.liststore2)-1:
                        d = -1
                iter = self.liststore2.get_iter(self.active_filament_type)
                self.liststore2.remove(iter)
                self.active_filament_type += d
                self.comboboxentry.set_active(self.active_filament_type)

        def OnComboboxChanged(self, window):
                self.active_filament = self.combobox.get_active()
                if self.active_filament == -1: return
                elif self.active_filament >= 0:
                        shrink = self.shrinks[self.active_filament]
                        shrink = self.SetStartStopLineOnShrinkImage(shrink.copy())
                        self.shrinkbuf = gtk.gdk.pixbuf_new_from_array(shrink, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                        self.shrink.set_from_pixbuf(self.shrinkbuf)
                        self.RefreshDisplay()
                        self.RedrawAllFilaments()
                if len(self.filament_types) > 0: self.comboboxentry.set_active(self.filament_types[self.active_filament])

        def OnComboboxEntryChanged(self, entry):
                if self.comboboxentry.get_active() >= 0:
                        self.active_filament_type = self.comboboxentry.get_active()
                iter = self.liststore2.get_iter(self.active_filament_type)
                self.liststore2.set(iter, 0, entry.get_text())
                self.comboboxentry.set_model(self.liststore2)
                if self.active_filament is not None: self.filament_types[self.active_filament] = self.active_filament_type

        def SetStartStopLineOnShrinkImage(self, shrink):
                r, g, b = (1, 0, 0)
                colarray = np.array([int(r*255), int(g*255), int(b*255)]*len(shrink[0])).reshape(len(shrink[0]), 3)
                shrinkheight = len(shrink)
                if self.start_stop[self.active_filament][0] != None:
                        shrink[int(float(self.start_stop[self.active_filament][0])/len(self.displayables[self.active_filament])*shrinkheight), :] = colarray
                else: shrink[-1, :] = colarray
                if self.start_stop[self.active_filament][1] != None:
                        tmp = int(float(self.start_stop[self.active_filament][1])/len(self.displayables[self.active_filament])*shrinkheight)
                        shrink[tmp, :] = colarray
                else: shrink[-1, :] = colarray
                return shrink
        
        def RefreshDisplay(self):
                cr = self.drarea.window.cairo_create()
                cr.set_source_pixbuf(self.embuf, 0, 0)
                cr.paint()

        def RedrawAllFilaments(self):
                if len(self.all_positions) > 0:
                        cr = self.drarea.window.cairo_create()
                        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.img.shape[1], self.img.shape[0])
                        context = cairo.Context(surface)
                        context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
                        context.set_font_size(28)
                        context.set_line_width(3)
                        for i in range(len(self.all_positions)):
                                initial = True
                                if self.filament_numbers[i] == self.filament_numbers[self.active_filament]:
                                        r, g, b = (1, 0, 0)
                                else: r, g, b = (1, 1, 1)
                                x, y = self.all_positions[i]
                                context.set_source_rgba(r, g, b, 0.9)
                                spos = (int(x[0]*self.scale/100.*BIN), int(y[0]*self.scale/100.*BIN))
                                intv = SPLINE_ACCURACY * 20 # place points at the interval of 5 pixels
                                for n, j in enumerate(range(0, x.size-intv, intv)):
                                        spos = (int(x[j]*self.scale/100.*BIN), int(y[j]*self.scale/100.*BIN))
                                        npos = (int(x[j+intv]*self.scale/100.*BIN), int(y[j+intv]*self.scale/100.*BIN))
                                        if self.displayables[i][j] and n%2==0:
                                                if initial:
                                                        context.move_to(spos[0], spos[1])
                                                        context.show_text(str(self.filament_numbers[i]+1))
                                                        initial = False
                                                context.move_to(spos[0], spos[1])
                                                context.line_to(npos[0], npos[1]) 
                                context.stroke()
                        cr.set_source_surface(context.get_target())
                        cr.paint()
        
        def OnDraw(self, window, event):
                cr = self.drarea.window.cairo_create()
                cr.set_source_pixbuf(self.embuf, 0, 0)
                cr.paint()
                self.RedrawAllFilaments()

        def OnImageRightClicked(self, window, event):
                self.OnButton3Pressed(window)

        def OnImageLeftClicked(self, window, event):
                self.togglebutton1.set_active(False) # in case togglebutton1 is pressed before self.shrink_clicks counts twice
                vpos = window.get_vadjustment().get_value()
                y = event.y + vpos
                hpos = window.get_hadjustment().get_value()
                x = event.x + hpos
                self.DisplayMessage("Searching from position (%d, %d)."%(x, y))
                origx, origy = int(x*(100./self.scale)), int(y*(100./self.scale))
                positions = self.FindFilament(origx, origy)
                if len(positions) == 0: 
                        self.DisplayMessage("Error: This filament is too short!")
                        return
                if positions[0][1] > positions[-1][1]: positions.reverse()
                cr = self.drarea.window.cairo_create()
                cr.set_source_rgb(1, 1, 1)
                if self.params["debug"]:
                        for j in range(len(positions)):
                                spos = (int(positions[j][0]*self.scale/100.*BIN), int(positions[j][1]*self.scale/100.*BIN))
                                cr.arc(spos[0], spos[1], 2, 0, 2.0 * pi)
                        cr.stroke()
                        return
                positions = createSplineCurveForAllLength(positions) # given as tuple (np.array([x0, x1, ...]), np.array([y0, y1, ...]))
                self.DisplayMessage("\tCalculating spline curve ...")
                if positions is None: 
                        self.DisplayMessage("Error: This filament is too short!")
                        return
                n = positions[0].size
                if self.active_filament is None: prev_active = None
                else: prev_active = float(self.active_filament)
                if len(self.filament_numbers) > 0: 
                        self.active_filament = len(self.filament_numbers)
                        self.filament_numbers.append(max(self.filament_numbers)+1)
                else:
                        self.active_filament = 0
                        self.filament_numbers.append(0)
                self.filament_types.append(self.active_filament_type)
                # draw on a screen
                self.all_positions.append(positions)
                self.start_stop.append([0, None])
                self.displayables.append([True]*n)
                res = self.DrawShrink(positions, -1)
                # as filaments are drawn in self.OnComboboxChanged method, nothing about drawing is done here
                if res is None:
                        self.all_positions.pop()
                        self.start_stop.pop()
                        self.displayables.pop()
                        self.filament_numbers.pop()
                        self.filament_types.pop()
                        if prev_active is None: self.active_filament = None
                        else: self.active_filament = int(prev_active)
                        return
                # add a row to combobox
                #r, g, b = COLORS[self.filament_numbers[self.active_filament]%len(COLORS)]
                #colorarray = np.zeros(10*10*3).reshape((10, 10, 3))
                #colorarray[:, :, 0] = int(r*255); colorarray[:, :, 1] = int(g*255); colorarray[:, :, 2] = int(b*255);
                #colorarray = colorarray.astype(np.uint8)
                #colorbuf = gtk.gdk.pixbuf_new_from_array(colorarray, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                self.liststore1.append(["Filament %d"%(max(self.filament_numbers)+1)])
                self.combobox.set_model(self.liststore1)
                self.combobox.set_active(self.active_filament)

        def OnImageMiddleClicked(self, window, event):
                self.OnCtrlNPressed()

        def OnImageClicked(self, window, event):
                if event.button == 1: self.OnImageLeftClicked(window, event)
                elif event.button == 2: self.OnImageMiddleClicked(window, event)
                elif event.button == 3: self.OnImageRightClicked(window, event)

        def DrawShrink(self, positions, index, loading=False):
                if not loading: self.DisplayMessage("\tCreating shrink image...")
                shrink_orig = straighten(positions, self.binimg)
                if shrink_orig is None:
                        self.DisplayMessage("Error: An internal error occured creating the shrink image. Filament may be too short! Please try again.")
                        return None
                if self.params["debug"]: save(shrink_orig, "shrink_orig.mrc")
                # make shrink image easy-to-interpret
                shrink_orig = flatten(shrink_orig)
                shrink_orig = cutoffNoise(shrink_orig, 1.2) # probably value around 1.2 is appropriate ...
                shrink_orig = toUINT8(shrink_orig)
                shrink = np.zeros((shrink_orig.shape[0]/SHRINK, shrink_orig.shape[1]))
                for i in range(shrink_orig.shape[0]/SHRINK):
                        shrink[i] = np.sum(shrink_orig[i*SHRINK:(i+1)*SHRINK], axis=0)
                shrink = toUINT8(shrink)
                shrink = shrink.repeat(3).reshape(shrink.shape[0], shrink.shape[1], 3)
                if index == -1: self.shrinks.append(shrink.copy())
                else: self.shrinks[index] = shrink.copy()
                shrink_disp = self.SetStartStopLineOnShrinkImage(shrink)
                self.shrinkbuf = gtk.gdk.pixbuf_new_from_array(shrink_disp, colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample=8)
                self.shrink.set_from_pixbuf(self.shrinkbuf)
                if not loading: self.DisplayMessage("Done.")
                self.precise_positions.append([(positions[0][i]*BIN, positions[1][i]*BIN) for i in range(positions[0].size)]) # precise_positions is not binned
                return True

        def FindFilament(self, origx, origy):
                # initial angular search
                r = 1.0
                width = int(FILAMENT_DIAMETER/BIN/APIX*r + self.FDAGK0.shape[0]) + 1
                x, y = np.meshgrid(np.arange(-width/2, width/2), np.arange(-HEIGHT_PIX/BIN/2, HEIGHT_PIX/BIN/2))
                sums = [None]*(180/ANGSTEP)
                #ds = [None] * (180/ANGSTEP)
                mins = np.arange(180./ANGSTEP)
                maxs = np.arange(180./ANGSTEP)
                for i in range(0, 180, ANGSTEP):
                        xx = origx/BIN + x*np.cos(np.deg2rad(i)) - y*np.sin(np.deg2rad(i))
                        yy = origy/BIN + x*np.sin(np.deg2rad(i)) + y*np.cos(np.deg2rad(i))
                        rot = interpolate2d(self.binimg, [yy, xx], order=0).reshape(xx.shape)
                        rot = (255 - rot) / 255. # invert contrast
                        rot_filtered = fftconvolve(rot, self.FDAGK0, mode="valid")
                        if self.params["debug"]:
                                EMNumPy.numpy2em(rot_filtered).write_image("rot_filtered.mrcs", -1)
                                EMNumPy.numpy2em(rot).write_image("rot.mrcs", -1)
                        size = rot_filtered.shape[0]
                        isum = np.sum(rot_filtered, axis=0) / float(size)
                        mins[i/ANGSTEP] = np.min(isum)
                        maxs[i/ANGSTEP] = np.max(isum)
                        #ds[i/ANGSTEP] = isum.std()
                        sums[i/ANGSTEP] = isum
                ds = maxs-mins
                if self.params["debug"]: print ds
                ang_index = np.where(ds == np.max(ds))[0][0]
                if ang_index%2 == 0: ang_index_norm = ang_index
                else:
                        cands = [ang_index+1, ang_index-1]
                        cands = [i-i/(180/ANGSTEP) for i in cands]
                        if ds[cands[0]] > ds[cands[1]]: ang_index_norm = cands[0]
                        else: ang_index_norm = cands[1]
                angle = ang_index*ANGSTEP
                angle_norm = ang_index_norm*ANGSTEP
                sum_best = sums[ang_index] # unfiltered image
                self.DisplayMessage("\tThe angle of this particle is %d degrees." % angle)
                # two unit vectors for initial search direction
                uvector1 = (np.cos((angle-90)/180.*np.pi), np.sin((angle-90)/180.*np.pi))
                uvector1 = np.array([uvector1[0]/sqrt(uvector1[0]**2+uvector1[1]**2), uvector1[1]/sqrt(uvector1[0]**2+uvector1[1]**2)])
                uvector2 = np.array([-uvector1[0], -uvector1[1]])
                # calculate shift
                limit = FILAMENT_DIAMETER/BIN/APIX*r
                shift, corval = findShiftToCenter(sum_best, limit=limit, ref=None) #find shift using unfiltered image
                size = sum_best.size
                radius = floor(FILAMENT_DIAMETER/APIX/BIN/2)
                xnew = np.arange(-radius, radius+1)
                xnew = (size-1)/2. + shift + xnew
                self.ref = np.interp(xnew, np.arange(size), sum_best)
                self.corval_initial = np.max(np.correlate(sum_best, self.ref, mode="valid"))
                # re-pick clicked particle at the exact center by larger box size (multiplied by sqrt(2))
                origx /= BIN; origy /= BIN
                origx -= np.cos(np.deg2rad(180-angle))*shift
                origy += np.sin(np.deg2rad(180-angle))*shift
                # it's time to search!
                # note that the clipped particle orientation is not the same as that of displayed one (vertically flipped).
                self.DisplayMessage("\tSearching...")
                positions1 = self.SearchToThisDirection((origx, origy), uvector1, angle_norm)
                positions2 = self.SearchToThisDirection((origx, origy), uvector2, angle_norm)
                try:
                        if positions1[-1][1] <= positions2[-1][1]:
                                positions = [positions1[-i-1] for i in range(len(positions1))] + [(int(origx), int(origy))] + positions2
                        else:
                                positions = [positions2[-i-1] for i in range(len(positions2))] + [(int(origx), int(origy))] + positions1
                except:
                        positions = [] # error should occur
                return positions

        def SearchToThisDirection(self, pos, vec, angle):
                npos = (int(pos[0]+DIST/APIX/BIN*vec[0]), int(pos[1]+DIST/APIX/BIN*vec[1]))
                positions = []
                margin = STACK_BOXSIZE/BIN/2
                while margin < npos[0] and npos[0] < self.binwidth-margin and margin < npos[1] and npos[1] < self.binheight-margin:
                        pos, vec, angle = self.FindBestPositionAndDirection(pos, vec, angle)
                        if vec[0] == 0 and vec[1] == 0: break
                        positions.append(pos)
                        npos = (int(pos[0]+DIST/APIX/BIN*vec[0]), int(pos[1]+DIST/APIX/BIN*vec[1]))
                        if len(positions) >= 3 and positions[-1] == positions[-3]: break
                return positions

        def FindBestPositionAndDirection(self, pos, vec, angle):
                # initialize
                npos = [int(pos[0]+DIST/APIX/BIN*vec[0]), int(pos[1]+DIST/APIX/BIN*vec[1])]
                mins = [None] * (ANGRANGE*2+1)
                maxs = [None] * (ANGRANGE*2+1)
                #ds = [None] * (ANGRANGE*2+1)
                # crop the area to be searched from the microgrpah
                angles = [angle + i*ANGSTEP for i in range(-ANGRANGE, ANGRANGE+1)]
                width = int(FILAMENT_DIAMETER/BIN/APIX + DIST/APIX/BIN*np.sin(np.deg2rad(ANGRANGE*ANGSTEP))*2 + self.FDAGK0.shape[1]) + 1
                x, y = np.meshgrid(np.arange(-width/2, width/2), np.arange(-int(HEIGHT_PIX/BIN/2), int(HEIGHT_PIX/BIN/2)))
                for i, angle in enumerate(angles):
                        xx = npos[0] + x*np.cos(np.deg2rad(angle)) - y*np.sin(np.deg2rad(angle))
                        yy = npos[1] + x*np.sin(np.deg2rad(angle)) + y*np.cos(np.deg2rad(angle))
                        ptcl = interpolate2d(self.binimg, [yy, xx], order=0).reshape(xx.shape)
                        #if self.params["debug"]: EMNumPy.numpy2em(ptcl).write_image("ptcl.mrcs", -1)
                        ptcl = (255 - ptcl) / 255. # invert contrast
                        # apply FDAGK to the particle to determine the angle
                        ptcl_f = fftconvolve(ptcl, self.FDAGK0, mode="valid")
                        #if self.params["debug"]: EMNumPy.numpy2em(ptcl_f).write_image("%d_%d_f.mrcs"%(npos[0], npos[1]), -1)
                        #if self.params["debug"]: EMNumPy.numpy2em(ptcl).write_image("%d_%d_nf.mrcs"%(npos[0], npos[1]), -1)
                        isum = np.sum(ptcl_f, axis=0) / float(ptcl_f.shape[0])
                        #ds[i] = isum.std()
                        mins[i] = np.min(isum)
                        maxs[i] = np.max(isum)
                mins = np.array(mins)
                maxs = np.array(maxs)
                ds = maxs-mins
                # stop processing?
                # determine the best angle and re-extract the particle with the angle
                #ds_expanded = np.array([0] + list(ds) + [0])
                #peak_indices = argrelmax(ds_expanded)[0]-1
                #peak_indices_abs_min = np.min(np.abs(peak_indices-ANGRANGE))
                #candidates = [ANGRANGE-peak_indices_abs_min, ANGRANGE+peak_indices_abs_min]
                #if candidates[0] not in peak_indices: ang_index = candidates[1]
                #elif candidates[1] not in peak_indices: ang_index = candidates[0]
                #else:
                #        if ds[candidates[0]] > ds[candidates[1]]: ang_index = candidates[0]
                #        else: ang_index = candidates[1]
                ang_index = np.where(ds == np.max(ds))[0][0]
                best_angle = angles[ANGRANGE] + ANGSTEP * (ang_index - ANGRANGE)
                xx = npos[0] + x*np.cos(np.deg2rad(best_angle)) - y*np.sin(np.deg2rad(best_angle))
                yy = npos[1] + x*np.sin(np.deg2rad(best_angle)) + y*np.cos(np.deg2rad(best_angle))
                best_ptcl = interpolate2d(self.binimg, [yy, xx], order=0).reshape(xx.shape)
                best_ptcl = (255 - best_ptcl) / 255.
                #if out_of_mic >= 0: best_ptcl = best_ptcl[out_of_mic:, :]
                # find shift to center
                best_ptcl_filtered = fftconvolve(best_ptcl, self.FDAGK0, mode="valid")
                limit = ceil(DIST / APIX / BIN * np.sin(np.deg2rad(ANGSTEP*ANGRANGE))) * 2
                best_ptcl_filtered_for_sum = np.sum(best_ptcl_filtered, axis=0) / float(best_ptcl_filtered.shape[0])
                size = best_ptcl_filtered_for_sum.size
                #if self.params["debug"]:
                        #plt.plot(best_ptcl_filtered_for_sum)
                        #plt.savefig("%d_%d.png"%(npos[0], npos[1]))
                        #plt.clf()
                shift, corval = findShiftToCenter(best_ptcl_filtered_for_sum, limit=limit, ref=self.ref) # find shift 
                if corval < self.corval_initial * self.pickthreshold: return (0, 0), (0, 0), 0
                #if self.params["debug"]:
                        #plt.plot(self.ref)
                        #plt.savefig("%s.png"%time.time())
                        #plt.clf()
                # re-calculated new position
                npos[0] -= np.cos(np.deg2rad(180-best_angle))*shift
                npos[1] += np.sin(np.deg2rad(180-best_angle))*shift
                # filament angle and search vector
                true_angle = best_angle
                vec1 = np.array([np.cos(np.deg2rad(true_angle+90)), np.sin(np.deg2rad(true_angle+90))])
                vec2 = -1 * vec1
                #the vector which is close to the previous search vector is chosen as the next search vector
                if np.dot(vec1, vec) > np.dot(vec2, vec):
                        vec = vec1
                else: vec = vec2
                vec = vec/np.sqrt(np.sum(vec**2))
                return npos, vec, true_angle

        def OnValue1Changed(self, window):
                self.black = self.vscale1.get_value()
                self.MaskAndChangeScale()

        def OnValue2Changed(self, window):
                self.white =  self.vscale2.get_value()
                self.MaskAndChangeScale()

        def OnValue3Changed(self, window):
                self.scale = self.vscale3.get_value()
                self.needtomovescale = True
                self.MaskAndChangeScale()

        def OnValue4Changed(self, window):
                self.pickthreshold = self.vscale4.get_value()

        def MaskAndChangeScale(self):
                if self.needtomovescale:
                        self.img = zoom(self.binimg, self.scale/float(MAX_ZOOM), order=1)
                        self.needtomovescale = False
                        self.imgbu = np.copy(self.img)
                else:
                        self.img = np.copy(self.imgbu)
                high_values_indices = self.img > self.white
                low_values_indices = self.img < self.black
                self.img[low_values_indices] = self.black
                self.img[high_values_indices] = self.white
                self.img = toUINT8(self.img)
                self.Display()

        def Display(self):
                h, w = self.img.shape
                self.drarea.set_size_request(w, h)
                self.embuf = gtk.gdk.pixbuf_new_from_array(self.img.repeat(3).reshape(h, w, 3), colorspace = gtk.gdk.COLORSPACE_RGB, bits_per_sample = 8)
                cr = self.drarea.window.cairo_create()
                cr.set_source_pixbuf(self.embuf, 0, 0)
                cr.paint()
                self.RedrawAllFilaments()

class ConfigWindow:
        def __init__(self):
                builder = gtk.Builder()
                builder.add_from_file("%s/sub_windows.glade"%os.path.dirname(__file__))
                self.window = builder.get_object("window2")
                self.APIX = builder.get_object("entry17")
                self.APIX.set_text(str(APIX))
                self.HEIGHT = builder.get_object("entry1")
                self.HEIGHT.set_text(str(HEIGHT))
                self.DIST = builder.get_object("entry3")
                self.DIST.set_text(str(DIST))
                self.SEGSTEP = builder.get_object("entry4")
                self.SEGSTEP.set_text(str(SEGSTEP))
                self.STACK_BOXSIZE = builder.get_object("entry7")
                self.STACK_BOXSIZE.set_text(str(STACK_BOXSIZE))
                self.FILAMENT_DIAMETER = builder.get_object("entry9")
                self.FILAMENT_DIAMETER.set_text(str(FILAMENT_DIAMETER))
                self.ANGSTEP = builder.get_object("entry16")
                self.ANGSTEP.set_text(str(ANGSTEP))
                self.ANGRANGE = builder.get_object("entry18")
                self.ANGRANGE.set_text(str(ANGRANGE))
                self.button4 = builder.get_object("button4")
                self.button4.connect("button-press-event", self.OnSaveButtonPressed)
                self.button5 = builder.get_object("button5")
                self.button5.connect("button-press-event", self.OnCloseButtonPressed)
                self.window.show_all()

        def OnSaveButtonPressed(self, window, event): # Save button
                global APIX, HEIGHT, HEIGHT, DIST, SEGSTEP, STACK_BOXSIZE, FILAMENT_DIAMETER, ANGSTEP, ANGRANGE
                APIX = float(self.APIX.get_text())
                HEIGHT = int(self.HEIGHT.get_text())
                DIST = float(self.DIST.get_text())
                SEGSTEP = float(self.SEGSTEP.get_text())
                STACK_BOXSIZE = int(self.STACK_BOXSIZE.get_text())
                FILAMENT_DIAMETER = int(self.FILAMENT_DIAMETER.get_text())
                HEIGHT_PIX = int(HEIGHT / APIX) # pixelize
                if HEIGHT % 2 == 1: HEIGHT += 1 # HEIGHT should be an even number
                ANGSTEP = int(self.ANGSTEP.get_text())
                ANGRANGE = int(self.ANGRANGE.get_text())
                writeConfig()
                loadGlobalVariables()

        def OnCloseButtonPressed(self, window, event): # Close button
                self.window.destroy()

class SaveconfWindow:
        def __init__(self):
                builder = gtk.Builder()
                builder.add_from_file("%s/sub_windows.glade"%os.path.dirname(__file__))
                self.window = builder.get_object("window6")
                self.MRCS = builder.get_object("checkbutton_win6_1")
                self.STAR = builder.get_object("checkbutton_win6_2")
                self.BOX = builder.get_object("checkbutton_win6_3")
                if SAVE[0]: self.MRCS.set_active(True)
                if SAVE[1]: self.STAR.set_active(True)
                if SAVE[2]: self.BOX.set_active(True)
                self.window.show_all()
                self.button1 = builder.get_object("button_win6_1")
                self.button1.connect("button-press-event", self.OnSaveButtonPressed)
                self.button2 = builder.get_object("button_win6_2")
                self.button2.connect("button-press-event", self.OnCloseButtonPressed)

        def OnSaveButtonPressed(self, window, event):
                global SAVE
                SAVE[0] = self.MRCS.get_active()
                SAVE[1] = self.STAR.get_active()
                SAVE[2] = self.BOX.get_active()
                writeConfig()
                loadGlobalVariables()

        def OnCloseButtonPressed(self, window, event):
                self.window.destroy()
                

class DispconfWindow:
        def __init__(self):
                builder = gtk.Builder()
                builder.add_from_file("%s/sub_windows.glade"%os.path.dirname(__file__))
                self.window = builder.get_object("window3")
                self.INITIAL_ZOOM = builder.get_object("entry11")
                self.INITIAL_ZOOM.set_text(str(INITIAL_ZOOM))
                self.SHRINK = builder.get_object("entry13")
                self.SHRINK.set_text(str(SHRINK))
                self.PERCENT_CONTRAST = builder.get_object("entry14")
                self.PERCENT_CONTRAST.set_text(str(PERCENT_CONTRAST))
                self.BIN = builder.get_object("entry401")
                self.BIN.set_text(str(BIN))
                self.button6 = builder.get_object("button6")
                self.button6.connect("button-press-event", self.OnSaveButtonPressed)
                self.button7 = builder.get_object("button7")
                self.button7.connect("button-press-event", self.OnCloseButtonPressed)
                self.window.show_all()

        def OnSaveButtonPressed(self, window, event):
                global INITIAL_ZOOM, SHRINK, PERCENT_CONTRAST, BIN
                INITIAL_ZOOM = int(self.INITIAL_ZOOM.get_text())
                SHRINK = int(self.SHRINK.get_text())
                PERCENT_CONTRAST = float(self.PERCENT_CONTRAST.get_text())
                BIN = int(self.BIN.get_text())
                writeConfig()
                loadGlobalVariables()

        def OnCloseButtonPressed(self, window, event):
                self.window.destroy()

class AdvoptWindow:
        def __init__(self):
                builder = gtk.Builder()
                builder.add_from_file("%s/sub_windows.glade"%os.path.dirname(__file__))
                self.window = builder.get_object("window4")
                self.LOCAL_FILAMENT_LENGTH = builder.get_object("entry15")
                self.LOCAL_FILAMENT_LENGTH.set_text(str(LOCAL_FILAMENT_LENGTH))
                self.INITIAL_THRES = builder.get_object("entry405")
                self.INITIAL_THRES.set_text(str(INITIAL_THRES))
                self.button8 = builder.get_object("button8")
                self.button8.connect("button-press-event", self.OnSaveButtonPressed)
                self.button9 = builder.get_object("button9")
                self.button9.connect("button-press-event", self.OnCloseButtonPressed)
                self.window.show_all()

        def OnSaveButtonPressed(self, window, event):
                global LOCAL_FILAMENT_LENGTH, INITIAL_THRES
                LOCAL_FILAMENT_LENGTH = int(self.LOCAL_FILAMENT_LENGTH.get_text())
                INITIAL_THRES = float(self.INITIAL_THRES.get_text())
                writeConfig()
                loadGlobalVariables()

        def OnCloseButtonPressed(self, window, event):
                self.window.destroy()

class DescriptionWindow:
        def __init__(self):
                builder = gtk.Builder()
                builder.add_from_file("%s/sub_windows.glade"%os.path.dirname(__file__))
                self.window = builder.get_object("window5")
                self.label = builder.get_object("label46")
                self.button = builder.get_object("button10")
                self.button.connect("button-press-event", self.OnCloseButtonPressed)
                self.button.set_label("Close")
                self.label.set_text(DESCRIPTION)
                self.window.show_all()

        def OnCloseButtonPressed(self, window, event):
                self.window.destroy()

