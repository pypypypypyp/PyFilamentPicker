#!/usr/bin/env python
#coding: utf-8

DOCUMENT = """
**********************
        StarRW
**********************

Description:
        This module helps your reading and writing of RELION data.star files which contain information about particles.

Usage:
        >>> star = Star("test.star")
        >>> star.count_particles()
        2647
        >>> star.get_micrograph_names()
        ['Micrographs/001.mrc', 'Micrographs/002.mrc', 'Micrographs/003.mrc', 'Micrographs/004.mrc']
        >>> star.remove_this_micrograph("Micrographs/001.mrc")
        >>> star.get_micrograph_names()
        ['Micrographs/002.mrc', 'Micrographs/003.mrc', 'Micrographs/004.mrc']
        >>> star.write("test2.star")

Required environment:
        EMAN2 (for Python 2.7.x)
"""

import copy
import warnings
from EMAN2 import *

ALL = "micrograph_all"

def doc():
        print DOCUMENT

def isfloat(str):
        try:
                float(str)
                return True
        except ValueError:
                return False

class Star:
        def __init__(self, filename):
                content = [i.rstrip() for i in open(filename, "r").readlines()]
                if "data_optics" in content: # relion v.3.1 or above
                        optics_start = content.index("data_optics")
                        optics_end = content.index("data_micrographs")
                        optics_lines = [cont for cont in content[optics_start:optics_end] if cont != "" and not cont.startswith("#")]
                        optics_lines_header = [cont.split()[0] for cont in optics_lines if cont.startswith("_rln") or cont.startswith("data_") or cont.startswith("loop_")] 
                        optics_lines_content = [cont.split() for cont in optics_lines if not cont.startswith("_rln") and not cont.startswith("data_") and not cont.startswith("loop_")] 
                        micrographs_lines = content[optics_end:]
                else: 
                        optics_lines_header = []
                        optics_lines_content = []
                        micrographs_lines = content
                header_last_index = 0
                for i in range(len(micrographs_lines)):
                        if micrographs_lines[i].startswith("_rln"): header_last_index = i
                header = [cont for cont in micrographs_lines[:header_last_index+1] if cont != "" and not cont.startswith("#")]
                content = micrographs_lines[header_last_index+1:]
                for i in range(len(header)):
                        if header[i] != "": header[i] = header[i].split()[0]
                content = [cont for cont in content if cont != ""]
                for i in range(len(content)):
                        content[i] = content[i].split()
                        for j in range(len(content[i])):
                                if content[i][j].isdigit(): content[i][j] = int(content[i][j])
                                elif isfloat(content[i][j]): content[i][j] = float(content[i][j])
                self.header = header
                self.content = content

        def find_parameter_index(self, parameter_name):
                params = []
                for i in self.header:
                        if i.startswith("_rln"): params.append(i)
                index = params.index(parameter_name)
                return index

        def remove_this_parameter(self, parameter_index):
                newcontent = self.content[:]
                for i in range(len(self.content)):
                        del newcontent[i][parameter_index]
                self.content = newcontent
                rlnStart = 0
                for i in range(len(self.header)):
                        if self.header[i].startswith("_rln"):
                                rlnStart = i
                                break
                del self.header[parameter_index+rlnStart]

        def remove_this_micrograph(self,  micrograph_name):
                index = self.find_parameter_index("_rlnMicrographName")
                newcontent = [self.content[i] for i in range(len(self.content)) if self.content[i][index] != micrograph_name]
                self.content = newcontent

        def split_by_image_name(self):
                spl = [[]]
                index = self.find_parameter_index("_rlnImageName")
                current_image_name = self.content[0][index][6:]
                for cont in self.content:
                        if current_image_name != cont[index][6:]:
                                spl.append([])
                        spl[len(spl)-1].append(cont)
                        current_image_name = cont[index][6:]
                return spl

        def split_by_micrograph_name(self):
                spl = [[]]
                index = self.find_parameter_index("_rlnMicrographName")
                current_micrograph_name = self.content[0][index]
                for cont in self.content:
                        if current_micrograph_name != cont[index]:
                                spl.append([])
                        spl[len(spl)-1].append(cont)
                        current_micrograph_name = cont[index]
                return spl

        def sort_by_micrograph_name(self):
                index = self.find_parameter_index("_rlnMicrographName")
                self.content = sorted(self.content, key=lambda x:x[index])
        
        def sort_by_this_parameter(self, parameter_name):
                index = self.find_parameter_index(parameter_name)
                self.content = sorted(self.content, key=lambda x:x[index])

        def arrange_angle(self, angle):
                if angle >= 360: return self.arrange_angle(angle-360)
                elif angle < 0: return self.arrange_angle(angle+360)
                else: return angle

        def arrange_angles(self):
                rotindex = self.find_parameter_index("_rlnAngleRot")
                tiltindex = self.find_parameter_index("_rlnAngleTilt")
                psiindex = self.find_parameter_index("_rlnAnglePsi")
                for i in range(len(self.content)):
                        rot = self.content[i][rotindex]
                        tilt = self.content[i][tiltindex]
                        psi = self.content[i][psiindex]
                        self.content[i][rotindex] = self.arrange_angle(float(rot))
                        self.content[i][tiltindex] = self.arrange_angle(float(tilt))
                        self.content[i][psiindex] = self.arrange_angle(float(psi))

        def rotate(self, micrograph_name, drot, dtilt, dpsi):
                rotindex = self.find_parameter_index("_rlnAngleRot")
                tiltindex = self.find_parameter_index("_rlnAngleTilt")
                psiindex = self.find_parameter_index("_rlnAnglePsi")
                micrographnameindex = self.find_parameter_index("_rlnMicrographName")
                micrographnameindex = self.find_parameter_index("_rlnMicrographName")
                if micrograph_name == ALL:
                        for i in range(len(self.content)):
                                rot = self.content[i][rotindex]
                                tilt = self.content[i][tiltindex]
                                psi = self.content[i][psiindex]
                                orig = Transform({"type":"spider", "phi":rot, "theta":tilt, "psi":psi})
                                T = Transform({"type":"spider", "phi":drot, "theta":dtilt, "psi":dpsi})
                                new = T*orig
                                newparams = new.get_params("spider")
                                self.content[i][rotindex] = newparams["phi"]
                                self.content[i][tiltindex] = newparams["theta"]
                                self.content[i][psiindex] = newparams["psi"]
                else:
                        for i in range(len(self.content)):
                                if self.content[i][micrographnameindex] == micrograph_name:
                                        rot = self.content[i][rotindex]
                                        tilt = self.content[i][tiltindex]
                                        psi = self.content[i][psiindex]
                                        orig = Transform({"type":"spider", "phi":rot, "theta":tilt, "psi":psi})
                                        T = Transform({"type":"spider", "phi":drot, "theta":dtilt, "psi":dpsi})
                                        new = T*orig
                                        newparams = new.get_params("spider")
                                        self.content[i][rotindex] = newparams["phi"]
                                        self.content[i][tiltindex] = newparams["theta"]
                                        self.content[i][psiindex] = newparams["psi"]

        def get_micrograph_names(self):
                spl = self.split_by_micrograph_name()
                index = self.find_parameter_index("_rlnMicrographName")
                return [i[0][index] for i in spl]
        
        def get_image_names(self):
                index = self.find_parameter_index("_rlnImageName")
                return [cont[index] for cont in self.content]
                        
        def get_parameter_names(self):
                params = []
                for i in self.header:                
                        if i.startswith("_rln"): params.append(i)
                return params

        def angles_to_prior(self):
                rotindex = self.find_parameter_index("_rlnAngleRot")
                tiltindex = self.find_parameter_index("_rlnAngleTilt")
                psiindex = self.find_parameter_index("_rlnAnglePsi")
                try:
                        rotpriorindex = self.find_parameter_index("_rlnAngleRotPrior")
                except ValueError:
                        self.header.append("_rlnAngleRotPrior")
                        for i in range(len(self.content)):
                                self.content[i].append(0)
                        rotpriorindex = len(self.content[i])-1
                try:
                        tiltpriorindex = self.find_parameter_index("_rlnAngleTiltPrior")
                except ValueError:
                        self.header.append("_rlnAngleTiltPrior")
                        for i in range(len(self.content)):
                                self.content[i].append(0)
                        tiltpriorindex = len(self.content[i])-1
                try:
                        psipriorindex = self.find_parameter_index("_rlnAnglePsiPrior")
                except ValueError:
                        self.header.append("_rlnAnglePsiPrior")
                        for i in range(len(self.content)):
                                self.content[i].append(0)
                        psipriorindex = len(self.content[i])-1
                for i in range(len(self.content)):
                        self.content[i][rotpriorindex] = self.content[i][rotindex]
                        self.content[i][tiltpriorindex] = self.content[i][tiltindex]
                        self.content[i][psipriorindex] = self.content[i][psiindex]

        def remove_this_image(self, image_name):
                index = self.find_parameter_index("_rlnImageName")
                content = [self.content[i] for i in range(len(self.content)) if self.content[i][index][6:] != image_name]
                self.content = content

        def remove_outside_image(self):
                xindex = self.find_parameter_index("_rlnCoordinateX")
                yindex = self.find_parameter_index("_rlnCoordinateY")
                shxindex = self.find_parameter_index("_rlnOriginX")
                shyindex = self.find_parameter_index("_rlnOriginY")
                micindex = self.find_parameter_index("_rlnMicrographName")
                contcopy = copy.deepcopy(self.content)
                tmp = 0
                h, w = EMNumPy.em2numpy(EMData(self.content[0][micindex])).shape
                for i in range(len(contcopy)):
                        centx = contcopy[i][xindex]-contcopy[i][shxindex]
                        centy = contcopy[i][yindex]-contcopy[i][shyindex]
                        if centx < 0 or centx > w or centy < 0 or centy > h:
                                del self.content[i-tmp]
                                tmp += 1
                return True

        def count_particles(self):
                return len(self.content)
        
        def extract_only_these_parameters(self, parameters):
                indexes = [self.find_parameter_index(parname) for parname in parameters]
                newcontent = []
                temp = []
                for i in range(self.count_particles()):
                        for index in indexes:
                                temp.append(self.content[i][index])
                        newcontent.append(temp)
                        temp = []
                return newcontent

        def add_parameter(self, parameter_name):
                for i in range(len(self.content)):
                        self.content[i].append(0)
                last = 0
                for i in range(len(self.header)):
                        if self.header[i].startswith("_rln"): last = i
                self.header = self.header[:last+1] + [parameter_name] + self.header[last+1:]

        def write(self, filename):
                content = "\n".join(["\t".join([str(i) for i in cont]) for cont in self.content])
                open(filename, "w").write("\n".join(self.header)+"\n\n"+content)

def create_header(parameters):
        return ["data_","loop_",""]+parameters

def write_star(filename, header, content):
        content = "\n".join(["\t".join([str(i) for i in cont]) for cont in content])
        open(filename, "w").write("\n".join(header)+"\n\n"+content)

def join_star(star1, star2):
        if set(star1.get_parameter_names()) == set(star2.get_parameter_names()):
                joined = copy.deepcopy(star1)
                joined.content += star2.content
                return joined
        else:
                parnames = list(set(star1.get_parameter_names()) & set(star2.get_parameter_names()))
                star1_extracted_content = star1.extract_only_these_parameters(parnames)
                star2_extracted_content = star2.extract_only_these_parameters(parnames)
                header = create_header(parnames)
                content = star1_extracted_content + star2_extracted_content
                joined = copy.deepcopy(star1)
                joined.header = header
                joined.content = content
                return joined
