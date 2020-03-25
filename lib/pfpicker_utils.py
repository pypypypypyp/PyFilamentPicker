#!/usr/bin/python
#coding: utf-8

import numpy as np
from math import *
from scipy.fftpack import fft2
from scipy.misc import imsave
from scipy.interpolate import splrep, splprep, splev, interp1d
from scipy.ndimage import map_coordinates as interpolate2d
from scipy.ndimage import zoom, rotate, interpolation
from scipy.ndimage.filters import gaussian_laplace as log
from scipy.ndimage import gaussian_filter as gauss
from scipy.signal import *
from skimage.feature import hessian_matrix, hessian_matrix_eigvals
from skimage.filters import hessian, laplace
from skimage.restoration import denoise_bilateral
from progressbar import ProgressBar
from StarRW import *
from matplotlib import pyplot as plt
import gtk
import multiprocessing
import optparse
import shutil
import time

sys.path.append(os.getcwd())
import config

def loadGlobalVariables():
        global APIX, HEIGHT, ANGSTEP, DIST, SEGSTEP, ANGRANGE, STACK_BOXSIZE, FILAMENT_DIAMETER, INITIAL_ZOOM, INITIAL_THRES, SHRINK, LOCAL_FILAMENT_LENGTH, BIN, PERCENT_CONTRAST, MAX_ZOOM, HEIGHT_PIX, SPLINE_ACCURACY, SAVE
        reload(config)
        APIX = config.APIX
        HEIGHT = config.HEIGHT
        ANGSTEP = config.ANGSTEP
        DIST = config.DIST
        SEGSTEP = config.SEGSTEP
        ANGRANGE = config.ANGRANGE
        STACK_BOXSIZE = config.STACK_BOXSIZE
        FILAMENT_DIAMETER = config.FILAMENT_DIAMETER
        INITIAL_ZOOM = config.INITIAL_ZOOM
        INITIAL_THRES = config.INITIAL_THRES
        SHRINK = config.SHRINK
        LOCAL_FILAMENT_LENGTH = config.LOCAL_FILAMENT_LENGTH
        BIN = config.BIN
        PERCENT_CONTRAST = config.PERCENT_CONTRAST
        SAVE = config.SAVE
        MAX_ZOOM = 100/BIN
        HEIGHT_PIX = int(HEIGHT/APIX)
        if HEIGHT_PIX % 2 == 1: HEIGHT_PIX += 1 # HEIGHT_PIX should be an even number
        if LOCAL_FILAMENT_LENGTH % 2 == 1: LOCAL_FILAMENT_LENGTH += 1 # LOCAL_FILAMENT_LENGTH should be an even number
        SPLINE_ACCURACY = 2 # points along spline curve is plotted at the distance of 1/SPLINE_ACCURACY pixels

def toUINT8(array): # expand the signal to 0-255
        array = ((array - array.min()) / (array.ptp() / 255.)).astype(np.uint8)
        return array

def save(data, name): # to save images for debug
        EMNumPy.numpy2em(data.astype(np.float32)).write_image(name)

def circularMask(img, radius):
        size = img.shape[0]
        y, x = np.ogrid[-size/2:size/2, -size/2:size/2]
        mask = x**2 + y**2 > radius**2
        img[mask] = 0
        return img

def findShiftToCenter(data, limit=None, ref=None):
        data_pad = np.pad(data, (data.size/2, data.size/2), "constant", constant_values=(0, 0))
        if ref is None:
                cor = correlate(data, -data[::-1], mode="same")
        else:
                cor = correlate(data, ref, mode="same")
        if limit is not None and cor.size/2-limit >= 0:
                cor = cor[cor.size/2-limit:cor.size/2+limit]
        if ref is None:
                shift = (np.argmax(cor)-cor.size/2.)/2.
        else:
                shift = np.argmax(cor)-cor.size/2.
        return shift, np.max(cor)

def tubularMask(ref, od, id): # for white tube with black background
        x, y = np.meshgrid(np.arange(-ref.shape[0]/2, ref.shape[0]/2), np.arange(-ref.shape[1]/2, ref.shape[1]/2))
        mask_inner = np.abs(x) < id/APIX/BIN/2
        ref[mask_inner] = 0
        mask_outer = np.abs(x) > od/APIX/BIN/2
        ref[mask_outer] = 0
        return ref

def oneDimensionalReference(ref): # make the reference image to 1-dimensional image
        newref = np.zeros(ref.shape)
        refsum = np.sum(ref, axis=0)
        for i in range(refsum.size):
                newref[:, i] = refsum[i]
        return newref

def rfftshift(data): # np.fft.fftshift(data) for rfft
        h, w = data.shape
        ndata = np.zeros((h, w)).astype(data.dtype)
        if ndata.shape[0] % 2 == 0:
                ndata[:h/2, :] = data[h/2:, :]
                ndata[h/2:, :] = data[:h/2, :]
        else:
                ndata[:h/2, :] = data[h/2+1:, :]
                ndata[h/2:, :] = data[:h/2+1, :]
        return ndata

def lowpassFilter(img, res):
        WIDTH_FMASK_EDGE = 10
        rows, cols = img.shape
        size = 2 ** int(ceil(log(max(rows, cols), 2)))
        img = img - np.mean(img)
        nimg = np.zeros((size, size))
        nimg[:rows, :cols] = img
        foutput = np.fft.fft2(nimg)
        fshift = np.fft.fftshift(foutput)
        x, y = np.meshgrid(np.arange(-size/2, size/2), np.arange(-size/2, size/2))
        r = np.sqrt(x**2+y**2)
        radius = int(size / 2 * APIX * BIN / res)
        radius -= WIDTH_FMASK_EDGE / 2
        radius_p = radius + WIDTH_FMASK_EDGE
        mask_inner = r < radius
        mask_outer_cosedge = r >= radius
        mask_outer_zero = r > radius_p
        fmask = 0.5 + 0.5 * np.cos(np.pi*(r-radius)/WIDTH_FMASK_EDGE)
        fmask[mask_outer_zero] = 0
        fmask[mask_inner] = 1
        fshift *= fmask
        f_ishift = np.fft.ifftshift(fshift)
        img_back = np.fft.ifft2(f_ishift)
        img_back = img_back[:rows, :cols]
        return img_back

def spline(x, y, point, deg):
        tck, u = splprep([x, y], k=deg, s=0)
        u = np.linspace(0, 1, num=point, endpoint=True)
        spl = splev(u, tck)
        return spl[0], spl[1]

def straighten(positions, image):
        posx = positions[0]
        posy = positions[1]
        width = int(FILAMENT_DIAMETER/APIX/BIN*3)
        # calculate distance from the top
        diff = np.sqrt(np.diff(posx)**2+np.diff(posy)**2)
        diff = np.insert(diff, 0, 0)
        dist_from_top = np.array([np.sum(diff[:i+1]) for i in range(posx.size)])
        # calculate angles
        angles = [atan2(posy[i+1]-posy[i], posx[i+1]-posx[i])-np.pi/2 for i in range(posx.size-1)]
        angles.append(angles[-1])
        # create blank image with the same size as the straightened image
        straightened = np.zeros((int(dist_from_top[-1])+1, width))
        x, y = np.meshgrid(np.arange(-width/2, width/2), np.arange(0, 1))
        current_dist = 0
        for i in range(posx.size):
                if dist_from_top[i] >= current_dist:
                        xx = posx[i] + x*np.cos(angles[i]) - y*np.sin(angles[i])
                        yy = posy[i] + x*np.sin(angles[i]) - y*np.cos(angles[i])
                        straightened[current_dist, :] = interpolate2d(image, [yy, xx], order=0).reshape(x.shape)
                        current_dist += 1
        return straightened

def createSplineCurve(x, y, popend = False):
        dist = DIST / APIX / BIN
        spl = spline(x, y, (x.size-1)*dist*SPLINE_ACCURACY, 2)
        if popend:
                np.delete(spl[0], -1)
                np.delete(spl[1], -1)
        return spl[0], spl[1]

def createSplineCurveForAllLength(positions, seg=True):
        positions = np.array([np.array(i) for i in positions])
        x = positions[:, 0]
        y = positions[:, 1]
        if not seg:
                return createSplineCurve(x, y)
        if x.size < LOCAL_FILAMENT_LENGTH: return
        n = (x.size-1)/(LOCAL_FILAMENT_LENGTH-1)-1
        spline_locals = [createSplineCurve(x[i*(LOCAL_FILAMENT_LENGTH-1):(i+1)*(LOCAL_FILAMENT_LENGTH-1)+1], y[i*(LOCAL_FILAMENT_LENGTH-1):(i+1)*(LOCAL_FILAMENT_LENGTH-1)+1], popend = True) for i in range(n)]
        if n != 0:
                local_filament_length = spline_locals[0][0].size
        else:
                local_filament_length = 0
        spline_locals.append(createSplineCurve(x[n*(LOCAL_FILAMENT_LENGTH-1):], y[n*(LOCAL_FILAMENT_LENGTH-1):]))
        spline_all_x, spline_all_y = np.zeros((2, local_filament_length*n + spline_locals[-1][0].size))
        for i in range(n):
                spline_all_x[i*local_filament_length:(i+1)*local_filament_length] = spline_locals[i][0]
                spline_all_y[i*local_filament_length:(i+1)*local_filament_length] = spline_locals[i][1]
        spline_all_x[n*local_filament_length:] = spline_locals[n][0][:spline_all_x.size-n*local_filament_length]
        spline_all_y[n*local_filament_length:] = spline_locals[n][1][:spline_all_y.size-n*local_filament_length]
        return spline_all_x, spline_all_y


def createSplineCurveForAllLength2(positions): # at this point positions are gives as [(x0, y0), (x1, y1), ...]
        positions = np.array([np.array(i) for i in positions])
        x = positions[:, 0]
        y = positions[:, 1]
        # raise error if the filament is too short
        if x.size < LOCAL_FILAMENT_LENGTH: return
        # split the filament into overlapping segments and calculate spline curve for each segment
        n = x.size/(LOCAL_FILAMENT_LENGTH/2)-1
        spline_locals = [createSplineCurve(x[i*LOCAL_FILAMENT_LENGTH/2:(i+2)*LOCAL_FILAMENT_LENGTH/2], y[i*LOCAL_FILAMENT_LENGTH/2:(i+2)*LOCAL_FILAMENT_LENGTH/2]) for i in range(n)]
        local_filament_length = spline_locals[0][0].size
        edge_is_covered = x.size == (n+1)*LOCAL_FILAMENT_LENGTH/2
        if not edge_is_covered:
                spline_locals.append(createSplineCurve(x[n*LOCAL_FILAMENT_LENGTH/2:], y[n*LOCAL_FILAMENT_LENGTH/2:]))
        # combine
        spline_all_x, spline_all_y = np.zeros((2, spline_locals[0][0].size/2*n + (spline_locals[0][0].size/2 if edge_is_covered else spline_locals[-1][0].size)))
        numSeg = len(spline_locals)
        for i in range(numSeg-1):
                spline_all_x[i*local_filament_length/2:i*local_filament_length/2+spline_locals[i][0].size] = spline_locals[i][0]
                spline_all_y[i*local_filament_length/2:i*local_filament_length/2+spline_locals[i][1].size] = spline_locals[i][1]
        spline_all_x[(numSeg-1)*local_filament_length/2:] = spline_locals[numSeg-1][0][:spline_all_x.size-(numSeg-1)*local_filament_length/2]
        spline_all_y[(numSeg-1)*local_filament_length/2:] = spline_locals[numSeg-1][1][:spline_all_y.size-(numSeg-1)*local_filament_length/2]
        return spline_all_x, spline_all_y

def flatten(img):
        imgmean = img.mean()
        means = [img[i].mean() for i in range(img.shape[0])]
        img = [img[i]-(means[i]-imgmean) for i in range(img.shape[0])]
        return np.array(img)

def toText(strlist):
        s = "\n".join(["> "+i for i in strlist])
        return s

def prepareOutput(pos, displayables):
        dist_from_top = [0]+[sqrt((pos[i+1][0]-pos[i][0])**2+(pos[i+1][1]-pos[i][1])**2) for i in range(len(pos)-1)]
        dist_from_top = np.array([sum(dist_from_top[:i+1]) for i in range(len(pos))])
        n = 0
        indices = []
        precise_angles = [np.rad2deg(atan2(pos[i+1][1]-pos[i][1], pos[i+1][0]-pos[i][0])-np.pi/2) for i in range(len(pos)-1)]
        precise_angles.append(precise_angles[-1])
        for i in range(len(pos)):
                if dist_from_top[i] >= SEGSTEP/APIX*n:
                       indices.append(i)
                       n += 1
        pos = np.array(pos)
        precise_angles = np.array(precise_angles)
        positions = pos[indices]
        angles = precise_angles[indices]
        dist_from_top = dist_from_top[indices]
        displayables = np.array(displayables)
        displayable_start = np.where(displayables == True)[0][0] * DIST/APIX
        displayable_start = np.argmin(np.abs(dist_from_top - displayable_start))
        displayable_end = list(np.where(displayables == True)[0])[-1] * DIST/APIX + SEGSTEP/APIX
        displayable_end = np.argmin(np.abs(dist_from_top - displayable_end))
        displayables = [displayable_start <= i and i <= displayable_end for i in range(dist_from_top.size)]
        return positions, angles, displayables

def symmetrize(ref_str):
        ref_str_flipped = np.fliplr(ref_str)
        return (ref_str+ref_str_flipped)/2

def sortbyabs(array, axis=0, auxiliary=None):
        index = list(np.ix_(*[np.arange(i) for i in array.shape]))
        index[axis] = np.abs(array).argsort(axis)
        if auxiliary is None:
                return array[tuple(index)]
        else:
                return array[tuple(index)], auxiliary[tuple(index)]

def divide_nonzero(array1, array2, cval=1e-10):
        denom = np.copy(array2)
        denom[denom == 0] = cval
        return np.divide(array1, denom)

def FDAGK(sigma, theta, rho):
        height = int(4 * sigma * rho)
        width = int(4 * sigma)
        theta = np.deg2rad(theta)
        x, y = np.meshgrid(np.arange(-(width/2), width/2+1), np.arange(-(height/2), height/2+1))
        t1 = 1 * np.cos(theta) * (x*np.cos(theta) + y*np.sin(theta)) - rho**-2 * np.sin(theta) * (-x*np.sin(theta) + y*np.cos(theta))
        t2 = 1 * np.sin(theta) * (x*np.cos(theta) + y*np.sin(theta)) + rho**-2 * np.cos(theta) * (-x*np.sin(theta) + y*np.cos(theta))
        phi = x*t1 + y*t2
        AGK = 1./(2 * np.pi * rho * sigma**2) * np.exp(-phi/(2*sigma**2))
        FOAGK = -1 * ((x*np.cos(theta)+y*np.sin(theta)) / sigma**2) * AGK
        return FOAGK

def meijering(image, black_ridges=True):
        # sigma equals the radius of the ridge being emphasized
        sigma = (FILAMENT_DIAMETER - FILAMENT_INNER_DIAMETER) / APIX / BIN / 2
        image = image.astype(np.float64)
        if black_ridges: image = -image
        value = np.zeros(image.shape)
        Hxx, Hxy, Hyy = hessian_matrix(image, sigma=sigma, order="rc")
        b1 = Hxx + Hyy
        b2 = Hxx - Hyy
        d = np.sqrt(4*Hxy*Hxy + b2*b2)
        L1 = (b1 + 2*d) / 3.0
        L2 = (b1 - 2*d) / 3.0
        vect1x = b2 + d
        vect2x = b2 - d
        vecty = 2 * Hxy
        vectx = np.array([vect1x, vect2x])
        sortedL, sortedvectx = sortbyabs(np.array([L1, L2]), auxiliary=vectx)
        L = sortedL[1]
        vectx = sortedvectx[0]
        vectlen = np.sqrt(vectx**2 + vecty**2)
        vectx /= vectlen
        vecty /= vectlen
        valL = np.where(L > 0, 0, L)
        valL = divide_nonzero(valL, np.min(valL))
        vect = np.array([vectx, vecty])
        return valL, vect

def zoom_up(img):
        rows, cols = img.shape
        rows_bin, cols_bin = rows/BIN*BIN, cols/BIN*BIN
        img = img[:rows_bin, :cols_bin]
        return img.reshape(rows_bin/BIN, BIN, cols_bin/BIN, BIN).sum(3).sum(1)

def readMrcFile(files, i):
        e = EMData()
        e.read_image(files[i])
        img = EMNumPy.em2numpy(e)
        origimg = img.copy()
        #binimg = cutoffNoise(img, PERCENT_CONTRAST)
        #binimg = interpolation.zoom(img, 1./BIN, order=1)
        binimg = zoom_up(img)
        binimg = cutoffNoise(binimg, PERCENT_CONTRAST)
        binimg = lowpassFilter(binimg, 20.)
        binimg = np.real(binimg)
        binimg = toUINT8(binimg)
        directory = files[i][3:].replace(".", "-").replace("/", "_")
        emname = files[i]
        if not os.path.exists(directory): os.system("mkdir %s"%directory)
        return origimg, binimg, directory, emname

def cutoffNoise(array, percentile_contrast): #percentile-contrast
        imgtype = array.dtype
        high = np.percentile(array, 100-percentile_contrast)
        low = np.percentile(array, percentile_contrast)
        high_values = array > high
        low_values = array < low
        array[high_values] = high
        array[low_values] = low
        array = array.astype(imgtype)
        return array

loadGlobalVariables()
