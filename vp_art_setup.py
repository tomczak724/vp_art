#  This script sorts the raw data into the
#  various directories for the reductions.

import os
import glob as gl
import numpy as np
import pyfits as pf
import shutil as shu
import subprocess as sp
from astLib import astCoords



###  Reading TARGETS.list
try: targets = np.loadtxt('TARGETS.list',dtype=str)
except:
    raise IOError("can't find TARGETS.list")
    quit()


###  These lines makes the code obsolete after 1 Jan 2100
class rawdir:
    def __init__(self,name):
        self.name = name
        os.chdir(name)
        self.images = gl.glob('*fits')
        self.bias = []
        self.comp = []
        self.flat = []
        self.obj = []
        for im in self.images:
            imtype = pf.getval(im,'IMAGETYP')
            if imtype=='comp': self.comp.append(im)
            if imtype=='zero': self.bias.append(im)
            if imtype=='flat': self.flat.append(im)
            if imtype=='object': self.obj.append(im)
        os.chdir('../')
            
os.chdir('raw/')
raw_dirs = [ rawdir(r) for r in gl.glob('20*') ]
os.chdir('../')


###  Making reduction directories
try: os.mkdir('redux/')
except: pass
os.chdir('redux/')

try:
    for raw_dir in raw_dirs: os.mkdir(raw_dir.name)
except: pass

###  THIS IS WHERE OBJECTS ARE FOUND FROM TARGETS.list
###  ... THEY MUST BE WITHIN 0.04deg (2.4') OF THE TARGETS.list VALUE
try:
    for raw_dir in raw_dirs:
        targetstxt = open('targets.txt','w')
        for target in targets:
            made_dir = False
            target_ra = astCoords.hms2decimal(target[1],':')
            target_dec = astCoords.dms2decimal(target[2],':')
            for im in raw_dir.obj:
                obj_ra = astCoords.hms2decimal(pf.getval('../raw/'+raw_dir.name+'/'+im,'RA'),':')
                obj_dec = astCoords.dms2decimal(pf.getval('../raw/'+raw_dir.name+'/'+im,'DEC'),':')
                if astCoords.calcAngSepDeg(target_ra,target_dec,obj_ra,obj_dec)<0.04:
                    if not made_dir:
                        os.mkdir(raw_dir.name+'/'+target[0])
                        made_dir = True
                        targetstxt.write(target[0]+'\n')
                    os.chdir(raw_dir.name+'/'+target[0])
                    ln_cmd = 'ln -s ../../../raw/'+raw_dir.name+'/'+im+' ./'
                    sp.call( [ln_cmd], shell=True )
                    os.chdir('../../')
        targetstxt.close()
        shu.move('targets.txt',raw_dir.name)
except: pass


try:
    for raw_dir in raw_dirs:
        os.mkdir(raw_dir.name+'/calib/')
        os.mkdir(raw_dir.name+'/fiber_trace/')
        os.mkdir(raw_dir.name+'/lambda_solve/')
        os.mkdir(raw_dir.name+'/lambda_solve/solutions/')
        os.mkdir(raw_dir.name+'/lambda_solve/comp_spectra/')
        os.mkdir(raw_dir.name+'/lambda_solve/figures/')
except: pass
try:
    for raw_dir in raw_dirs:
        os.mkdir(raw_dir.name+'/calib/bias/')
        os.chdir(raw_dir.name+'/calib/bias/')
        for im in raw_dir.bias:
            ln_cmd = 'ln -s ../../../../raw/'+raw_dir.name+'/'+im+' ./'
            sp.call( [ln_cmd], shell=True )
        os.chdir('../../../')
except: pass
try:
    for raw_dir in raw_dirs:
        os.mkdir(raw_dir.name+'/calib/comp/')
        os.chdir(raw_dir.name+'/calib/comp/')
        for im in raw_dir.comp:
            ln_cmd = 'ln -s ../../../../raw/'+raw_dir.name+'/'+im+' ./'
            sp.call( [ln_cmd], shell=True )
        os.chdir('../../../')
except: pass
try:
    for raw_dir in raw_dirs:
        os.mkdir(raw_dir.name+'/calib/flat/')
        os.chdir(raw_dir.name+'/calib/flat/')
        for im in raw_dir.flat:
            ln_cmd = 'ln -s ../../../../raw/'+raw_dir.name+'/'+im+' ./'
            sp.call( [ln_cmd], shell=True )
        os.chdir('../../../')
except: pass

