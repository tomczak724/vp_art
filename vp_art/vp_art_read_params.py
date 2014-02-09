###  This script reads in the parameter-file "vp_art.param"
###  and prepares a python-ready structure array.

import pylab as pl

###  This is the structure definition
class params:
    def __init__(self, paramfile):
        for row in paramfile:

            if row[0]=='BIAS_COMBINE': self.BIAS_COMBINE = row[1]
            if row[0]=='COMP_COMBINE': self.COMP_COMBINE = row[1]
            if row[0]=='FLAT_COMBINE': self.FLAT_COMBINE = row[1]
            if row[0]=='FIBER_TRACE': self.FIBER_TRACE = row[1]
            if row[0]=='BIAS_SUBTRACT': self.BIAS_SUBTRACT = row[1]
            if row[0]=='LAMBDA_SOLVE': self.LAMBDA_SOLVE = row[1]
            if row[0]=='COSMIC_RAYS': self.COSMIC_RAYS = row[1]
            if row[0]=='WRITE_SPECTRA': self.WRITE_SPECTRA = row[1]
            
            if row[0]=='FIBERS_TOTAL': self.FIBERS_TOTAL = int(row[1])
            if row[0]=='FIBERS_EXCLUDE':
                fibs = pl.array( [ int(i) for i in row[1].split(',') ] )
                fibs.sort()
                fibs = fibs[ pl.find( 0<fibs ) ]
                fibs = fibs[ pl.find( fibs<self.FIBERS_TOTAL+1 ) ]
                self.FIBERS_EXCLUDE = fibs
            if row[0]=='FIBER_WIDTH': self.FIBER_WIDTH = int(row[1])
            if row[0]=='FIBER_SEP': self.FIBER_SEP = int(row[1])
            if row[0]=='LO_BUFFER': self.LO_BUFFER = int(row[1])
            if row[0]=='HI_BUFFER': self.HI_BUFFER = int(row[1])
            if row[0]=='COMP_HEADER': self.COMP_HEADER = row[1]
            if row[0]=='COMP_NAME': self.COMP_NAME = row[1]
            if row[0]=='POLYFIT_ORDER': self.POLYFIT_ORDER = int(row[1])
            if row[0]=='CR_SCAN_DX': self.CR_SCAN_DX = int(row[1])

###  This routine spits out the parameter array when you're ready
def get_params(): return params( pl.loadtxt('vp_art.param',dtype=str) )


