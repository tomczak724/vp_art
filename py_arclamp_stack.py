
#  DESCRIPTION:
#    I wrote this script to take an image from VIRUS-P and (crudely)
#    align the fibers along the wavelength axis. Then, once the spectra
#    are aligned, to stack them all together. This script was designed
#    to make stacks of arclamp spectra.
#
#  INPUTS:
#    (1) A high S/N flat image to trace out the fibers
#    (2) The data image (arclamp, as mentioned above)
#
#  OUTPUTS:
#    (1) An image of all the aligned 1d spectra
#    (2) A txt of the final stacked spectrum
#
#  CALLING SEQUENCE:
#    [shell] >>>  python py_arclamp_stack.py data_im trace_im
#        data_im  = name of the data image
#        trace_im = name of the flat image
#
#  NOTES:
#    (1) Due to the positioning of fibers of VIRUS
#        (when I wrote this), fibers 1,2,3 & 4 are
#        either partially or entirely off of the CCD.
#        Thus, I only consier the last 200 spectra.

import sys
import time as tm
import pylab as pl
import pyfits as pf

def reverse(arr): return arr[::-1]
def timestamp(): return str(tm.localtime().tm_hour)+':'+str(tm.localtime().tm_min)+':'+str(tm.localtime().tm_sec)

#
#dat = pf.getdata(sys.argv[1])
#trc = pf.getdata(sys.argv[2])
dat = pf.getdata('arclamp_Ne_medi.fits')
trc = pf.getdata('stack_twiflat.fits')

class fiber:
    def __init__(self, xaxis):
        self.xaxis = pl.array(xaxis)*0
        self.y0 = -1
        self.xpeaks = []
        self.flux = pl.arange(len(xaxis))*0
        self.s,self.m = 0,0
fibers = [ fiber(range(len(dat))) for i in range(200) ]


for x in range(2048):
    dat_col = dat[:,x]
    trc_col = trc[:,x]

# FIRST WE FIND THE WHICH Y-COORD A FIBER IS AT FOR GIVEN X-COORD
    trc_above_medi = pl.find( trc_col > pl.median(trc_col) )
    trc_above_medi_rev = reverse(trc_above_medi)

    fib_cnt = 0
    for i in range(len(trc_above_medi_rev)):

# HERE WE LOOK AT THE Y-COORDS AROUND A FIBER...
        ind = trc_above_medi_rev[i]
        slope0 = trc_col[ind] - trc_col[ind-1]
        slope1 = trc_col[ind+1] - trc_col[ind]
# ... AND FIND THE PEAK BY WHERE dFLUX/dY TURNS OVER
        if slope0>=0 and slope1<=0:
            fib = fibers[fib_cnt]
            flx = dat_col[ind]
            if fib.y0==-1: fib.y0 = ind
            else:
                dy = ind - fib.y0
                dr = pl.sqrt( 1**2 + dy**2 )
                fib.xaxis[x] = fib.xaxis[x-1] +1# + dr
                #                if fib_cnt==3: print x, dy
                fib.y0 = ind
            fib.flux[x] = flx
            fib_cnt += 1
            if fib_cnt==200: break

    if x % (len(dat[0])/10) == 0: print 'PROGRESS: ', timestamp(),'--', x,'/', len(dat[0])


outer_pass1 = []
outer_pass1.append( pl.arange(len(fib.flux))*0 )
outer_pass1.append( pl.arange(len(fib.flux))*0 )
outer_pass1.append( pl.arange(len(fib.flux))*0 )
for f in fibers: outer_pass1.append( f.flux )
outer_pass1.append( pl.arange(len(fib.flux))*0 )
outer_pass1.append( pl.arange(len(fib.flux))*0 )
outer_pass1.append( pl.arange(len(fib.flux))*0 )

try: pf.writeto(sys.argv[1][:-5]+'_spec_pass1.fits', pl.array(outer_pass1) )
except: pass







outer_pass2 = []
outer_pass2.append( pl.arange(len(fib.flux))*0 )
outer_pass2.append( pl.arange(len(fib.flux))*0 )
outer_pass2.append( pl.arange(len(fib.flux))*0 )

n = 100
tmp0 = fibers[0]
tmpn = fibers[n]
fibers[0] = tmpn
fibers[n] = tmp0

stack_x = fibers[0].xaxis
stack_f = fibers[0].flux

for f in fibers:

# FINDS THE EMISSION LINES, WHERE FLUX IS ABOVE THE AVERAGE
    flx_above_medi0 = pl.find( f.flux > pl.average(f.flux) )
    flx_above_medi1 = []

# FINDS THE EMISSION THAT HAVE MORE THAN 5 PIXELS ABOVE
# AVERAGE FLUX. THEN RECORDS THE AVERAGE X-COORD FOR THE LINE.
    tmp = [ flx_above_medi0[0] ]
    for i in range(1,len(flx_above_medi0)):
        if flx_above_medi0[i] - tmp[-1] < 2: tmp.append( flx_above_medi0[i] )
        else:
            if len(tmp)>=5: f.xpeaks.append( pl.average(tmp) )
            tmp = [ flx_above_medi0[i] ]

    if f==fibers[0]: continue

# NOW WE NEED TO FIND THE RESCALING FACTOR (s) AND OFFSET (m)
    f.s = pl.average( (pl.array(fibers[0].xpeaks)-fibers[0].xpeaks[0]+1) / (pl.array(f.xpeaks)-f.xpeaks[0]+1) )
    f.m = fibers[0].xpeaks[0] - f.xpeaks[0]*f.s
#   pl.fill_between( f.xaxis*f.s+f.m, f.flux, 0, color='k', alpha=0.01 )

    outer_pass2.append( pl.interp( f.xaxis , f.xaxis*f.s+f.m , f.flux ) )
    stack_f += pl.interp( f.xaxis , f.xaxis*f.s+f.m , f.flux )

stack_f /= len(fibers)
outer_txt = open('outer_stack.txt','w')
for i in range(len(stack_f)): outer_txt.write( str(stack_x[i]) +'  '+ str(stack_f[i]) +'\n' )
outer_txt.close()

outer_pass2.append( pl.arange(len(fib.flux))*0 )
outer_pass2.append( pl.arange(len(fib.flux))*0 )
outer_pass2.append( pl.arange(len(fib.flux))*0 )

try: pf.writeto(sys.argv[1][:-5]+'_spec_pass2.fits', pl.array(outer_pass2) )
except: pass

