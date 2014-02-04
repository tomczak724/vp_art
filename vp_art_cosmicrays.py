#  This script identifies cosmic rays in individual images.
#  Images must be trace-aligned.
#

import cosmics
import pyfits as pf


def clean_cosmicrays( image, params ):

    imdat, imhead = cosmics.fromfits(image, verbose=False)

    # DEFAULT HEADER PARAMETERS FOR VIRUS-P AS OF 2012
    try: gain = imhead['GAIN1']
    except: gain = 2.0
    try: rdnoise = imhead['RDNOISE1']
    except: rdnoise = 2.5

    c = cosmics.cosmicsimage( imdat, gain=gain, readnoise=rdnoise, \
                              sigclip=5., sigfrac=0.3, objlim=5.0 )
    c.run( maxiter=4 )

    cosmics.tofits(image[:-5]+'c.fits', c.cleanarray, imhead)
    cosmics.tofits(image[:-5]+'_mask.fits', c.mask, imhead)

