import glob as gl
import numpy as np



###  This script reads in the given list of sky-fiber
###  spectra and median-combines them.
def sky_combine(sky_list):

    dat = np.loadtxt(sky_list[0])
    lam_axis = dat[:,0]
    flux_axes = [ dat[:,1] ]

    for skyfib in sky_list[1:]:

        dat = np.loadtxt(skyfib)
        flux = np.interp( lam_axis, dat[:,0], dat[:,1] )
        flux_axes.append( flux )

    flux_axis = np.median( flux_axes, axis=0 )
    return np.array( [ [lam_axis[i],flux_axis[i]] for i in range(len(lam_axis)) ] )



###  This script takes in the data spectrum and the
###  median sky spectrum from sky_combine() and
###  returns the sky-subtracted spectrum. The data
###  wavelength-axis is interpolated to the sky spectrum.
def im_skysub(spec,skyspec):

    dat = np.loadtxt(spec)
    rest_lam = skyspec[:,0]

    rest_flux = np.interp( rest_lam, dat[:,0], dat[:,1] )

    skysub_flux = rest_flux - skyspec[:,1]

    return np.array( [ [rest_lam[i],skysub_flux[i]] for i in range(len(rest_lam)) ] )
    


def pauper_skysub(ims):

    for im in ims:
        im_specs = gl.glob('spec_'+im+'*dat')
        im_spec_dats = np.array( [ np.loadtxt(s) for s in im_specs ] )

        min_lam, max_lam = max(im_spec_dats[:,0][:,0]), min(im_spec_dats[:,-1][:,0])
        lam_ax = np.linspace( min_lam, max_lam, len(im_spec_dats[0]) )
        im_spec_interps = np.array( [ np.interp( lam_ax, s[:,0], s[:,1] ) for s in im_spec_dats ] )

        median_stk = np.median( im_spec_interps, axis=0 )
        for i in range(len(im_specs)):
            outer = open('skysub_'+im_specs[i],'w')
            for x,y in zip(lam_ax,im_spec_interps[i]-median_stk):
                outer.write( str(x)+' '+str(y)+'\n' )
            outer.close()

