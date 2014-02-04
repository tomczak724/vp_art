
import mypy
import glob as gl
import pylab as pl
import pyfits as pf


'''
lam = pl.loadtxt('correct_x_lam.txt')[:,1]

fig = pl.figure()
c = 0
for a in atoms:
    flux = pf.getdata( gl.glob('*'+a+'*pass1.fits')[0] )[100]
    dat = [ [ lam[i], flux[i] ] for i in range(len(lam)) ]
    pl.savetxt( 'spectrum_'+a+'.txt', dat )

    c+=1
    pl.subplot(3,2,c)
    pl.plot( lam, flux )
    pl.title( a )
'''



#  Returns the x value at the peak of a gaussian fit to input x,y data
def find_xpeak( xdat, ydat ):
    xfit, yfit = mypy.gaussfit( xdat, ydat )
    return xfit[ yfit.tolist().index(max(yfit)) ]

#  Takes a single linear array of the peaks and returns an
#  array of arrays of the peaks.
#    EXAMPLE:  [1,2,3,7,8,9] --> [ [1,2,3], [7,8,9] ]
#  If dx between adjacent values is >5, starts a new peak
def sep_peaks( inx, iny ):
    outer_x = []
    outer_y = []
    tmp_x = [ inx[0] ]
    for i in inx[1:]:
        if i-tmp_x[-1]<5:
            tmp_x.append(i)
        else:
            outer_x.append(tmp_x)
            outer_y.append( iny[pl.array(tmp_x)] )
            tmp_x = [i]
    outer_x.append(tmp_x)
    outer_y.append( iny[pl.array(tmp_x)] )
    return outer_x, outer_y

#  Returns the wavelength axis
necd_lines = pl.loadtxt( 'NeCd.dat' )
def get_wvlnth_soln( xax, xpeaks ):
    f = pl.polyfit( xpeaks, necd_lines, 4 )
    def dldx(xn): return 4*f[0]*xn**3 + 3*f[1]*xn**2 + 2*f[2]*xn + f[3]
    lam = [ f[-1] ]
    for a in range(len(xax)-1): lam.append( lam[-1] + dldx(xax[a])*(xax[a+1]-xax[a]) )
    return lam



atoms = [ 'Ar','Cd','He','Hg','Na','Ne' ]

c=0
fig = pl.figure()
wvlnth_soln_dat = pf.getdata( 'NeCd.fits' )[3:-3]
for a in atoms:
    c += 1
    pl.subplot(3,2,c)

    restlam = pl.arange(4555,6820)
    restflux = pl.arange(4555,6820)*0.0
    imdat = pf.getdata( gl.glob('*'+a+'*pass1.fits')[0] )[3:-3]
    outerim = []

    for i in range(len(imdat)):

        wvlnth_flux = wvlnth_soln_dat[i]
        lines_x = pl.find( wvlnth_flux > pl.median(wvlnth_flux)+mypy.nmad(wvlnth_flux)*15 )
        lines_x, lines_f = sep_peaks( lines_x, wvlnth_flux )
        xpeaks = []
        for j in range(len(lines_x)): xpeaks.append( find_xpeak( lines_x[j], lines_f[j] ) )

        lam = get_wvlnth_soln( pl.arange(len(wvlnth_flux))+1, xpeaks )
        fiber_flux = imdat[i]

        interpflux = pl.interp( restlam, lam, fiber_flux )
        restflux += interpflux / len(imdat)
        outerim.append(interpflux)

    pf.writeto( 'arclamp_'+a+'_medi_spec_pass2.fits', pl.array(outerim) )

    pl.savetxt( 'spectrum_'+a+'.dat', [ [restlam[i],restflux[i]] for i in range(len(restlam)) ] )

    pl.plot( restlam, restflux )
    pl.title( a )
    print 'DONE WITH', a
