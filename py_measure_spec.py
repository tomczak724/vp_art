
import mypy
import glob as gl
import pylab as pl



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




atoms = [ 'Ar','Cd','He','Hg','Na','Ne' ]
#atoms = [ 'Hg' ]

for a in atoms:

    dat = pl.loadtxt('spectrum_'+a+'.dat')

    x1 = pl.find( dat[:,1] > pl.median(dat[:,1])+5*mypy.nmad(dat[:,1]) )
    l1 = dat[:,0][ pl.find( dat[:,1] > pl.median(dat[:,1])+50*mypy.nmad(dat[:,1]) ) ]
    f1 = dat[:,1][ pl.find( dat[:,1] > pl.median(dat[:,1])+50*mypy.nmad(dat[:,1]) ) ]

    xp,fp = sep_peaks( x1, dat[:,1] )
    lp = [ dat[:,0][ij] for ij in xp ]

    peaks_l = []
    peaks_f = []
    for i in range(len(xp)):
        line_l = lp[i]
        line_f = fp[i]

        if len(line_l)<5: continue

        fit_l, fit_f = mypy.gaussfit( line_l, line_f )

        peaks_l.append( fit_l[ fit_f.tolist().index(max(fit_f)) ] )
        peaks_f.append( max(fit_f)-pl.median(dat[:,1]) )

        outer = open('peaks_'+a+'.dat','w')
        for i in range(len(peaks_l)):
            if peaks_f[i]/max(peaks_f) > 0.02:
                outer.write( str(peaks_l[i])[:7]+'  '+str(peaks_f[i]/max(peaks_f))[:4]+'\n' )
        outer.close()

    fig = pl.figure( figsize=(24,6) )
    pl.plot( dat[:,0], dat[:,1], color='k' )
    pl.title( 'Spectrum of '+a )
    pl.xlabel( 'wavelength ($\AA$)' )
    pl.ylabel( 'counts' )

    if a=='Ne':
        peaks_l.append(6382.92)
        peaks_f.append(17000.0)
    for i in range(len(peaks_l)):
        l = peaks_l[i]
        if peaks_f[i]/max(peaks_f) > 0.02:
            pass
            #pl.axvline( l, color='k', ls=':' )

    d = max(dat[:,1])-pl.median(dat[:,1])
    pl.axis( [ 4500, 7000, max(dat[:,1])-d*1.05, max(dat[:,1]) + d*1.0/5 ] )
    mypy.tickthick()
    #pl.ylim( 100, pl.axis()[-1] )
    pl.savefig('spectrum_'+a+'.png')
    pl.clf()
