###  Once "vp_art_setup.py" has been run this script
###  should be run in the "working_dir"
###
###  NOTE TO PROGRAMMER: Always return to redux/
###                      after running each script.

import os
import time
import glob as gl
import numpy as np
import pylab as pl
import pyfits as pf
import shutil as shu
import decimal as deci
import subprocess as sp


###  Removing compiled scripts, in case there were changes made
pycs = gl.glob('scripts/*pyc')
for pyc in pycs: os.remove(pyc)

###  Adding vp_art scripts to PYTHONPATH
spcall = sp.call(['export PYTHONPATH="'+os.getenv('PYTHONPATH')+':'+os.getcwd()+'/scripts/"'],shell=True)

import vp_art_im_stack as im_stack
import vp_art_trace_fibers as trc
import vp_art_lambda_solve as lam
import vp_art_record_spectrum as spec
import vp_art_skysub as skysub
import vp_art_cosmicrays as cosmicrays

months = [ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' ]
t = time.gmtime()
h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
if len(h)==1: h='0'+h
if len(m)==1: m='0'+m
if len(s)==1: s='0'+s
print '\n\tBEGINNING VP_ART REDUCTION   -   '+months[t.tm_mon-1]+' '+str(t.tm_mday)+' '+str(t.tm_year)+'   GMT ('+h+':'+m+':'+s+')\n\n'



######################################################
###  Reading in vp_art.param... from the working_dir
######################################################
import vp_art_read_params as read_params
params = read_params.get_params()
master_targets = []
master_targets_txt = open('TARGETS.list')
for line in master_targets_txt.readlines():
    l = line.split()
    if l[0][0]!='#': master_targets.append(l[0])
master_targets_txt.close()



sky_txt = open('sky_fibers.dat','r')
sky_fibers = []
for line in sky_txt:
    l = line.split()
    if l[0][0]!='#':
        name = l[0]
        skyfibs = np.array( [ int(i) for i in l[1].split(',') ] )
        sky_fibers.append( [ name, skyfibs ] )


#########################
###  Reading LINE_LIST
#########################
try: linesdat = np.loadtxt(params.LINE_LIST)
except:
    raise IOError("can't find "+params.LINE_LIST)
    quit()
linesdat.sort(axis=0)



os.chdir('redux/')
nights = gl.glob('*')
for night in nights:


    ######################################################
    ###  switching into night's directory
    ######################################################
    os.chdir(night)
    night_targets = np.loadtxt('targets.txt',dtype=str)
    if np.size(night_targets)==1: night_targets = np.array([night_targets])
    redux_targets = []
    for nf in night_targets:
        if nf=='targets': continue
        for mf in master_targets:
            if nf==mf: redux_targets.append(nf)


    ######################################################
    ### bias combine: combining
    ######################################################
    if params.BIAS_COMBINE=='y':
        os.chdir('calib/bias/')
        bias_ims = gl.glob('*fits')
        bias_stk = im_stack.bias_stack(bias_ims)
        bias_exists = len(gl.glob('../../'+bias_stk))==1
        if bias_exists:
            print '\n\t   PREVIOUS BIAS-STACK ALREADY EXISTS... CLOBBERING\n'
            os.remove('../../'+bias_stk)
        shu.move(bias_stk,'../../')
        bias_exists = True
        os.chdir('../../')
    else:
        t = time.gmtime()
        h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
        if len(h)==1: h='0'+h
        if len(m)==1: m='0'+m
        if len(s)==1: s='0'+s
        print 'GMT ('+h+':'+m+':'+s+')  ---  SKIPPING BIAS COMBINE\n\n\n'
        bias_stk = im_stack.bias_stk
        bias_exists = len(gl.glob(bias_stk))==1


    ######################################################
    ### comp combine: combining
    ######################################################
    if params.COMP_COMBINE=='y':
        os.chdir('calib/comp/')
        comp_ims = gl.glob('*fits')
        comp_stk = im_stack.comp_stack(comp_ims,params.COMP_HEADER,params.COMP_NAME)
        comp_exists = len(gl.glob('../../'+comp_stk))==1
        if comp_exists:
            print '\n\t   PREVIOUS COMP-STACK ALREADY EXISTS... CLOBBERING\n'
            os.remove('../../'+comp_stk)
        comp_exists = True
        shu.move(comp_stk,'../../')
        os.chdir('../../')
    else:
        t = time.gmtime()
        h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
        if len(h)==1: h='0'+h
        if len(m)==1: m='0'+m
        if len(s)==1: s='0'+s
        print 'GMT ('+h+':'+m+':'+s+')  ---  SKIPPING COMP COMBINE\n\n\n'
        comp_stk = im_stack.comp_stk
        comp_exists = len(gl.glob(comp_stk))==1


    ######################################################
    ### flat combine: combining
    ######################################################
    if params.FLAT_COMBINE=='y':
        os.chdir('calib/flat/')
        flat_ims = gl.glob('*fits')
        try: flat_stk = im_stack.flat_stack(flat_ims)
        except: pass
        flat_exists = len(gl.glob('../../'+flat_stk))==1
        if flat_exists:
            print '\n\t   PREVIOUS FLAT-STACK ALREADY EXISTS... CLOBBERING\n'
            os.remove('../../'+flat_stk)
        flat_exists = True
        shu.move(flat_stk,'../../')
        os.chdir('../../')
    else:
        t = time.gmtime()
        h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
        if len(h)==1: h='0'+h
        if len(m)==1: m='0'+m
        if len(s)==1: s='0'+s
        print 'GMT ('+h+':'+m+':'+s+')  ---  SKIPPING FLAT COMBINE\n\n\n'
        flat_stk = im_stack.flat_stk
        flat_exists = len(gl.glob(flat_stk))==1


    ######################################################
    ### performing bias subtraction
    ######################################################
    if params.BIAS_SUBTRACT=='y' and not bias_exists:
        raise IOError('Master bias stack not found')

    if params.BIAS_SUBTRACT=='y':

        t = time.gmtime()
        h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
        if len(h)==1: h='0'+h
        if len(m)==1: m='0'+m
        if len(s)==1: s='0'+s
        print 'GMT ('+h+':'+m+':'+s+')  ---  PERFORMING BIAS SUBTRACTION\n'

        biasdat = pf.getdata(bias_stk)
        av = str(round(np.average(biasdat),2))
        st = str(round(np.std(biasdat),2))
        if len( av[av.index('.'):] )==2: av+='0'
        if len( st[st.index('.'):] )==2: st+='0'
        while len(av)<len(st): av = ' '+av
        while len(st)<len(av): st = ' '+st
        print '    BIAS STATS:'
        print '      mean    = '+ av
        print '      stddev  = '+ st +'\n'

        for target in redux_targets:
            os.chdir(target)
            ims = gl.glob('vp_art*fits')
            imstxt = open('images.list','w')
            for im in ims:
                imstxt.write(im[:im.index('.')]+'\n')
                fits = pf.open(im)
                biassub = fits[0].data - biasdat
                pf.writeto( im[:-5]+'_b.fits', biassub, fits[0].header, clobber=True )
                print '    writing to... ', target+'/'+im[:-5]+'_b.fits'
            imstxt.close()
            print ''
            os.chdir('../')

        t = time.gmtime()
        h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
        if len(h)==1: h='0'+h
        if len(m)==1: m='0'+m
        if len(s)==1: s='0'+s
        print '\nGMT ('+h+':'+m+':'+s+')  ---  COMPLETED BIAS SUBTRACTION\n\n\n'


    ######################################################
    ### tracing fibers from flat stack
    ######################################################
    fiber_folder = 'fiber_trace'
    if params.FIBER_TRACE=='y':
        fibers = trc.trace_fibers(flat_stk,params)

        trace_files = gl.glob(fiber_folder+'/fibertrace_id*')
        if len(trace_files)==0: fibers = trc.calc_trace_align(flat_stk,bias_stk,fibers,fiber_folder,params)
        else: fibers = trc.read_fibers_from_txt(trace_files,fiber_folder+'/fibertrace_ys.txt')
        print ''

        print '\twriting to... ', trc.do_trace_align(comp_stk, fibers, params) , '\n'
        for target in redux_targets:
            os.chdir(target)
            ims = gl.glob('*_b.fits')
            for im in ims: print '\twriting to... ', target+'/'+trc.do_trace_align(im,fibers,params)
            print ''
            os.chdir('../')
    else:
        trace_files = gl.glob(fiber_folder+'/fibertrace_id*')
        fibers = trc.read_fibers_from_txt(trace_files,fiber_folder+'/fibertrace_ys.txt')


    ######################################################
    ### calculating wavelength solution
    ######################################################
    if params.LAMBDA_SOLVE=='y':

        lambdafiles = gl.glob('lambda_solve/lambda_soln_*')

        if len(lambdafiles)==0:
            print '\n\tCALCULATING WAVELENGTH SOLUTIONS\n'
            lam.find_lines(comp_stk[:-5]+'_t.fits',fibers,linesdat,params)
        else:
            print '\n\tREADING WAVELENGTH SOLUTIONS\n'
            fibers = lam.read_lamsoln_from_txt('lambda_solve/solutions/lambda_soln_id',fibers)

        print '\n\tEVALUATING WAVELENGTH SLOUTIONS'
        print '\n\t   Calculating scatter in pixel coordinates of peaks' \
              '\n\t   found in the spectrum to inferred pixel locations' \
              '\n\t   from the polynomial fit.'
        print '\n\tFIBER   AVE       STD'
        for fib in fibers:
            xaxis = np.array( [ xy[0]*1.0 for xy in fib.xy ] )
            line_fits = np.array(fib.lines_lam_x)
            fib.lamaxis = lam.calc_wavlength_soln(comp_stk[:-5]+'_t.fits', \
                                                  fib, line_fits, xaxis, \
                                                  params.POLYFIT_ORDER)

    else:
        print '\n\tREADING WAVELENGTH SOLUTIONS\n'
        print '\n\tEVALUATING WAVELENGTH SLOUTIONS'
        print '\n\t   Calculating scatter in pixel coordinates of peaks' \
              '\n\t   found in the spectrum to inferred pixel locations' \
              '\n\t   from the polynomial fit.'
        print '\n\tFIBER   AVE       STD'
        fibers = lam.read_lamsoln_from_txt('lambda_solve/solutions/lambda_soln_id',fibers)
        for fib in fibers:
            xaxis = np.array( [ xy[0]*1.0 for xy in fib.xy ] )
            line_fits = np.array(fib.lines_lam_x)
            fib.lamaxis = lam.calc_wavlength_soln(comp_stk[:-5]+'_t.fits', \
                                                  fib, line_fits, xaxis, \
                                                  params.POLYFIT_ORDER)



    ######################################################
    ### cleaning out cosmic rays
    ######################################################
    if params.COSMIC_RAYS=='y':

        print '\n\tGENERATING COSMIC RAY MASKS'

        for target in redux_targets:
            os.chdir(target)
            data_ims = gl.glob('*_bt.fits')
            for im in data_ims: cosmicrays.clean_cosmicrays(im, params)
            os.chdir('../')


    ######################################################
    ### writing spectra to data files
    ######################################################
    if params.WRITE_SPECTRA=='y':

        t = time.gmtime()
        h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
        if len(h)==1: h='0'+h
        if len(m)==1: m='0'+m
        if len(s)==1: s='0'+s
        print 'GMT ('+h+':'+m+':'+s+')  ---  WRITING 1D SPECTRA\n\n\n'

        for target in redux_targets:
            os.chdir(target)
            try: os.mkdir('spectra/')
            except: pass


#  TRY READING IN COSMIC RAY CLEANED IMAGES, ELSE TAKE BIAS+TRACED
            data_ims = gl.glob('*_btc.fits')
            if len(data_ims)==0: data_ims = gl.glob('*_bt.fits')

            for data_im in data_ims:
                specfiles = spec.write_spectra(data_im, '../trace_pass2.fits', fibers)
                for s in specfiles: shu.move( s, 'spectra/' )
                print '\twrote spectra for', data_im, 'to', target +'/spectra/'
            os.chdir('../')


    ######################################################
    ### performing sky-subtraction
    ######################################################
    if params.SKY_SUBTRACT=='y':
        print '\n\tBEGINNING SKY SUBTRACTION\n'
        for target in redux_targets:
            os.chdir(target)
            ims = np.loadtxt('images.list',dtype=str)
            if np.size(ims)==1: ims = np.array([ims])

            os.chdir('spectra')
            skysub.pauper_skysub(ims)
            os.chdir('../../')


    '''
        print '\n\tBEGINNING SKY SUBTRACTION\n'
        for target in redux_targets:
            for skydat in sky_fibers:
                name = skydat[0]
                fibs = skydat[1]
                if name==target:
                    os.chdir(target)
                    ims = gl.glob('vp_*_bt.fits')
                    os.chdir('spectra/')
                    for im in ims:
                        skyfibs = [ 'spec_'+im[:-5]+'_fib'+spec.filenum(n,3)+'.dat' for n in fibs ]
                        median_skydat = skysub.sky_combine(skyfibs)
                        allspecs = gl.glob('spec_'+im[:-5]+'_fib*.dat')
                        for spe in allspecs:
                            imskysub = skysub.im_skysub(spe,median_skydat)
                            np.savetxt(spe.replace('bt','bts'),imskysub)
                        print '\twriting to... ', target+'/spec_'+im[:-5]+'s_fib*.dat'
                    os.chdir('../../')
    '''


                    

    ######################################################
    ### moving out of this night's directory
    ######################################################
    print ''
    os.chdir('../')


