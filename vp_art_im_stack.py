import time
import numpy as np
import pyfits as pf

def nmad( array ):
    """
  #  Returns the normalized median absolute
  #  deviation of a the given array
    """
    return 1.483*np.median( abs( np.array(array) - np.median(array) ) )



bias_stk = 'stk_bias.fits'
comp_stk = 'stk_comp.fits'
flat_stk = 'stk_flat.fits'

def bias_stack(images):

    t = time.gmtime()
    h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
    if len(h)==1: h='0'+h
    if len(m)==1: m='0'+m
    if len(s)==1: s='0'+s
    print 'GMT ('+h+':'+m+':'+s+')  ---  COMBINING BIAS FRAMES\n'

    bias_data_cube = np.array( [ pf.getdata(im) for im in images ] )
    bias_head = pf.getheader(images[0])

    namelen = max( [ len(im) for im in images ] ) +7
    if namelen-7<len(bias_stk): namelen = len(bias_stk) +7

    av = str(round(np.average(bias_data_cube[0]),2))
    st = str(round(np.std(bias_data_cube[0]),2))
    if len( av[av.index('.'):] )==2: av+='0'
    if len( st[st.index('.'):] )==2: st+='0'
    me = str(round(np.median(bias_data_cube[0]),2))
    nm = str(round(nmad(bias_data_cube[0]),2))
    if len( me[me.index('.'):] )==2: me+='0'
    if len( nm[nm.index('.'):] )==2: nm+='0'
    statlen = max( len(av), len(st), len(me), len(nm) ) +1
    if statlen<4: statlen=4

    title = 'FRAME'
    while len(title)<namelen: title+=' '
    title+='  AVE'
    for i in range(statlen-3+2): title+=' '
    title+='STD'
    for i in range(statlen-3+2): title+=' '
    title+='MEDI'
    for i in range(statlen-4+2): title+=' '
    print '    '+title+'NMAD'

    for i in range(len(images)):
        av = str(round(np.average(bias_data_cube[i]),2))
        st = str(round(np.std(bias_data_cube[i]),2))
        if len( av[av.index('.'):] )==2: av+='0'
        if len( st[st.index('.'):] )==2: st+='0'
        me = str(round(np.median(bias_data_cube[i]),2))
        nm = str(round(nmad(bias_data_cube[i]),2))
        if len( me[me.index('.'):] )==2: me+='0'
        if len( nm[nm.index('.'):] )==2: nm+='0'

        while len(av)<statlen: av+=' '
        while len(st)<statlen: st+=' '
        while len(me)<statlen: me+=' '
        while len(nm)<statlen: nm+=' '

        s0 = images[i]+'  '
        while len(s0)<namelen: s0+='-'
        s0+='  '
        s0 += av+'  '+st+'  '
        s0 += me+'  '+nm
        print '    '+s0


    stack = np.median(bias_data_cube,axis=0)

    av = str(round(np.average(stack),2))
    st = str(round(np.std(stack),2))
    if len( av[av.index('.'):] )==2: av+='0'
    if len( st[st.index('.'):] )==2: st+='0'
    me = str(round(np.median(stack),2))
    nm = str(round(nmad(stack),2))
    if len( me[me.index('.'):] )==2: me+='0'
    if len( nm[nm.index('.'):] )==2: nm+='0'

    while len(av)<statlen: av+=' '
    while len(st)<statlen: st+=' '
    while len(me)<statlen: me+=' '
    while len(nm)<statlen: nm+=' '

    s0 = bias_stk+'  '
    while len(s0)<namelen: s0+='-'
    s0+='  '
    s0 += av+'  '+st+'  '
    s0 += me+'  '+nm
    print '\n    '+s0

    pf.writeto(bias_stk, stack, header=bias_head)


    t = time.gmtime()
    h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
    if len(h)==1: h='0'+h
    if len(m)==1: m='0'+m
    if len(s)==1: s='0'+s
    print '\nGMT ('+h+':'+m+':'+s+')  ---  COMPLETED BIAS STACKING\n\n\n'
    return bias_stk






def comp_stack(images,comp_header,comp_name):

    t = time.gmtime()
    h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
    if len(h)==1: h='0'+h
    if len(m)==1: m='0'+m
    if len(s)==1: s='0'+s
    print 'GMT ('+h+':'+m+':'+s+')  ---  COMBINING COMP FRAMES\n'

    objects = [ pf.getheader(im)[comp_header] for im in images ]
    comp_images = np.array(images)[np.where(np.array(objects)==comp_name)[0]]
    title = 'FOUND '+str(len(comp_images))+' '+comp_name+' FRAMES'
    print '    '+title

    comp_data_cube = np.array( [ pf.getdata(im) for im in comp_images ] )
    comp_head = pf.getheader(comp_images[0])
    stack = np.median(comp_data_cube,axis=0)
    pf.writeto(comp_stk, stack, header=comp_head)

    t = time.gmtime()
    h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
    if len(h)==1: h='0'+h
    if len(m)==1: m='0'+m
    if len(s)==1: s='0'+s
    print '\nGMT ('+h+':'+m+':'+s+')  ---  COMPLETED COMP STACKING\n\n\n'
    return comp_stk





def flat_stack(images):

    t = time.gmtime()
    h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
    if len(h)==1: h='0'+h
    if len(m)==1: m='0'+m
    if len(s)==1: s='0'+s
    print 'GMT ('+h+':'+m+':'+s+')  ---  COMBINING FLAT FRAMES\n'

    flat_data_cube = np.array( [ pf.getdata(im) for im in images ] )
    flat_head = pf.getheader(images[0])
    stack = np.median(flat_data_cube,axis=0)
    pf.writeto(flat_stk, stack, header=flat_head)

    t = time.gmtime()
    h,m,s = str(t.tm_hour), str(t.tm_min), str(t.tm_sec)
    if len(h)==1: h='0'+h
    if len(m)==1: m='0'+m
    if len(s)==1: s='0'+s
    print '\nGMT ('+h+':'+m+':'+s+')  ---  COMPLETED FLAT STACKING\n\n\n'
    return flat_stk


