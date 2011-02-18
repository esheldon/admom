from __future__ import print_function
import sys
from sys import stdout
import numpy

from . import _admomf
from . import util

flagmap={0:'ok',
         2**0: 'sum0',
         2**1: 'shift',
         2**2: 'sum0',
         2**3: 'mom0',
         2**4: 'detm',
         2**5: 'detn',
         2**6: 'w0',
         2**7: 'maxit',
         2**8: 'detw0',
         2**9: 'badnsub'}


def admom(image, row, col, 
          sky=0.0, sigsky=1.0, 
          guess=None,
          shiftmax=5.0,
          nsub=4):
    """
    Name:
        admom
    Purpose:
        Calculate adaptive moments on the input image.
    Calling Sequence:
        res = admom(image, row, col, sky=0.0, sigsky=None, guess=None,
                    shiftmax=5.0)

    Inputs:
        image: 
            A two-dimensional image.  It will be converted if it is not 4-byte
            float or fortran order.
        row, col: 
            The center or centers for which moments are to be measured.
    Keywords:
        sky: 
            The sky value in the image.  Can be a scalar or array of the same
            length as the row,col centers. Default is 0.0
        sigsky: 
            The sigma of the sky in the image.  Can be a scalar or array of the
            same length as the row,col centers.  This is only used for error
            calculations.  Default is 1.0

        guess: 
            A guess for the second order moments, Ixx or Iyy.  The same
            value is used as a guess for both.  This would be ~sigma**2

        shiftmax: 
            Maximum allowed shift of the centroid in pixels.
        nsub:
            If objects are not well resolved, sigma <= 2.3 pixels,
            pixelization corrections will be applied by studying
            nsubXnsub sub-pixel grids.  You can turn off sub-pixel
            corrections by setting nsub=1

    Outputs:
        A dictionary containing
            Irr,Irc,Icc,e1,e2,uncer,rho4,whyflag,wrow,wcol
        As well as most of the inputs.
    """

    if len(image.shape) != 2:
        raise ValueError("image must be 2-dimensional")

    is_scalar=numpy.isscalar(row)

    row = numpy.array(row, ndmin=1, copy=True, dtype='f4')
    col = numpy.array(col, ndmin=1, copy=True, dtype='f4')
    # add 1 for fortran index
    row += 1
    col += 1
    
    if col.size != row.size:
        raise ValueError("row and col must be same size")

    sky,sigsky = get_sky_sigsky(row.size, sky, sigsky)

    Irr = numpy.zeros(row.size, dtype='f4')
    Irc = numpy.zeros(row.size, dtype='f4')
    Icc = numpy.zeros(row.size, dtype='f4')
    rho4 = numpy.zeros(row.size, dtype='f4')
    uncer = numpy.zeros(row.size, dtype='f4')

    numiter = numpy.zeros(row.size, dtype='i4')

    wrow = numpy.zeros(row.size, dtype='f4') - 9999.
    wcol = numpy.zeros(row.size, dtype='f4') - 9999.

    interpolated=numpy.zeros(row.size, dtype='i2')

    whyflag = numpy.zeros(row.size, dtype='i4')
    whystr = numpy.zeros(row.size, dtype='S5')

    if guess is not None:
        Irr[:] = guess
        Icc[:] = guess
            
    _admomf.ad_mom_sub(image,sky,sigsky,row,col,shiftmax,nsub,
                       Irr,Irc,Icc,rho4,wrow,wcol,uncer,numiter,interpolated,whyflag)

                     
    row -= 1
    col -= 1
    w,=numpy.where((wrow >= 1) & (wcol >= 1))
    if w.size > 0:
        wrow[w] -= 1
        wcol[w] -= 1

    e1 = numpy.zeros(row.size, dtype='f4') + -9999.0
    e2 = numpy.zeros(row.size, dtype='f4') + -9999.0
    a4 = numpy.zeros(row.size, dtype='f4') + -9999.0
    s2 = numpy.zeros(row.size, dtype='f4') + -9999.0
    w,=numpy.where(whyflag == 0)
    if w.size > 0:
        a4[w] = rho4[w]/2.0 - 1.0

        s2[w] = (Irr[w]+Icc[w])*(1-a4[w])/(1+a4[w])

        T = Irr + Icc
        w2, = numpy.where(T[w] > 0)
        if w2.size > 0:
            w2 = w[w2]
            e1[w2] = (Icc[w2] - Irr[w2])/T[w2]
            e2[w] = 2.0*Irc[w2]/T[w2]

    for i in xrange(row.size):
        whystr[i] = flagmap[whyflag[i]]

    out={'row':row,
         'col':col,
         'Irr':Irr,
         'Irc':Irc,
         'Icc':Icc,
         'e1':e1,
         'e2':e2,
         'rho4':rho4,
         'a4':a4, 
         's2':s2,
         'uncer':uncer,
         'numiter':numiter,
         'wrow':wrow,
         'wcol':wcol,
         'interp':interpolated,
         'whyflag':whyflag,
         'whystr':whystr,
         'shiftmax':shiftmax,
         'sky':sky,
         'sigsky':sigsky}

    if guess is not None:
        out['guess'] = numpy.array(guess, ndmin=1, copy=False, dtype='f4')
    if is_scalar:
        for key in out:
            if key != 'shiftmax':
                out[key] = out[key][0]
    return out

def get_sky_sigsky(n, sky, sigsky):
    if sky is None:
        sky = numpy.zeros(n, dtype='f4')
    else:
        sky = numpy.array(sky, ndmin=1, copy=False, dtype='f4')
        if sky.size == 1:
            tmpsky = sky[0]
            sky = numpy.zeros(n, dtype='f4')
            sky[:] = tmpsky
        elif sky.size != n:
            raise ValueError("sky must be scalar or same size as row/col")

    if sigsky is None:
        sigsky = numpy.ones(n, dtype='f4')
    else:
        sigsky = numpy.array(sigsky, ndmin=1, copy=False, dtype='f4')
        if sigsky.size == 1:
            tmpsigsky = sigsky[0]
            sigsky = numpy.zeros(n, dtype='f4')
            sigsky[:] = tmpsigsky
        elif sigsky.size != n:
            raise ValueError("sigsky must be scalar or same size as row/col")

    return sky, sigsky


def test(e=0.4, sigma=10.0, ntest=1000, doplot=False):
    """
    Test accuracy with random orientations
    """
    import images

    T = 2*sigma**2
    emeas=numpy.empty(ntest, dtype='f4')
    Tmeas=numpy.empty(ntest, dtype='f4')
    angle=numpy.empty(ntest, dtype='f4')

    ntot=0
    i=0
    while i < ntest:
        tangle = numpy.random.random()*360.0

        e1 = e*numpy.cos(2*tangle*numpy.pi/180.0)
        e2 = e*numpy.sin(2*tangle*numpy.pi/180.0)

        Irr, Irc, Icc = util.ellip2mom(T, e1, e2)

        res, image = test1(Irr,Irc,Icc)

        if res['whyflag'] != 0:
            #raise ValueError("unexpected failure: whystr = %s" % res['whystr'])
            print("angle:",angle,"whystr:",res['whystr'])
            #images.view(image)
            #tmp=raw_input('hit a key: ')
            #if tmp == 'q': return
        else:
            e1m=res['e1']
            e2m=res['e2']
            emeas[i] = numpy.sqrt(e1m**2 + e2m**2)
            Tmeas[i] = res['Irr'] + res['Icc']
            angle[i] = tangle
            i += 1

        ntot += 1

    print("ntest:",ntest,"total trials:",ntot)

    efdiff = emeas/e-1.0
    efdiff_mean = efdiff.mean()
    efdiff_std  = efdiff.std()
    stdout.write("<(e-etrue)/etrue> = %s +/- %s\n" % (efdiff_mean,efdiff_std))

    Tfdiff = Tmeas/T-1.0
    Tfdiff_mean = Tfdiff.mean()
    Tfdiff_std  = Tfdiff.std()
    stdout.write("<(T-Ttrue)/Ttrue> = %s +/- %s\n" % (Tfdiff_mean,Tfdiff_std))


    if doplot:
        import biggles
        import esutil
        tab=biggles.Table(2,1)

        eplt = biggles.FramedPlot()

        be = esutil.stat.Binner(angle % 90, efdiff)
        be.dohist(nperbin=100)
        be.calc_stats()

        #ep=biggles.Points(angle, efdiff, type='filled circle')
        ep=biggles.Points(be['xmean'], be['ymean'], type='filled circle')

        eplt.add(ep)
        eplt.xlabel='angle'
        eplt.ylabel=r'$(e-e_{true})/e_{true}$'

        Tplt = biggles.FramedPlot()

        bT = esutil.stat.Binner(angle % 90, Tfdiff)
        bT.dohist(nperbin=100)
        bT.calc_stats()

        Tp=biggles.Points(bT['xmean'], bT['ymean'], type='filled circle')

        Tplt.add(Tp)
        Tplt.xlabel='angle'
        Tplt.ylabel=r'$(T-T_{true})/T_{true}$'

        tab[0,0] = eplt
        tab[1,0] = Tplt

        tab.show()

def test_admom_residuals(image, cen, guess, sky=0.0, sigsky=1.0):
    """
    Fit adaptive moments to the input image and display the residuals
    """
    import images
    import fimage
    import biggles
    res = admom(image, cen[0], cen[1], guess=guess, sky=sky, sigsky=sigsky,
                nsub=8)

    counts = image.sum()
    wcen=[res['wrow'],res['wcol']]
    fake = fimage.model_image('gauss',image.shape,wcen,res['Irr'],res['Irc'],res['Icc'],
                              counts=counts)

    resid = fake-image

    maxval = max( image.max(), fake.max() )
    minval = 0.0

    levels=7
    tab=biggles.Table(2,3)
    #tab=biggles.Table(3,2)
    implt=images.view(image, levels=levels, show=False, min=minval, max=maxval)
    fakeplt=images.view(fake, levels=levels, show=False, min=minval, max=maxval)
    residplt=images.view(resid, show=False, min=minval, max=maxval)

    sigma = numpy.sqrt((res['Irr']+res['Icc'])/2.0)
    lab = biggles.PlotLabel(0.1,0.9,r'$\sigma$: %0.2f' % sigma, fontsize=4, halign='left')
    fakeplt.add(lab)

    implt.title='original'
    fakeplt.title='gaussian model'
    residplt.title='residuals'


    # cross-sections
    imrows = image[:,wcen[1]]
    imcols = image[wcen[0],:]
    fakerows = fake[:,wcen[1]]
    fakecols = fake[wcen[0],:]
    resrows = resid[:,wcen[1]]
    rescols = resid[wcen[0],:]

    himrows = biggles.Histogram(imrows, color='blue')
    himcols = biggles.Histogram(imcols, color='blue')
    hfakerows = biggles.Histogram(fakerows, color='orange')
    hfakecols = biggles.Histogram(fakecols, color='orange')
    hresrows = biggles.Histogram(resrows, color='red')
    hrescols = biggles.Histogram(rescols, color='red')

    himrows.label = 'image'
    hfakerows.label = 'model'
    hresrows.label = 'resid'
    key = biggles.PlotKey(0.1,0.9,[himrows,hfakerows,hresrows]) 
    rplt=biggles.FramedPlot()
    rplt.add( himrows, hfakerows, hresrows,key )

    cplt=biggles.FramedPlot()
    cplt.add( himcols, hfakecols, hrescols )

    rplt.aspect_ratio=1
    cplt.aspect_ratio=1


    tab[0,0] = implt
    tab[0,1] = fakeplt
    tab[0,2] = residplt
    tab[1,0] = rplt
    tab[1,1] = cplt

    #tab[0,0] = implt
    #tab[0,1] = fakeplt
    #tab[1,0] = residplt
    #tab[1,1] = rplt
    #tab[2,0] = cplt


    tab.show()


def test_errors(e=0.4, sigma=10.0, theta='random', ntest=1000, counts=5000, doplot=False):
    """
    Test the error estimates.

    REMEMBER, for larger images you need higher counts!

    sigma is in pixels.  sigma=1.07 pix is like 1'' seeing in the SDSS

    """
    import images
    from pprint import pprint
    T = 2*sigma**2

    sigsky = 5.5

    etrue=numpy.zeros(ntest,dtype='f4') + e
    emeas=numpy.zeros(ntest,dtype='f4')

    e1true=numpy.zeros(ntest,dtype='f4')
    e2true=numpy.zeros(ntest,dtype='f4')
    e1meas=numpy.zeros(ntest,dtype='f4')
    e2meas=numpy.zeros(ntest,dtype='f4')
    uncer=numpy.zeros(ntest,dtype='f4')
    whyflag=numpy.zeros(ntest,dtype='i4')

    for i in xrange(ntest):
        if theta == 'random':
            angle = numpy.random.random()*360.0
        else:
            angle = theta

        e1 = e*numpy.cos(2*angle*numpy.pi/180.0)
        e2 = e*numpy.sin(2*angle*numpy.pi/180.0)

        Irr, Irc, Icc = util.ellip2mom(T, e1, e2)

        res, image = test1(Irr,Irc,Icc,counts=counts,sigsky=sigsky)

        e1true[i] = e1
        e2true[i] = e2
        e1meas[i] = res['e1']
        e2meas[i] = res['e2']
        emeas[i] = numpy.sqrt( res['e1']**2 + res['e2']**2)

        uncer[i] = res['uncer']
        whyflag[i] = res['whyflag']

        '''
        if res['whyflag'] != 0:
            pprint(res)
            images.multiview(image)
            return
        '''

    stdout.write('\n')

    w,=numpy.where(whyflag == 0)
    stdout.write("Good measurements: %s/%s\n" % (w.size,ntest))
    if w.size != ntest:
        etrue=etrue[w]
        e1true=e1true[w]
        e2true=e2true[w]
        e1meas=e1meas[w]
        e2meas=e2meas[w]
        emeas=emeas[w]
        uncer=uncer[w]
        whyflag=whyflag[w]

    ediff = emeas-etrue
    efrac = (emeas-etrue)/etrue

    e1diff = e1meas-e1true
    e2diff = e2meas-e2true

    ediff_mean = ediff.mean()
    ediff_std  = ediff.std()

    efrac_mean = efrac.mean()
    efrac_std  = efrac.std()

    e1diff_mean = e1diff.mean()
    e1diff_std  = e1diff.std()
    e2diff_mean = e2diff.mean()
    e2diff_std  = e2diff.std()

    uncer_mean = uncer.mean()
    stdout.write("<e1-e1true> = %s +/- %s\n" % (e1diff_mean,e1diff_std))
    stdout.write("<e2-e2true> = %s +/- %s\n" % (e2diff_mean,e2diff_std))
    stdout.write("predicted uncertainty: %s\n" % uncer_mean)
    stdout.write("\n")
    stdout.write("<e-etrue>   = %s +/- %s\n" % (ediff_mean,ediff_std))
    stdout.write("<(e-etrue)/etrue>   = %s +/- %s\n" % (efrac_mean,efrac_std))

    if doplot:
        import biggles
        rng = [-4.0*uncer_mean, 4.0*uncer_mean]
        nbin=40
        e1h,e1edges = numpy.histogram(e1diff, range=rng, bins=nbin, normed=True)
        e2h,e2edges = numpy.histogram(e2diff, range=rng, bins=nbin, normed=True)

        p1 = biggles.Histogram(e1h, x0=e1edges[0], binsize=e1edges[1]-e1edges[0], color='blue')
        p1.label = r'$e_1-e_1^{true}'
        p2 = biggles.Histogram(e2h, x0=e2edges[0], binsize=e2edges[1]-e2edges[0], color='red')
        p2.label = r'$e_2-e_2^{true}'

        sigma = uncer_mean
        x=numpy.linspace(rng[0],rng[1],nbin)
        gauss = numpy.exp(-0.5*x**2/sigma**2 )
        # use same weak normalization method as numpy.histogram 

        db = numpy.array(numpy.diff(e1edges), float)
        gauss = gauss/(gauss*db).sum()


        pg = biggles.Curve(x, gauss,color='black')
        pg.label = 'Expected'

        key=biggles.PlotKey(0.1,0.9,[p1,p2,pg])

        theta_lab = biggles.PlotLabel(0.95,0.9,'theta=%s' % theta, halign='right')
        e_lab     = biggles.PlotLabel(0.95,0.85,'e=%s' % e, halign='right')
        T_lab     = biggles.PlotLabel(0.95,0.80,'T=%s' % T, halign='right')

        plt=biggles.FramedPlot()

        plt.add(p1,p2,pg,key,theta_lab, e_lab, T_lab)

        plt.xlabel = r'$e-e^{true}$'

        plt.show()


        

def test1(Irr,Irc,Icc,counts=1, sigsky=None, nsub=4):
    import pprint
    import fimage

    T = Irr+Icc
    sigma = numpy.sqrt(T/2)
    imsize = int( 2*4.5*sigma )
    if (imsize % 2) == 0:
        imsize += 1
    imdims=[imsize]*2
    cen = [(imsize-1)/2., (imsize-1)/2.]

    e1 = (Icc-Irr)/(Irr+Icc)
    e2 = 2*Irc/(Irr+Icc)

    g = fimage.model_image('gauss', imdims, cen, Irr, Irc, Icc, counts=counts, nsub=16)

    if sigsky is not None:
        g[:,:] += numpy.random.normal(size=g.size, scale=sigsky).reshape(g.shape)
    res = admom(g, cen[0], cen[1], sigsky=sigsky, guess=T/2, nsub=nsub)

    res['Irr_true'] = Irr
    res['Irc_true'] = Irc
    res['Icc_true'] = Icc
    res['e1_true'] = e1
    res['e2_true'] = e2

    return res, g

def view2(im1,im2):
    import biggles
    import images
    tab = biggles.Table(1,2)

    p1 = images.view(im1,show=False)
    p2 = images.view(im2,show=False)

    tab[0,0] = p1
    tab[0,1] = p2
    tab.show()
