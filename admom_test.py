from __future__ import print_function

import numpy
from numpy import linspace,array,zeros,sqrt,log10,sin
from .wrappers import admom

import os
import esutil as eu
from esutil.ostools import path_join

from . import util

import time

def size_meas_vs_size_input(model):
    """
    
    Test the measured Ixx+Iyy compared to the input Ixx+Iyy in the covariance
    matrix.  

    """

    import fimage
    import images

    if model not in ['exp','dev']:
        raise ValueError("models 'exp','dev'")
    if model == 'exp':
        sigfac = 7.0
    else:
        raise ValueError("Figure out dev sigfac")

    sigma_vals = linspace(1.0, 20.0, 20)

    data=zeros(sigma_vals.size, dtype=[('sigma_input','f4'),
                                       ('sigma_meas','f4')])

    for i in xrange(sigma_vals.size):
        sigma = sigma_vals[i]

        dim = int( numpy.ceil(2.*sigfac*sigma ) )
        if (dim % 2) == 0:
            dim += 1
        dims=[dim,dim]
        cen=[(dim-1)/2]*2
        print("sigma:",sigma,"dims:",dims)

        print("Creating image...")
        im=fimage.model_image(model,dims,cen,[sigma**2,0.0,sigma**2],nsub=8)

        print("Running admom...")
        res = admom(im, cen[0], cen[1], guess=sigma, nsub=4)

        data['sigma_input'][i] = sigma
        data['sigma_meas'][i] = sqrt(res['Irr'])

        implt = images.multiview(im, levels=7, show=False)
        epsfile = sizemeas_sigma_image_file(model, sigma)
        print("Writing image file:",epsfile)
        implt.write_eps(epsfile)

    f=sizemeas_file(model)
    eu.io.write(f, data, verbose=True, clobber=True)

'''
def fitfunc(p,x):
    return p[0] + p[1]*x
def errfunc(p,x,y):
    return fitfunc(p,x)-y
'''

def plot_size_meas_vs_size_input(model,show=False):
    #from scipy import optimize

    import biggles
    from biggles import PlotLabel,FramedPlot,Table,Curve,PlotKey,Points
    from pcolors import rainbow
    import pprint

    data = eu.io.read(sizemeas_file(model))

    # Target function

    fitfunc = lambda p, x: (p[0] + p[1]*x)
    errfunc = lambda p, x, y: (fitfunc(p, x) - y) # Distance to the target function

    x,y = data['sigma_input'].copy(),data['sigma_meas'].copy()

    pars=numpy.polyfit(x,y,1)

    plt=FramedPlot()

    p = Points(x,y, type='filled circle')
    plt.add(p)
    plt.xlabel=r'$\sigma$ input'
    plt.ylabel=r'$\sigma$ measured'

    lab=PlotLabel(0.1,0.9,model)
    plt.add(lab)


    yfit=pars[0]*x + pars[1]
    cfit=Curve(data['sigma_input'], yfit, color='red')
    cfit.label = r'$%0.2f \sigma_{in} + %0.2f$' % tuple(pars)
    plt.add( cfit )
    key=PlotKey(0.95,0.07,[cfit], halign='right')
    plt.add(key)
    if show:
        plt.show()

    epsfile=sizemeas_file(model, 'eps')
    print("Writing eps file:",epsfile)
    plt.write_eps(epsfile)

def sizemeas_outdir():
    dir=os.environ.get('REGAUSSIM_DIR',None)
    if dir is None:
        raise ValueError("REGAUSSIM_DIR must be set")
    dir=path_join(dir,'admom-size-meas-vs-input')
    return dir

def sizemeas_file(model, type='fits'):
    dir=sizemeas_outdir()
    f='size-meas-vs-input-%s' % model
    f=path_join(dir,f+'.'+type)
    return f
def sizemeas_sigma_image_file(model, sigma):
    dir=sizemeas_outdir()
    f='size-meas-vs-input-%s-sigma%05.2f.eps' % (model,sigma)
    f=path_join(dir,f)
    return f






def test_sub_pixel_many():
    ellipvals = linspace(0.0,0.7,7+1)
    thetavals=[0.,15.,30.,45.]

    for e in ellipvals:
        for theta in thetavals:
            print('-'*70)
            test_sub_pixel(e,theta)

def plot_sub_pixel_many():
    ellipvals = linspace(0.0,0.7,7+1)
    thetavals=[0.,15.,30.,45.]

    for e in ellipvals:
        for theta in thetavals:
            print('-'*70)
            plot_sub_pixel(e,theta)


def test_sub_pixel(ellip,theta):
    """
    Round objects for now
    """
    import fimage
    import images

    sigma_vals = linspace(1.0,3.0,10)
    nsub_vals = array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],dtype='i4')

    data = subpixel_struct(sigma_vals.size, nsub_vals.size)

    ct=0
    amt=0
    for j in xrange(sigma_vals.size):
        sigma = sigma_vals[j]
        data['sigma'][j] = sigma


        Irr,Irc,Icc=util.ellip2mom(2*sigma**2, e=ellip, theta=theta)
        print("sigma: %0.2f Irr: %0.2f Irc: %0.2f Icc: %0.2f" % (sigma,Irr,Irc,Icc))

        d=int( numpy.ceil( 4.5*sigma*2 ) )
        if (d % 2) == 0:
            d += 1
        cen=[(d-1)/2]*2
        dims = [d,d]

        # create image with super good subpixel integration 
        ct0=time.time()
        im=fimage.model_image('gauss',dims,cen,[Irr,Irc,Icc],nsub=64)
        ct += time.time()-ct0

        amt0=time.time()
        for i in xrange(nsub_vals.size):
            nsub = nsub_vals[i]
            res = admom(im, cen[0], cen[1], guess=sigma**2, nsub=nsub)
            if res['whyflag'] != 0:
                print("    ** Failure:'%s'" % res['whystr'])
                images.multiview(im, levels=7)
                #key=raw_input('hit a key: ')
                #if key == 'q': return

            data['Irr'][j,i] = res['Irr']
            data['Irc'][j,i] = res['Irc']
            data['Icc'][j,i] = res['Icc']
            data['a4'][j,i] = res['a4']
            data['nsub'][j,i] = nsub

            #print("  nsub:",nsub)

        amt += time.time()-amt0

    #print("time for creation:",ct)
    #print("time for admom:",amt)

    outf=subpixel_file(ellip,theta,'fits')
    eu.io.write(outf,data,verbose=True,clobber=True)
    plot_sub_pixel(ellip,theta)

def plot_sub_pixel(ellip,theta, show=False):
    import biggles
    from biggles import PlotLabel,FramedPlot,Table,Curve,PlotKey,Points
    from pcolors import rainbow

    f=subpixel_file(ellip,theta,'fits')
    data = eu.io.read(f)
    colors = rainbow(data.size,'hex')

    pltSigma = FramedPlot()
    pltSigma.ylog=1
    pltSigma.xlog=1

    curves=[]
    for j in xrange(data.size):
        sigest2 = (data['Irr'][j,:] + data['Icc'][j,:])/2

        pdiff = sigest2/data['sigma'][j]**2 -1
        nsub=numpy.array(data['nsub'][j,:])

        #pc = biggles.Curve(nsub, pdiff, color=colors[j])
        pp = Points(data['nsub'][j,:], pdiff, type='filled circle',color=colors[j])

        pp.label = r'$\sigma: %0.2f$' % data['sigma'][j]
        curves.append(pp)
        pltSigma.add(pp)
        #pltSigma.add(pc)
        #pltSigma.yrange=[0.8,1.8]
        #pltSigma.add(pp)


    c5 = Curve(linspace(1,8, 20), .005+zeros(20))
    pltSigma.add(c5)

    key=PlotKey(0.95,0.95,curves,halign='right',fontsize=1.7)
    key.key_vsep=1

    pltSigma.add(key)
    pltSigma.xlabel='N_{sub}'
    pltSigma.ylabel=r'$\sigma_{est}^2  /\sigma_{True}^2 - 1$'

    lab=PlotLabel(0.05,0.07,r'$\epsilon: %0.2f \theta: %0.2f$' % (ellip,theta),halign='left')
    pltSigma.add(lab)

    pltSigma.yrange = [1.e-5,0.1]
    pltSigma.xrange = [0.8,20]
    if show:
        pltSigma.show()

    epsfile=subpixel_file(ellip,theta,'eps')

    print("Writing eps file:",epsfile)
    pltSigma.write_eps(epsfile)

def subpixel_outdir():
    dir=os.environ.get('REGAUSSIM_DIR',None)
    if dir is None:
        raise ValueError("REGAUSSIM_DIR must be set")
    dir=path_join(dir,'admom-subpixel')
    return dir

def subpixel_file(ellip, theta, type='fits'):
    dir=subpixel_outdir()
    f='subpixel-accuracy-e%0.2f-theta%0.1f' % (ellip,theta)
    f=path_join(dir,f+'.'+type)
    return f

def subpixel_struct(nsigma, n_nsub):
    dt=[('sigma','f4'),
        ('nsub','f4',n_nsub),
        ('Irr','f4',n_nsub),
        ('Irc','f4',n_nsub),
        ('Icc','f4',n_nsub),
        ('a4','f4',n_nsub)]
    return zeros(nsigma, dtype=dt)


def testfit():
    import biggles 
    from biggles import FramedPlot,Points,Curve
    import scipy
    from scipy.optimize import leastsq

    ## Parametric function: 'v' is the parameter vector, 'x' the independent varible
    fp = lambda v, x: v[0]/(x**v[1])*sin(v[2]*x)

    ## Noisy function (used to generate data to fit)
    v_real = [1.5, 0.1, 2.]
    fn = lambda x: fp(v_real, x)

    ## Error function
    e = lambda v, x, y: (fp(v,x)-y)

    ## Generating noisy data to fit
    n = 30
    xmin = 0.1
    xmax = 5
    x = linspace(xmin,xmax,n)
    y = fn(x) + scipy.rand(len(x))*0.2*(fn(x).max()-fn(x).min())

    ## Initial parameter value
    v0 = [3., 1, 4.]

    ## Fitting
    v, success = leastsq(e, v0, args=(x,y), maxfev=10000)

    print('Estimater parameters: ', v)
    print('Real parameters: ', v_real)
    X = linspace(xmin,xmax,n*5)
    plt=FramedPlot()
    plt.add(Points(x,y))
    plt.add(Curve(X,fp(v,X),color='red'))

    plt.show()

