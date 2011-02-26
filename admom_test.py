'''
Testing programs

test_sub_pixel: test sub-pixel effects as a function of sub-pixel correction
    factor resolution, nsub

CovarVsInput: A class for measuring the covariance matrix and comparing to
    input.
run_covar_vs_input: shortcut function

'''
from __future__ import print_function

import numpy
from numpy import linspace,array,zeros,sqrt,log10,sin
from .wrappers import admom

import os
import esutil as eu
from esutil.ostools import path_join
from esutil.numpy_util import where1

from . import util

import time

def run_covar_vs_input(model):
    cvi = CovarVsInput(model)
    cvi.measure_covar()

class CovarVsInput:
    '''
    THERE IS ALSO ONE FOR unweighted moms

    Compare the input covariance matrix to measured using
    adaptive moments for non-gaussian models.
    '''
    def __init__(self, model):
        if model not in ['exp','dev']:
            raise ValueError("models 'exp','dev'")
        self.model = model

        if model == 'exp':
            self.sigfac = 7.0
        else:
            raise ValueError("Figure out dev sigfac")

    def sigma_vals(self):
        return linspace(1.0, 20.0, 20)
    def ellip_vals(self):
        return linspace(0.0,0.7,7+1)

    def struct(self, n):
        data=zeros(n, dtype=[('sigma_index','i4'),
                             ('ellip_index','i4'),
                             ('Irr_input','f4'),
                             ('Irc_input','f4'),
                             ('Icc_input','f4'),
                             ('Irr_meas','f4'),
                             ('Irc_meas','f4'),
                             ('Icc_meas','f4')])
        return data

    def dir(self):
        dir=os.environ.get('REGAUSSIM_DIR',None)
        if dir is None:
            raise ValueError("REGAUSSIM_DIR must be set")
        dir=path_join(dir,'admom-covar-meas-vs-input')

        if not os.path.exists(dir):
            os.makedirs(dir)
        return dir

    def read(self):
        return eu.io.read(self.file())

    def file(self):
        dir=self.dir()
        f='covar-meas-vs-input-%s.rec' % self.model
        f=path_join(dir,f)
        return f

    def epsfile(self, type):
        dir=self.dir()
        epsfile = '%s-%s.eps' % (type,self.model)
        return path_join(dir,epsfile)


    def measure_covar(self):
        """
        
        Test the measured covariance matrix vs the input

        """

        import fimage

        f=self.file()

        sigma_vals = self.sigma_vals()
        ellip_vals = self.ellip_vals()
        data = self.struct(sigma_vals.size*ellip_vals.size)

        ii=0
        for i in xrange(sigma_vals.size):
            sigma=sigma_vals[i]
            for j in xrange(ellip_vals.size):
                ellip = ellip_vals[j]

                Irr,Irc,Icc = util.ellip2mom(2*sigma**2, e=ellip, theta=0.0)

                dim = int( numpy.ceil(2.*self.sigfac*sigma ) )
                if (dim % 2) == 0:
                    dim += 1
                dims=[dim,dim]
                cen=[(dim-1)/2]*2
                print("sigma:",sigma,"ellip:",ellip,"dims:",dims)

                im=fimage.model_image(self.model,dims,cen,[Irr,Irc,Icc],nsub=8)
                res = admom(im, cen[0], cen[1], guess=sigma, nsub=4)

                data['sigma_index'][ii] = i
                data['ellip_index'][ii] = j
                data['Irr_input'][ii] = Irr
                data['Irc_input'][ii] = Irc
                data['Icc_input'][ii] = Icc
                data['Irr_meas'][ii] = res['Irr']
                data['Irc_meas'][ii] = res['Irc']
                data['Icc_meas'][ii] = res['Icc']

                ii+=1

        hdr={'model':self.model}
        eu.io.write(f, data, delim=' ', verbose=True, clobber=True, header=hdr)

    def plot_ellip_vs_input(self, show=False):
        '''
        Plot the measured ellip as a function of the input for sigma_index=0
        which is a reasonably large object
        '''
        import biggles
        from biggles import PlotLabel,FramedPlot,Table,Curve,PlotKey,Points
        from pcolors import rainbow
        import pprint

        data = self.read()
        w=where1(data['sigma_index'] == 10)
        data = data[w]
        
        e1_input, e2_input, Tinput = util.mom2ellip(data['Irr_input'],
                                                    data['Irc_input'],
                                                    data['Icc_input'])
        e1_meas, e2_meas, Tinput = util.mom2ellip(data['Irr_meas'],
                                                  data['Irc_meas'],
                                                  data['Icc_meas'])

        einput = sqrt(e1_input**2 + e2_input**2)
        emeas = sqrt(e1_meas**2 + e2_meas**2)

        plt=FramedPlot()

        p = Points(einput,emeas, type='filled circle')
        plt.add(p)
        plt.xlabel=r'$\epsilon_{in}$'
        plt.ylabel=r'$\epsilon_{meas}$'

        sig=sqrt((data['Irr_meas'][0]+data['Icc_meas'][0])/2)
        lab1=PlotLabel(0.1,0.9,self.model, halign='left')
        lab2=PlotLabel(0.1,0.8,r'$\sigma: %0.2f$' % sig, halign='left')
        plt.add(lab1,lab2)


        einput.sort()
        c = Curve(einput, einput, color='red')
        c.label = r'$\epsilon_{input} = \epsilon_{meas}$'
        key=PlotKey(0.95,0.07,[c], halign='right')
        plt.add(c)
        plt.add(key)
        if show:
            plt.show()

        epsfile=self.epsfile('ellip-vs-input')
        print("Writing eps file:",epsfile)
        plt.write_eps(epsfile)



    def plot_size_vs_input(self, show=False):
        '''
        Plot recovered size vs input for ellip=0, which is ellip_index=0
        '''

        import biggles
        from biggles import PlotLabel,FramedPlot,Table,Curve,PlotKey,Points
        from pcolors import rainbow
        import pprint

        data = self.read()
        w=where1(data['ellip_index'] == 0)
        data = data[w]

        siginput = sqrt(data['Irr_input'])
        sigmeas  = sqrt(data['Irr_meas'])


        pars=numpy.polyfit(siginput, sigmeas, 1)
        print("offset:",pars[1])
        print("slope: ",pars[0])
        print("IGNORING OFFSET")

        plt=FramedPlot()

        p = Points(siginput,sigmeas, type='filled circle')
        plt.add(p)
        plt.xlabel=r'$\sigma_{in}$'
        plt.ylabel=r'$\sigma_{meas}$'

        lab=PlotLabel(0.1,0.9,self.model)
        plt.add(lab)

        yfit2=pars[0]*siginput
        cfit2=Curve(siginput, yfit2, color='steel blue')
        cfit2.label = r'$%0.2f \sigma_{in}$' % pars[0]

        plt.add( cfit2 )

        key=PlotKey(0.95,0.07,[cfit2], halign='right')
        plt.add(key)
        if show:
            plt.show()

        epsfile=self.epsfile('size-vs-input')
        print("Writing eps file:",epsfile)
        plt.write_eps(epsfile)




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

