from __future__ import print_function

import numpy
from numpy import cos, sin, pi as PI

import admom

try:
    import fimage
except:
    print("Could not import fimage")
try:
    import images
except:
    print("Could not import images module")
try:
    import scipy.signal
except:
    print("Could not import scipy.signal")

from . import unweighted
from pprint import pprint

class ReGauss(dict):
    def __init__(self, image, row, col, psf, **keys):
        '''
        The psf must have odd dimention in pixels, e.g. [31,31] and the center must
        be at [(nrow-1)/2, (ncol-1)/2]

        The image sky value must be zero.
        The psf sky value must be zero.
        '''

        self.image = image
        self.row=row
        self.col=col
        self.psf = psf

        self.conv = keys.get('conv','fft')
        self.sigsky = keys.get('sigsky', 1.0)
        self.guess=keys.get('guess', None)
        self.guess_psf=keys.get('guess_psf',None)
    
        self.detf0_tol = keys.get('detf0_tol',1.e-5)

        self.verbose=keys.get('verbose',False)
        self.debug=keys.get('debug',False)

        self.check_init()

    def check_init(self):
        if self.conv not in ['fft','real']:
            raise ValueError("conv should be one of: "+str(['fft','real']))

        if len(self.image.shape) != 2:
            raise ValueError("image must be 2-d")
        if len(self.psf.shape) != 2:
            raise ValueError("psf must be a 2-d image")

        if (self.psf.shape[0] % 2) != 1 or (self.psf.shape[1] % 2) != 1:
            raise ValueError("psf must have odd number of pixels in each "
                             "dimension for FFT convolution")

        if (not numpy.isscalar(self.row) or not numpy.isscalar(self.col)
                or not numpy.isscalar(self.sigsky)):
            raise ValueError("regauss requires scalar row,col,sigsky")



    def do_all(self):
        self.do_admom()
        self.do_psf_admom()
        self.do_basic_corr()
        self.do_regauss()
        self.do_rg_corr()
        self.do_unweighted()

    def do_unweighted(self):
        res = unweighted.correct(self.image, self.psf)
        self['uwcorrstats'] = res

    def do_admom(self):
        '''
        row, col:
            row,column in image.  Note numpy arrays are indexed [row,col]
        sky:
            sky value at the location of row,col
        sigsky:
            sky variance at row,col
        '''


        out = admom.admom(self.image, 
                          self.row, 
                          self.col, 
                          sky=0.0,
                          sigsky=self.sigsky,
                          guess=self.guess)

        self['imstats'] = out
        if self.verbose:
            print("image stats:")
            pprint(out)


    def do_psf_admom(self):
        row = (self.psf.shape[0]-1)/2
        col = (self.psf.shape[1]-1)/2

        out = admom.admom(self.psf, 
                          row, 
                          col, 
                          sky=0, 
                          guess=self.guess_psf)

        self['psfstats'] = out
        if self.verbose:
            print("psf stats:")
            pprint(out)

    def do_basic_corr(self):
        if 'imstats' not in self or 'psfstats' not in self:
            raise ValueError("run do_psf_admom() and do_admom() first")

        if self['psfstats']['whyflag'] != 0 or self['imstats']['whyflag'] != 0:
            self['corrstats'] = None
            return

        # apply basic correction
        ims = self['imstats']
        psfs = self['psfstats']
        e1,e2,R,flags = admom.correct(ims['Irr']+ims['Icc'],
                                      ims['e1'],ims['e2'],ims['a4'],
                                      psfs['Irr']+psfs['Icc'],
                                      psfs['e1'],psfs['e2'],psfs['a4'])

        self['corrstats'] = {'e1':e1,'e2':e2,'R':R,'flags':flags}
        if self.verbose:
            print("corrstats:")
            pprint(self['corrstats'])


    def do_regauss(self):
        self['rgstats'] = None

        if 'imstats' not in self or 'psfstats' not in self:
            raise ValueError("run admom on image and psf first")

        if self['imstats']['whyflag'] != 0:
            if self.verbose:
                print("admom failed, cannot run regauss")
            return

        self.make_f0()
        if self.f0 == None:
            self['rgstats'] = {'f0flags':2**0}
            return

        self.make_epsilon()
        self.make_f0conv()

        self.iprime = self.image - self.f0conv

        if self.debug:
            images.compare_images(self.image, self.iprime, label1='original',label2='minus f0conv')
            #k=raw_input('hit a key: ')

        guess = (self['imstats']['Irr'] + self['imstats']['Irr'])/2
        wrow = self['imstats']['wrow']
        wcol = self['imstats']['wcol']
        out = admom.admom(self.iprime,
                          wrow,
                          wcol,
                          sky=0.0,
                          sigsky=self.sigsky,
                          guess=guess)

        out['f0flags'] = 0
        self['rgstats'] = out

        if self.verbose:
            print("rgstats:")
            pprint(self['rgstats'])
   
    def make_f0(self):
        if 'imstats' not in self or 'psfstats' not in self:
            raise ValueError("run admom on image and psf first")

        Irr=self['imstats']['Irr']
        Irc=self['imstats']['Irc']
        Icc=self['imstats']['Icc']

        Irr_psf=self['psfstats']['Irr']
        Irc_psf=self['psfstats']['Irc']
        Icc_psf=self['psfstats']['Icc']

        # construct the f0 parameters
        Irr_f0 = Irr - Irr_psf
        Irc_f0 = Irc - Irc_psf
        Icc_f0 = Icc - Icc_psf
        det_f0 = Irr_f0*Icc_f0 - Irc_f0**2

        imcounts = self.image.sum()

        self.f0 = None
        if det_f0 > self.detf0_tol:
            wrow = self['imstats']['wrow']
            wcol = self['imstats']['wcol']
            self.f0 = fimage.model_image('gauss',self.image.shape,[wrow,wcol],
                                         [Irr_f0, Irc_f0, Icc_f0],
                                         counts=imcounts)

            if self.debug:
                plt=images.multiview(self.f0,show=False,levels=7)
                plt.title='f0'
                plt.show()
                #k=raw_input('hit a key: ')
        else:
            if self.verbose:
                print("Found det(f0) less than tolerance:",det_f0,"<",self.detf0_tol)

    def make_epsilon(self):
        """
        make a model image for the fit gaussian and subtract it from the
        psf.
        
        This becomes a convolution kernel on our simplified model for the galaxy
        
        Note psf and the subtracted gaussian are both normalized to 1, and
        epsilon thus integrates to zero
        """


        if 'psfstats' not in self:
            raise ValueError("run admom on psf first")

        Irr = self['psfstats']['Irr']
        Irc = self['psfstats']['Irc']
        Icc = self['psfstats']['Icc']

        row = (self.psf.shape[0]-1)/2
        col = (self.psf.shape[1]-1)/2

        gauss = fimage.model_image('gauss',self.psf.shape,[row,col],
                                   [Irr,Irc,Icc],counts=1)
        
        # need both our gaussian and the psf to be normalized
        tpsf = self.psf/self.psf.sum()

        if self.debug:
            images.compare_images(tpsf, gauss, label1='psf',label2='gauss')
            #k=raw_input('hit a key: ')

        epsilon = tpsf - gauss

        self.epsilon = epsilon


    def make_f0conv(self):
        """
        Make f0 convolved with epsilon.  It is on this image we will
        measure new adaptive moments.

        f0 will be expanded if epsilon is larger.  The center will remain
        the same in that expansion

        """
        if self.f0 is None or self.epsilon is None:
            raise ValueError("Create both f0 and epsilon first")
        if self.f0 is None:
            raise ValueError("f0 is None")

        f0 = self.f0
        ep = self.epsilon

        if self.verbose:
            if ep.shape[0] > f0.shape[0] or ep.shape[1] > f0.shape[0]:
                print("epsilon is larger than image:")
                print("  epsize:",ep.shape)
                print("  f0size:",f0.shape)
        # if f0 is bigger, it is returned unchanged
        f0_expand = images.expand(f0, ep.shape)

        if self.conv == 'fft':
            f0conv = scipy.signal.fftconvolve(f0_expand, ep, mode='same')
        else:
            f0conv = scipy.signal.convolve2d(f0_expand, ep, old_behavior=False,mode='same')

        # trim back to the original size
        if (f0conv.shape[0] > self.image.shape[0] 
                or f0conv.shape[1] > self.image.shape[1]):
            f0conv = f0conv[ 0:self.image.shape[0], 0:self.image.shape[1] ]
            if self.debug:
                print("trimming back")

        self.f0conv = f0conv

        if self.debug:
            print("f0.shape",self.f0.shape)
            print("f0conv.shape",self.f0conv.shape)
            plt=images.view(self.f0,show=False,levels=7)
            plt.title='f0'
            plt.show()
            cplt=images.view(self.f0conv, show=False,levels=7)
            cplt.title='f0conv'
            cplt.show()
            #images.compare_images(self.f0, self.f0conv, label1='f0',label2='f0conv')
            #k=raw_input('hit a key: ')


    def do_rg_corr(self):
        self['rgcorrstats'] = None
        if 'rgstats' not in self or 'psfstats' not in self:
            return
        if self['psfstats'] is None or self['rgstats'] is None:
            return
        if len(self['rgstats']) == 1:
            # this is when we couldn't make f0
            return

        if self['psfstats']['whyflag'] != 0 or self['rgstats']['whyflag'] != 0:
            return

        # apply basic correction
        rgs = self['rgstats']
        psfs = self['psfstats']

        #e1,e2,R,flags = admom.correct(rgs['Irr']+rgs['Icc'],
        #                              rgs['e1'],rgs['e2'],rgs['a4'],
        #                              psfs['Irr']+psfs['Icc'],
        #                              psfs['e1'],psfs['e2'],psfs['a4'])

        # the PSF is supposed to be gaussian now, put in a4=0?
        e1,e2,R,flags = admom.correct(rgs['Irr']+rgs['Icc'],
                                      rgs['e1'],rgs['e2'],rgs['a4'],
                                      psfs['Irr']+psfs['Icc'],
                                      psfs['e1'],psfs['e2'],0.0)

        self['rgcorrstats'] = {'e1':e1,'e2':e2,'R':R,'flags':flags}
        if self.verbose:
            print("rgcorrstats:")
            pprint(self['rgcorrstats'])



def test_regauss(gal_model='gauss',
                 gal_T=16.0,
                 gal_e1=0.0,
                 gal_e2=0.0,
                 psf_T=4.0,
                 psf_e1=0.0,
                 psf_e2=0.0,
                 show=True):

    from pprint import pprint
    import fimage

    nsigma=5

    psf_sigma=numpy.sqrt(psf_T/2.)
    psf_imsize = int( round(2*nsigma*psf_sigma) )
    if (psf_imsize % 2) == 0:
        psf_imsize+=1

    psf_dims=[psf_imsize]*2
    psf_cen=[(psf_imsize-1)/2.]*2
    psf_moms = admom.ellip2mom(psf_T,e1=psf_e1,e2=psf_e2)
    psf = fimage.model_image('gauss', 
                             psf_dims,
                             psf_cen,
                             psf_moms,
                             nsub=16)

    gal_moms = admom.ellip2mom(gal_T,e1=gal_e1,e2=gal_e2)


    # for the image size, make it big enough to easily fit
    # the *convolved* galaxy
    gal_convolved_sigma = admom.mom2sigma(gal_T+psf_T)
    gal_imsize = int( round(2*nsigma*gal_convolved_sigma) )
    if gal_imsize < psf_imsize:
        gal_imsize = psf_imsize

    gal_cen  = [ (gal_imsize-1)/2. ]*2
    gal_dims = [int(gal_imsize)]*2

    gal = fimage.model_image(gal_model,
                             gal_dims,
                             gal_cen,
                             gal_moms,
                             nsub=16)

    levels=7

    # now convolve with the PSF
    imconv = scipy.signal.fftconvolve(gal, psf, mode='same')

    if show:
        images.multiview(psf, levels=levels)
        images.multiview(gal, levels=levels)
        images.multiview(imconv, levels=levels)

    rg = ReGauss(imconv, gal_cen[0], gal_cen[1], psf, guess=gal_T/2, guess_psf=psf_T/2)

    rg.do_admom()
    print("admom stats")
    pprint(rg['imstats'])
    rg.do_psf_admom()
    print("psf admom stats")
    pprint(rg['psfstats'])
    rg.do_basic_corr()
    print("admom corrected stats")
    pprint(rg['corrstats'])
    rg.do_regauss()
    print("rg stats")
    pprint(rg['rgstats'])
    rg.do_rg_corr()
    print("rg corrected stats")
    pprint(rg['rgcorrstats'])
