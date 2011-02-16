import numpy
from numpy import where, sqrt, arctanh, exp, log

BADE      = 0x1  # 2^0 Failed first sqrt test in CompEA4
BADEMULT1 = 0x2  # 2^1 Failed in first shearmult
BADEMULT2 = 0x4  # 2^2 Failed in second shearmult
BADERED1  = 0x8  # 2^3 Failed sqrt e1red test
BADERED2  = 0x10 # 2^4 Failed sqrt e1red test
BADRVAL   = 0x20 # 2^5 R==0.0



def correct(To, e1o, e2o, a4o,
            Tp, e1p, e2p, a4p):
    """
    Inputs:
        To: Irr + Icc for the object
        e1o, e2o: ellipticity for object
        a4o: rho4/2-1 for the object
        Tp, e1p, e2p, a4p: same for the psf
    Outputs:
        A tuple:
            (e1, e2, R)
        e1: psf-corrected e1
        e2: psf-corrected e2
        R: Resolution parameter.
        flags: processing flags
    """

    isscalar = numpy.isscalar(To)

    try:
        To,e1o,e2o,a4o=_check4arr(To,e1o,e2o,a4o)
        Tp,e1p,e2p,a4p=_check4arr(Tp,e1p,e2p,a4p)
    except:
        raise ValueError("All entries must be the same length")

    if (To.size != Tp.size):
        raise ValueError("All entries must be the same length")

    flags = numpy.zeros(To.size, dtype='i4')
    e1 = numpy.empty(To.size, dtype='f8')
    e2 = numpy.empty(To.size, dtype='f8')
    R = numpy.empty(To.size, dtype='f8')
    e1[:] = -9999.0
    e2[:] = -9999.0
    R[:] = -9999.0

    for i in xrange(e1o.size):
        te1, te2, tR, tflags = _compea41(To[i],e1o[i],e2o[i],a4o[i],
                                         Tp[i],e1p[i],e2p[i],a4p[i])
        e1[i] = te1
        e2[i] = te2
        R[i] = tR
        flags[i] = tflags

    if isscalar:
        e1=e1[0]
        e2=e2[0]
        R=R[0]
        flags=flags[0]

    return e1,e2,R,flags


def _compea41(To,e1o,e2o,a4o,
              Tp,e1p,e2p,a4p):

    e1=-9999.
    e2=-9999.
    R=-9999.
    flags=0

    Tratio = Tp/To

    # Take us to sig2ratio = sigma2(P)/sigma2(O) since this is
    # shear-invariant */
    tmp1 = 1-e1p**2-e2p**2
    tmp2 = 1-e1o**2-e2o**2

    if tmp1 <= 0.0 or tmp2 <= 0.0:
        flags = BADE
        return e1, e2, R, flags

    coshetap = 1./sqrt(tmp1)
    coshetao = 1./sqrt(tmp2)
    sig2ratio = Tratio * coshetao/coshetap # since sigma2 = T / cosh eta
    e1red,e2red=_shearmult1(e1o,e2o,-e1p,-e2p)
    if e1red < -9990:
        flags = BADEMULT1
        return e1, e2, R, flags

    # compute resolution factor and un-dilute

    tmp1 = e1red*e1red+e2red*e2red
    if tmp1 <= 0.0:
        flags = BADERED1
        return e1, e2, R, flags

    e = sqrt(tmp1)
    eta = arctanh(e)
    a2 = exp(-eta)*sig2ratio  # fraction of major axis variance from PSF
    b2 = exp(eta)*sig2ratio   # fraction of minor axis variance from PSF

    A = 1-a2; B = 1-b2  # fractions from intrinsic image
    ca4p = 0.375*(a2**2+b2**2)+0.25*a2*b2
    ca4i = 0.375*(A**2+B**2)+0.25*A*B
    a4i = (a4o - ca4p*a4p) / ca4i
    Ti = (A-B) * (-2+1.5*(A+B))
    Tpp = (a2-b2) * (-2+1.5*(a2+b2))
    deltaeta = Ti * a4i + Tpp * a4p

    # 4th moment correction for R: must find etai
    EI = e;
    etai = 0.5 * log( (1./a2-1) / (1./b2-1) )

    tmp2 = 1-e1red**2-e2red**2

    if tmp2 <= 0.0:
        flags = BADERED2
        return e1, e2, R, flags

    coshetao = 1./sqrt(tmp2)
    deltamu = (-1.5*A**2 - A*B - 1.5*B**2 +2*(A+B)) * a4i + (-1.5*a2**2 - a2*b2 - 1.5*b2**2 + 2*(a2+b2))*a4p
    deltamu *= 0.5
    deltaeta *= -1.0

    # This is equation B16 Hirata & Seljak */
    Rtmp = ( 1 - 2*deltamu - deltaeta*EI - sig2ratio/coshetao ) / (-deltaeta/EI + 1-2*deltamu )

    if Rtmp == 0.0:
        flags = BADRVAL
        return e1, e2, R, flags

    e1red /= Rtmp
    e2red /= Rtmp

    e1,e2 = _shearmult1(e1red,e2red,e1p,e2p)
    if e1 < -9990:
        flags = BADEMULT2
        return e1, e2, R, flags

    R = Rtmp
    return e1, e2, R, flags

def shearmult(e1a, e2a, e1b, e2b):

    isscalar = numpy.isscalar(e1a)
    e1a,e2a,e1b,e2b=_check4arr(e1a,e2a,e1b,e2b)
    
    e1out = numpy.empty(e1a.size, dtype='f8')
    e2out = numpy.empty(e1a.size, dtype='f8')

    for i in xrange(e1a.size):
        e1,e2 = _shearmult1(e1a[i], e2a[i], e1b[i], e2b[i])
        e1out[i] = e1
        e2out[i] = e2

    if isscalar:
        e1out=e1out[0]
        e2out=e2out[0]
    return e1out, e2out

def _shearmult1(e1a, e2a, e1b, e2b):
    dotp = e1a*e1b + e2a*e2b

    tmp = 1-e1b*e1b-e2b*e2b

    e1out=-9999.0
    e2out=-9999.0

    if tmp > 0.0:
        factor = (1.-sqrt(tmp)) / (e1b*e1b + e2b*e2b)

        tmp = 1+dotp

        if tmp != 0:
            e1out = (e1a + e1b + e2b*factor*(e2a*e1b - e1a*e2b))/(1+dotp)
            e2out = (e2a + e2b + e1b*factor*(e1a*e2b - e2a*e1b))/(1+dotp)

    return e1out, e2out

def _check4arr(arr1, arr2, arr3, arr4):
    arr1 = numpy.array(arr1, ndmin=1, copy=False) 
    arr2 = numpy.array(arr2, ndmin=1, copy=False) 
    arr3 = numpy.array(arr3, ndmin=1, copy=False) 
    arr4 = numpy.array(arr4, ndmin=1, copy=False) 

    if not (arr1.size == arr2.size == arr3.size == arr4.size):
        raise ValueError("all entries must be same size")
    return arr1, arr2, arr3, arr4

def compea4_test():
    e1o = 0.25
    e2o = 0.35

    e1p = 0.1
    e2p = 0.15

    # exp object, gaussian psf
    a4o=0.154
    a4p=0.0

    To = 10.0
    Tp = 5.0

    e1,e2,R,flags = correct(To,e1o,e2o,a4o,
                            Tp,e1p,e2p,a4p)
    print 'flags:',flags
    print 'e1: %0.6f e2: %0.6f R: %0.6f' % (e1,e2,R)
 
