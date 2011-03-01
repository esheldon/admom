try:
    import fimage
except:
    print("Could not import fimage")

from numpy import sqrt
def correct(image, psf):
    mom_image = fimage.stat.fmom(image)
    mom_psf = fimage.stat.fmom(psf)

    # simply subtract the covariance matrix
    
    cov = mom_image['cov'] - mom_psf['cov']

    T = cov[0] + cov[2]
    e1 = (cov[2]-cov[0])/T
    e2 = 2*cov[1]/T
    e=sqrt(e1**2 + e2**2)

    out = {'cen':mom_image['cen'],
           'cov':mom_image['cov'],
           'cen_psf':mom_psf['cen'],
           'cov_psf':mom_psf['cov'],
           'cov_corr':cov,
           'e1':e1,
           'e2':e2,
           'e':e}

    return out


