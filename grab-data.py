import pylab as plt
import numpy as np
from glob import glob
import os
from astrometry.util.file import file_size
from astrometry.util.fits import *
from legacypipe.survey import *
import fitsio
from astrometry.util.util import Tan

# NERSC
C = fits_table('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz')

K = np.flatnonzero((np.hypot(C.ra - 41.3, C.dec - -0.6) < 0.1) * np.array([c.strip() == 'S17' for c in C.ccdname]))

IG = K[C.filter[K] == 'g']
IR = K[C.filter[K] == 'r']
IZ = K[C.filter[K] == 'z']

G = C[IG]
G.cut(G.exptime == 175.)
R = C[IR]
R.cut(R.exptime == 150.)
Z = C[IZ]
Z.cut(Z.exptime == 200.)

G.cut((G.ccdzpt > 24.9) * (G.ccdzpt < 25.3))
R.cut((R.ccdzpt > 25.2) * (R.ccdzpt < 25.4))
Z.cut((Z.ccdzpt > 24.8) * (Z.ccdzpt < 25.1))

for nb in [1, 2, 4, 16, 100]:
    best = merge_tables([
        X[np.argsort(X.sig1 * X.fwhm)][:nb]
        for X in [G,R,Z]
    ])
    best.ccd_cuts[:] = 0
    fn = 'survey-ccds-best-%i.fits' % nb
    best.writeto(fn)
    print('Wrote', fn)

survey = LegacySurveyData('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9')

# ~550-by-550-pixel box:
r1,r2 = [41.40, 41.44]
d1,d2 = [-0.55, -0.51]

radecpoly = np.array([[r1,d1],[r2,d1],[r2,d2],[r1,d2]])

for band,Z in [('g',G), ('r',R),('z',Z)]:
    for i,t in enumerate(Z[np.argsort(Z.sig1 * Z.fwhm)][:100]):
        im = survey.get_image_object(t)
        print(im)
        tim = im.get_tractor_image(radecpoly=radecpoly, old_calibs_ok=True)
        print('Read', tim, tim.shape)
        H,W = tim.shape
        refx = W//2 + 1.
        refy = H//2 + 1.
        refr,refd = tim.subwcs.pixelxy2radec(refx, refy)
        r1,d1     = tim.subwcs.pixelxy2radec(W,    refy)
        r2,d2     = tim.subwcs.pixelxy2radec(refx, H)
        cosdec = np.cos(np.deg2rad(refd))
        cd11 = (r1 - refr) * cosdec / (W - refx)
        cd12 = (r2 - refr) * cosdec / (H - refy)
        cd21 = (d1 - refd) / (W - refx)
        cd22 = (d2 - refd) / (H - refy)
        tanwcs = Tan(refr, refd, refx, refy, cd11, cd12, cd21, cd22, float(W), float(H))

        hdr = fitsio.FITSHDR()
        tanwcs.add_to_header(hdr)
        hdr.add_record(dict(name='BAND', value=band))
        hdr.add_record(dict(name='CAMERA', value=im.camera))
        hdr.add_record(dict(name='EXPNUM', value=im.expnum))
        hdr.add_record(dict(name='CCDNAME', value=im.ccdname))
        hdr.add_record(dict(name='X0', value=tim.x0))
        hdr.add_record(dict(name='Y0', value=tim.y0))
    
        fn = '/global/cscratch1/sd/dstn/coadd/image-%s-%03i.fits' % (band, i)
        fitsio.write(fn, tim.getImage(), header=hdr, extname='IMAGE', clobber=True)
        fitsio.write(fn, tim.getInvError(), extname='INVERR')
        print('Wrote', fn)

