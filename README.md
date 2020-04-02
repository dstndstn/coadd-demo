# coadd-demo

This is a demo of coadding astronomical images to yield increased depth
or sensitivity.

Data are from one of the DES Supernova fields, near (41.42 -0.53).

Total number of exposures available: g=122, r=114, z=245.

Can produce jpegs with
python legacypipe/runbrick.py --survey-dir fake --radec 41.42 -0.53 --width 550 --height 550 --stage image_coadds --outdir f2 --force-all

For harder stretch, use sdss_rgb() with rgbscales multiplied by 5.  Results are in coadd-???.jpg.

Julia code demo (class notebook) is here
https://github.com/PerimeterInstitute/Computational-Physics-Course-Winter-2020/blob/master/class-2020-04-01/resampling.ipynb

Data (temporarily) here
https://portal.nersc.gov/project/cosmo/temp/dstn/coadd-data.tgz
