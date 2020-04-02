# coadd-demo

This is a demo of coadding astronomical images to yield increased depth
or sensitivity.

Data are from one of the DES Supernova fields, near (41.42 -0.53).

Total number of exposures available: g=122, r=114, z=245.

Can produce jpegs with
python legacypipe/runbrick.py --survey-dir fake --radec 41.42 -0.53 --width 550 --height 550 --stage image_coadds --outdir f2 --force-all

Julia code demo (class notebook) is here
https://github.com/PerimeterInstitute/Computational-Physics-Course-Winter-2020/blob/master/class-2020-04-01/resampling.ipynb

