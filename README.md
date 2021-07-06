# alnitak

**About:**

Takes a folder of .fits files, and, provided appropriate darks and flats exist, reduces and plate solves them.

**Dependencies:**

astropy

astroquery

photutils

numpy

scipy

barycorrpy

matplotlib

The code uses the plate solver api by nova.astrometry.net, so you need your own api key in order to use. Head to http://nova.astrometry.net/ and make an account (if you dont already have one. You can find your key at https://nova.astrometry.net/api_help. You can copy and paste your key in config.json next to API_Key

**How to use:**

Make sure you have your flatfield images in a folder called Flats, and your dark frames in a folder called Darks, both of which sit in the folder with your data images

**In command line:**

python /path/to/alnitak/data_reduction_run.py /path/to/folder/containing/images/to/reduce
