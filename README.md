# alnitak

**About:**

Takes a folder of .fits files, and, provided appropriate darks and flats exist, reduces and plate solves them.

**Dependencies:**

astropy

astroquery

photutils

glob

numpy

scipy


**How to use:**

Make sure you have your flatfield images in a folder called Flats, and your dark frames in a folder called Darks, both of which sit in the folder with your data images

**In command line:**

python /path/to/alnitak/data_reduction_run.py /path/to/folder/containing/images/to/reduce
