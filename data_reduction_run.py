from astropy.io import fits
from astropy.io.fits.hdu.compressed import CompImageHDU
from astropy.time import Time
from astroquery.astrometry_net import AstrometryNet
from astropy.table import QTable
from astropy import units as u
from astropy.stats import sigma_clipped_stats
from photutils import find_peaks, detect_threshold
from astropy.coordinates import SkyCoord, Angle
from source_tools import Star_Tools
from header_modifications import Header_Info
from reduction_utils import Calibration_Correct, Simple_Reduce
from net.client import Client
import threading
import logging

import numpy as np

import json
import glob
import os
import sys
import shutil

#Get coordinates for plate solving from http://simbad.u-strasbg.fr/simbad/
#Use icrs frame


def set_logger(name):
    #Creates a logger object that can be passed to threads
    formatter = logging.Formatter(fmt='[%(asctime)s] %(levelname)s: %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    try:
        os.remove('{}/log.txt'.format(path))
    except:
        pass
    handler = logging.FileHandler('{}/log.txt'.format(path), mode='w')
    handler.setFormatter(formatter)
    screen_handler = logging.StreamHandler(stream=sys.stdout)
    screen_handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    logger.addHandler(screen_handler)
    logger.info('Initialized logging')
    return logger


ast = AstrometryNet()
cli = Client()

ast.api_key = 'mxnnhvvdbxdwqrfa'
directory = os.path.dirname(sys.argv[0])
path = sys.argv[1]

logger = set_logger('logger')
calib = Simple_Reduce(directory, log=logger)

plate_solve = True
threading = True


for zz in sys.argv:
    if zz == '-ns':
        plate_solve = False
        logger.info('These science images will NOT be plate solved')
    if zz == '-nothreading':
        threading = False
        logger.info('These images will be reduced without threading. This will take longer')


max_threads = 5
thread_num = dict()
for x in range(0, max_threads): #Generates a dict the length of max_threads
    thread_num[f'Thread {x}'] = None
logger.info(f'Initialized {max_threads} threads')


#imports target coordinates if they exist as .json file
prop_mot = (None, None)
if os.path.exists(os.path.join(path + '/coords.json')):
    with open(os.path.join(path + '/coords.json')) as f:
        coords = json.load(f)
        if coords['pm_ra']:
            prop_mot = (coords['pm_ra'], coords['pm_dec'])

    g = coords['RA'] + ' ' + coords['DEC']
    radec_deg = SkyCoord(g, unit=(u.hourangle, u.deg), frame='icrs')


mdark_list = calib.mdark_create(path=path)

mflat_list = calib.mflat_create(path=path, mdarks=mdark_list)



'''--------SCIENCE--------'''
x = 0
img_savepath = path + '/reduced'
if os.path.exists(img_savepath):
    shutil.rmtree(img_savepath)
if not os.path.exists(img_savepath):#makes the reduced directory if not already exists
    os.mkdir(img_savepath)


#write if science reduce = None thing
for file in glob.glob(path + '/*.fit*'):
    science = fits.open(file, ignore_missing_end=True)
    try:
        ra_obj = science[0].header['RA_OBJ']
        dec_obj = science[0].header['DEC_OBJ']
        coords = SkyCoord(ra=ra_obj, dec=dec_obj, unit=u.degree, frame='icrs')
        radec_deg = coords
    except:
        pass
    science.close()
    break_loop = False
    if threading == True:
        while break_loop == False:
            for num in thread_num.keys():
                if thread_num[num] == None: # Starts new thread, should only run once per thread
                    new_thread = Calibration_Correct(config_path=directory, file=file, plate_solve=True,
                                                    flat_list=mflat_list, dark_list=mdark_list, coords=radec_deg, astrometrynet_instance=ast, log=logger)
                    new_thread.start()
                    thread_num[num] = new_thread
                    break_loop = True
                    break
                elif not thread_num[num].is_alive(): # Starts new thread once a previous one finishes
                    new_thread = Calibration_Correct(config_path=directory, file=file, plate_solve=True,
                                                    flat_list=mflat_list, dark_list=mdark_list, coords=radec_deg, astrometrynet_instance=ast, log=logger)
                    new_thread.start()
                    thread_num[num] = new_thread
                    break_loop = True
                    break
                else:
                    continue


    else:
        path, filename = os.path.split(file)

        calib.science_reduce(path=file, flat_list=mflat_list, dark_list=mdark_list,
                                        directory=img_savepath, coords=radec_deg, save=True, platesolve=True)



'''
    #plate solves the reduced image if specified
    wcs_header = None
    n = 0
    try_again = True
    if plate_solve:
        print('Plate solving {}'.format(filename))
        print("Determining background stats", flush=True)
        mean, median, std = sigma_clipped_stats(reduced, sigma=3.0)
        threshold = detect_threshold(reduced, nsigma=5)
        print("Finding sources", flush=True)
        sources = find_peaks(reduced, threshold=threshold, box_size=25, border_width=50)
        sources.sort('peak_value')

        stars, peaks, star_table = bad_pix(sources, reduced)
        star_table_out = fwhm(star_table, reduced, name=x)
        x+=1
        print('{} sources found, {} were not bad pixels.\n{} pass fwhm criteria'.format(len(sources), len(stars), len(star_table_out)))

        if sources and median <=8000:
            if len(star_table_out) < 10:
                print('Fewer than 10 sources found, performing standard reduction. Check to see if image is good')
                continue
            else:
                #print('{} sources found:\n{}'.format(len(sources), sources))
                sources.sort('peak_value')
                sources.reverse()
                #print(type(sources))
                submission_id = None
            while n <= 2 and try_again:
                try:
                    if not submission_id:
                        wcs_header = ast.solve_from_source_list(star_table_out['starx'], star_table_out['stary'], submission_id=submission_id,
                                                                center_ra=radec_deg.ra.degree, center_dec=radec_deg.dec.degree,
                                                                image_width=4096, image_height=4096,
                                                                solve_timeout=180, scale_units='arcsecperpix', scale_lower=0.36)
                        print('---This is the wcs header: {}---'.format(wcs_header))
                        if len(wcs_header) == 0:
                            n += 1
                            print('Solve for {} failed, trying again up to 2 times'.format(filename))
                    else:
                        wcs_header = ast.monitor_submission(submission_id=submission_id, timeout=500)
                except:
                    print('Solve timed out, retrying up to 2 times')
                    n += 1

            if wcs_header:
                print('---{}---'.format(wcs_header))
                # Code to execute when solve succeeds
            else:
                # Code to execute when solve fails
                print('---SOLVE TIMEOUT---')
                    #hdu = fits.PrimaryHDU(reduced, header = wcs_header)
                    #hdu.writeto(newfile_save)

        else:
            print('Median counts greater than 10,000, assuming cloudy')
    else:
        hdu = fits.PrimaryHDU(reduced, header = header_out)
        hdu.writeto(newfile_save)
    print('{} reduction complete! Moving on to the next'.format(filename))

    '''


print('Data reduction complete')








'''

FOR MAKING MDARK/MFLATS

for file in glob.glob('/Users/owenalfaro/Desktop/Astro_Images/M42/*.fit'):
    dlist.append(fits.getdata(file))

medlist = np.array(dlist)
outfile = np.median(medlist, axis = 0)
done = fits.PrimaryHDU(outfile)
done.writeto('/Users/owenalfaro/Desktop/test.fits')
'''
'''

make master flat dark
make master science dark
subtract flat dark from flat
make mflat
    MEDIAN
    NORMALIZE
reduced = (science - science dark)/(mflat/median()
'''
