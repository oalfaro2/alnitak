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
import threading
import logging

import numpy as np

import json
import glob
import os
import sys
import shutil
import time
import datetime

#Get coordinates for plate solving from http://simbad.u-strasbg.fr/simbad/
#Use icrs frame
now = datetime.datetime.now()
directory = os.path.dirname(sys.argv[0])
with open(os.path.join(directory, 'config.json')) as f:
    conf = json.load(f)

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
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    logger.addHandler(screen_handler)
    logger.info('Initialized logging')
    return logger


ast = AstrometryNet()

ast.api_key = conf['API_Key'] #insert api key from nova.astrometry.net
max_threads = int(conf['Max_Threads'])
path = sys.argv[1]

logger = set_logger('logger')
calib = Simple_Reduce(directory, log=logger, astrometry_net=ast)

plate_solve = True
threading = True


for zz in sys.argv:
    if zz == '-ns' or zz == '-nosolve':
        plate_solve = False
        logger.info('These science images will NOT be plate solved')
    if zz == '-nt' or zz == '-nothreading':
        threading = False
        logger.info('These images will be reduced without threading. This will take longer')


thread_num = dict()
for x in range(0, max_threads): #Generates a dict the length of max_threads
    thread_num[f'Thread {x}'] = None
logger.info(f'Initialized {max_threads} threads')

coords_in_json = False
#imports target coordinates if they exist as .json file
prop_mot = (None, None)
if os.path.exists(os.path.join(path + '/coords.json')):
    with open(os.path.join(path + '/coords.json')) as f:
        coords = json.load(f)
        coords_in_json = True
        if coords['pm_ra']:
            prop_mot = (coords['pm_ra'], coords['pm_dec'])

    g = coords['RA'] + ' ' + coords['DEC']
    radec_deg = SkyCoord(g, unit=(u.hourangle, u.deg), frame='icrs')


mdark_list = calib.mdark_create(path=path)

mflat_list = calib.mflat_create(path=path, mdarks=mdark_list)



img_savepath = path + '/reduced'
bad_img_path = path + '/bad'
if os.path.exists(img_savepath):
    shutil.rmtree(img_savepath)
if not os.path.exists(img_savepath):#makes the reduced directory if not already exists
    os.mkdir(img_savepath)

if os.path.exists(bad_img_path):
    shutil.rmtree(bad_img_path)
if not os.path.exists(bad_img_path):#makes the bad image directory if not already exists
    os.mkdir(bad_img_path)


#write if science reduce = None thing
for file in glob.glob(path + '/*.fit*'):
    science = fits.open(file, ignore_missing_end=True)
    if not coords_in_json:
        try:
            ra_obj = science[0].header['RA_OBJ']
            dec_obj = science[0].header['DEC_OBJ']
            radec_deg = SkyCoord(ra=ra_obj, dec=dec_obj, unit=u.degree, frame='icrs')
        except:
            radec_deg = None
            pass
    science.close()
    break_loop = False
    if threading == True:
        while break_loop == False:
            for num in thread_num.keys():
                if thread_num[num] == None: # Starts new thread, should only run once per thread
                    new_thread = Calibration_Correct(config_path=directory, file=file, plate_solve=plate_solve,
                                                    flat_list=mflat_list, dark_list=mdark_list, coords=radec_deg, astrometrynet_instance=ast, log=logger)
                    new_thread.start()
                    thread_num[num] = new_thread
                    break_loop = True
                    break
                elif not thread_num[num].is_alive(): # Starts new thread once a previous one finishes
                    new_thread = Calibration_Correct(config_path=directory, file=file, plate_solve=plate_solve,
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
                             coords=radec_deg, platesolve=plate_solve)

for n in thread_num.keys():
    try:
        while thread_num[n].is_alive():
            time.sleep(1)
    except:
        time.sleep(1)
        continue
#       Waits until every thread has finished to print complete
end = datetime.datetime.now()
total = end-now
logger.info(f'Data reduction complete. Took {total.total_seconds()}s')
