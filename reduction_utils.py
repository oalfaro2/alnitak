from astropy.io import fits
from astropy.io.fits.hdu.compressed import CompImageHDU
import astroquery

from source_tools import Star_Tools
from header_modifications import Header_Info


import threading
import numpy as np
import glob
import os
import json



class Calibration_Correct(threading.Thread):
    def __init__(self, config_path, file, plate_solve, flat_list, dark_list, coords,
                astrometrynet_instance, log, proper_motion=None):
        '''
        Description
        -----------
        A threading object which calibrates science fits images and plate solves them via astroquery.astrometry_net

        Parameters
        ----------
        config_path: STR
            The path to the config.json with information about observatory
        file: STR
            Filepath to science image to be reduced
        plate_solve: BOOL
            If True, will plate solve image
        dark_list: DICT
            A dict with the master darks. Ex: {'EXPOSURE TIME': mdark}
            Master darks are numpy arrays
        coords: SkyCoord obj
            Input coordinates of science image
        astrometrynet_instance:
            Instance of the astroquery package after login
        proper_motion: TUPLE
            Tuple of the proper motion in milli arcseconds/year

        '''

        self.head = Header_Info(config_path, log=log)
        self.star = Star_Tools(log=log, config_path=config_path)
        self.reduce = Simple_Reduce(config_path, log=log, astrometry_net=astrometrynet_instance)
        self.ast = astrometrynet_instance
        self.log = log

        self.error_message = '{"error": "no calibration data available for job'
        self.file = file
        self.timeout = 360
        self.platesolve = plate_solve
        self.flat_list = flat_list
        self.dark_list = dark_list
        self.coords = coords
        self.proper_motion = proper_motion
        super(Calibration_Correct, self).__init__()

    def run(self):

        self.reduce.science_reduce(path=self.file, dark_list=self.dark_list, flat_list=self.flat_list, coords=self.coords, proper_motion=self.proper_motion, platesolve=self.platesolve)





class Simple_Reduce():
    def __init__(self, config_path, log, astrometry_net):
        self.log = log
        self.head = Header_Info(config_path, log=log)
        self.ast = astrometry_net
        self.star = Star_Tools(log=log, config_path=config_path)
        self.error_message = '{"error": "no calibration data available for job'
        with open(os.path.join(config_path, 'config.json')) as f:
            self.config = json.load(f)


    def mflat_create(self, path, mdarks):
        '''
        Description
        -----------
        Creates a list of master flats. Dark subtracts them if there is a dark with the same exposure time as the flat

        :param path: STR
            path to folder containing flats
        :param mdarks: DICT
            dictionary of mdarks with keys as exposure times
        :return:
            mflats: DICT
            dictionary of master flats with keys as filters
        '''
        flat_list = dict()
        flat_header = dict()
        mflats = dict()
        temp_list = None

        for file in glob.glob(os.path.join(path, '*Flat*/') + '*.fit*'):
            fsave, filename = os.path.split(file)   #gets pathname without *s
            if str(filename)[0:len('mflat')] == 'mflat':
                continue
            flat_dark = False
            flat_data = fits.open(file, unit=False)
            self.log.info('Loaded in flat frame {}'.format(filename))
            if flat_data[0].header['FILTER'] not in flat_list.keys(): #If the exposure time is not in the dict, creates new key
                keyname = flat_data[0].header['FILTER']
                flat_list[keyname] = list()
                flat_header[keyname] = list()
                flat_header[flat_data[0].header['FILTER']] = flat_data[0].header

                self.log.info('New flat filter and exposure detected: {} {}s'.format(keyname, flat_data[0].header['EXPTIME']))
                if flat_data[0].header['SET-TEMP'] <= 15 and flat_data[0].header['SET-TEMP'] >= 0:
                    self.log.warning('Warning, flat is of poor quality and will be extremely noisy. \nThis might make plate solving more difficult \nCCD Set temp: {}'.format(flat_data[0].header['SET-TEMP']))

            for time in mdarks.keys():
                if flat_data[0].header['EXPTIME'] == time:
                    wrk_data = flat_data[0].data - mdarks[time]  #dark subtracts each flat
                    temp_list = flat_list[flat_data[0].header['FILTER']]
                    temp_list.append(wrk_data)
                    flat_list[flat_data[0].header['FILTER']] = temp_list
                    flat_dark = True
                    break

            if flat_dark:
                flat_data.close()
                continue
            else:
                temp_list = flat_list[flat_data[0].header['FILTER']]# If there are no flat darks it just adds data to dict
                temp_list.append(np.array(flat_data[0].data))
                flat_list[flat_data[0].header['FILTER']] = temp_list

            flat_data.close()

        for filter in flat_list.keys():
            mflat_ = np.median(flat_list[filter], axis = 0)
            self.log.info(f'Master Flat {filter} calculated')

            mflat = fits.PrimaryHDU(data=mflat_, header = flat_header[filter], scale_back=True)
            mflat.scale('int16')
            mflats[filter] = mflat.data # Adds to list of Mflats
            if not os.path.exists(os.path.join(fsave, ('mflat_{}.fits'.format(filter)))):
                mflat.writeto(os.path.join(fsave, ('mflat_{}.fits'.format(filter))))
                self.log.info('Mflat {} saved to {}'.format(filter, os.path.join(fsave, ('mflat_{}.fits'.format(filter)))))
            else:
                self.log.info('Mflat {} already exists, continuing'.format(filter))

        return mflats

    def mdark_create(self, path):
        '''
        Description
        -----------
        Creates a list of master darks in a dictionary with keys as exposure times

        :param path: STR
            Path to folder containing dark .fits images
        :return:
            mdarks: DICT
                List of master darks
        '''
        exp_times = dict()
        header_times = dict()
        mdarks = dict()
        temp_list = None
        #opens and makes master darks
        for file in glob.glob(os.path.join(path, '*Dark*/') + '*.fit*'):

            save, filename = os.path.split(file)

            if str(filename)[0:len('mdark')] == 'mdark':
                continue

            data = fits.open(file, unit=True)
            fitsdata = fits.getdata(file)
            self.log.info('Loaded in dark frame {}'.format(filename))

            if data[0].header['EXPTIME'] not in exp_times.keys(): #If the exposure time is not in the dict, creates new key
                keyname = data[0].header['EXPTIME']
                exp_times[keyname] = list()
                header_times[keyname] = list()
                mdarks[keyname] = list()
                self.log.info('New dark exposure detected: {}s'.format(keyname))
                if data[0].header['SET-TEMP'] <= 15 and data[0].header['SET-TEMP'] >= 0:
                    self.log.warning('Warning, dark is of poor quality and will be extremely noisy. \nThis might make plate solving more difficult \nCCD Set temp: {}'.format(data[0].header['SET-TEMP']))

            for key in exp_times.keys():
                if data[0].header['EXPTIME'] == key:#groups darks in lists by exposure time
                    temp_list = exp_times[key]
                    temp_list.append(fitsdata)
                    exp_times[key] = temp_list
                    if len(header_times[key]) == 0:
                        header_times[key] = data[0].header
                    break
            data.close()
        # saves the mdarks
        for xx in exp_times.keys():
            dark_list = exp_times[xx]
            med_dark = np.median(dark_list, axis = 0)

            self.log.info('Master dark {}s calculated'.format(xx))
            dark = fits.PrimaryHDU(data=med_dark, header=header_times[xx], scale_back=True)#.scale('int16', bzero=32768)
            dark.scale('int16')
            if len(mdarks[xx]) == 0:
                mdarks[xx] = dark.data #Creates list of mdarks with keys of exposure time

            if not os.path.exists(os.path.join(save, ('mdark_{}.fits'.format(xx)))):
                dark.writeto(os.path.join(save, ('mdark_{}.fits'.format(xx))))
                self.log.info('Mdark saved to {}'.format(os.path.join(save, ('mdark_{}.fits'.format(xx)))))
            else:
                self.log.info('Mdark {}s already exists, continuing'.format(xx))

        return mdarks # A list with the median darks based on exposure


    def science_reduce(self, path, dark_list, flat_list, coords, platesolve=False, timeout=None, proper_motion=(None, None)):
        '''
        Description
        -----------
        Dark subtracts and flat divides inputted science .fits image with appropriate flat and dark
        Passes image to plate solver for uploading to nova.astrometry.net
        Saves reduced, plate solved image

        :param path: STR
            path to science .fits image
        :param dark_list: DICT
            Dictionary of master darks. Keys are exposure times
        :param flat_list: DICT
            Dictionary of master flats. Keys are Filters
        :param coords: Astropy Coordinates Obj
            Initial coordinates of science .fits image
        :param platesolve: BOOL
            If true, plate solves image using nova.astrometry.net
        :param timeout: INT
            Integer time in seconds for astrometry.net timeout
        :param proper_motion: TUPLE
            Tuple of measured proper motion of host star, can be found at
            https://exofop.ipac.caltech.edu/tess/
            default (None, None). Not required
        :return:
        '''

        if timeout:
            pass
        else:
            timeout = 180
        img_savepath, filename = os.path.split(path)
        self.path = path
        img_savepath_ = os.path.join(img_savepath, 'reduced')
        bad_folder = os.path.join(img_savepath, 'bad')
        newfile_save = os.path.join(img_savepath_, ('reduced_' + filename))
        bad_newfile_save = os.path.join(bad_folder, ('reduced_' + filename))

        science = fits.open(path, unit=True)
        self.log.info('Loaded in science file {}'.format(filename))
        self.log.info(
            'Exposure time: {}, Filter: {}'.format(science[0].header['EXPTIME'], science[0].header['FILTER']))

        reduced = (science[0].data - dark_list[science[0].header['EXPTIME']]) / (
                    flat_list[science[0].header['FILTER']] / np.median(
                flat_list[science[0].header['FILTER']]))

        header_in = science[0].header

        if platesolve:
#       Plate solving here
            header_out, success = self.plate_solve(coords=coords, header_in=header_in,
                                                      data_in=reduced, filename=filename, timeout=timeout)
            if success:
#           If solve is successful, saves to newfile save
                header_out['HISTORY'] = 'Dark corrected with mdark_{}.fits'.format(science[0].header['EXPTIME'])
                header_out['HISTORY'] = 'Flat corrected with mflat_{}.fits'.format(science[0].header['FILTER'])
                header_out['HISTORY'] = 'Reduced using https://github.com/oalfaro2/alnitak'

                compressed = fits.PrimaryHDU(data=reduced, header=header_out, uint=True)
                compressed.scale('uint16')
                compressed.writeto(newfile_save, overwrite=True)  # Overwrites saved image with plate solved version
            else:
#           If not successful, saves to different directory. Performs what header calculations it can
                header_out = self.head.header_calculations(plate_solved=False, header=header_in, file=self.path)
                header_out['HISTORY'] = 'Could not be plate solved'

                compressed = fits.PrimaryHDU(data=reduced, header=header_out, uint=True)
                compressed.scale('uint16')
                compressed.writeto(bad_newfile_save, overwrite=True)
        else:
            header = self.head.header_calculations(file=str(path), coords=coords)

            header['HISTORY'] = 'Dark corrected with mdark_{}.fits'.format(science[0].header['EXPTIME'])
            header['HISTORY'] = 'Flat corrected with mflat_{}.fits'.format(science[0].header['FILTER'])
            header['HISTORY'] = 'Reduced using https://github.com/oalfaro2/alnitak'

            compressed = fits.PrimaryHDU(data=reduced, header=header, uint=True)
            compressed.scale('uint16')
            compressed.writeto(newfile_save,
                               overwrite=True)  # Saves reduced image. Can be left alone or plate solved


    def plate_solve(self, coords, header_in, data_in, timeout, filename):
        '''
        Description
        -----------

        Plate solves image that has been loaded into thread
        :returns
        self.header: Header object

        '''
        star_table = self.star.find_peaks(data=data_in)
        star_table_out, med_fwhm = self.star.bad_pix(source_list=star_table, image=data_in)
        self.log.info(f'{filename}: Found {len(star_table_out)} stars. Median FWHM: {med_fwhm} px')
        if len(star_table_out) <= int(self.config['Min_Stars']):
            self.log.warning(f'{filename}: Less than {self.config["Min_Stars"]} stars found, aborting solve')
            success = False
            header_out = None
            return header_out, success

        try_again = True
        submission_id = None
        success = True

        while try_again:
            try:
                if not submission_id:
                    self.log.info(f'{filename}: Beginning solve')
                    wcs_header = self.ast.solve_from_source_list(star_table_out['starx'],
                                                                 star_table_out['stary'],
                                                                 image_width=4096, image_height=4096,
                                                                 solve_timeout=timeout,
                                                                 submission_id=submission_id)
                else:
                    wcs_header = self.ast.monitor_submission(submission_id, solve_timeout=timeout)

            except astroquery.exceptions.TimeoutError as e:
                submission_id = e.args[1]
            else:
                # got a result, so terminate
                try_again = False

        if wcs_header:
            self.log.info(f'{filename}: Solution found')
            for name in wcs_header.cards:
                if name[0] == 'SIMPLE' or name[0] == 'BITPIX':
                    continue
                elif name[0] == 'COMMENT' and name[1][0:5] == 'Index':
                    continue
                else:
                    header_in[name[0]] = name[1]
        else:
            self.log.critical('No solution could be found for {}, please check for image quality'.format(filename))
            success = False
            header_out = None
            return header_out, success
            # Code to execute when solve fails

        header_out = self.head.header_calculations(plate_solved=True, header=header_in, file=self.path)

        return header_out, success
