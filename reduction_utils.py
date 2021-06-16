from astropy.io import fits
from astropy.io.fits.hdu.compressed import CompImageHDU
import astroquery

from source_tools import Star_Tools
from header_modifications import Header_Info
from net.client import Client


import threading
import numpy as np
import glob
import os



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

        self.cli = Client()
        self.head = Header_Info(config_path, log=log)
        self.star = Star_Tools(log=log)
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
        self.finished = threading.Event()
        super(Calibration_Correct, self).__init__()

    def run(self):

        img_savepath, self.filename = os.path.split(self.file)
        self.img_savepath = os.path.join(img_savepath, 'reduced')

        self.science = fits.open(self.file, unit=True)
        self.log.info('Loaded in science file {}'.format(self.filename))
        self.log.info('Exposure time: {}, Filter: {}'.format(self.science[0].header['EXPTIME'], self.science[0].header['FILTER']))

        self.reduced = (self.science[0].data - self.dark_list[self.science[0].header['EXPTIME']])/(self.flat_list[self.science[0].header['FILTER']]/np.median(self.flat_list[self.science[0].header['FILTER']]))


        self.newfile_save = os.path.join(self.img_savepath, ('reduced_' + self.filename))
        self.header_out = self.science[0].header
        self.header = self.head.header_calculations(file=str(self.file), coords=self.coords, proper_motion=self.proper_motion)

        self.header['HISTORY'] = 'Dark corrected with mdark_{}.fits'.format(self.science[0].header['EXPTIME'])
        self.header['HISTORY'] = 'Flat corrected with mflat_{}.fits'.format(self.science[0].header['FILTER'])
        self.header['HISTORY'] = 'Reduced using data_reduction.py'

        compressed = fits.PrimaryHDU(data=self.reduced, header=self.header, uint=True)
        compressed.scale('uint16')
        compressed.writeto(self.newfile_save, overwrite=True) #Saves reduced image. Can be left alone or plate solved


        if self.platesolve == True:
            self.header_out = self.plate_solve(coords=self.coords, header_in=self.header)
            compressed = fits.PrimaryHDU(data=self.reduced, header=self.header_out, uint=True)
            compressed.scale('uint16')
            compressed.writeto(self.newfile_save, overwrite=True) #Overwrites saved image with plate solved version



    def plate_solve(self, coords, header_in):
        '''
        Description
        -----------

        Plate solves image that has been loaded into thread
        :returns
        self.header: Header object

        '''
        self.star_table = self.star.find_peaks(data=self.reduced)
        self.star_table_out = self.star.bad_pix(source_list=self.star_table, image=self.reduced)
        self.log.info(f'{self.filename}: Found {len(self.star_table_out)} stars')

        try_again = True
        submission_id = None

        while try_again:
            try:
                if not submission_id:
                    self.log.info(f'{self.filename}: Beginning solve')
                    self.wcs_header = self.ast.solve_from_source_list(self.star_table_out['starx'], self.star_table_out['stary'],
                                                                image_width=4096, image_height=4096,
                                                                solve_timeout=self.timeout, submission_id=submission_id)
                else:
                    self.wcs_header = self.ast.monitor_submission(submission_id, solve_timeout=self.timeout)

            except astroquery.exceptions.TimeoutError as e:
                submission_id = e.args[1]
            else:
                # got a result, so terminate
                try_again = False

        if self.wcs_header:
            self.log.info(f'{self.filename}: Solution found')
            for name in self.wcs_header.cards:
                if name[0] == 'SIMPLE' or name[0] == 'BITPIX':
                    continue
                elif name[0] == 'COMMENT' and name[1][0:5] == 'Index':
                    continue
                else:
                    self.header[name[0]] = name[1]
        else:
            self.log.critical('No solution could be found for {}, please check for image quality'.format(self.filename))
            # Code to execute when solve fails


        return self.header


class Simple_Reduce():
    def __init__(self, config_path, log):
        self.log = log
        self.head = Header_Info(config_path, log=log)

    def mflat_create(self, path, mdarks):
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

        def science_reduce(self, path, dark_list, flat_list, coords, proper_motion=(None, None)):
            img_savepath, filename = os.path.split(path)

            science = fits.open(path)
            self.log.info('---Loaded in science file {}---'.format(filename))
            self.log.info(science[0].header['EXPTIME'], science[0].header['FILTER'] )
            reduced = (science[0].data - dark_list[science[0].header['EXPTIME']])/(flat_list[science[0].header['FILTER']]/np.median(flat_list[science[0].header['FILTER']]))

            try:
                reduced = (science[0].data - dark_list[science[0].header['EXPTIME']])/(flat_list[science[0].header['FILTER']]/np.median(flat_list[science[0].header['FILTER']]))
            except:
                self.log.warning('Could not find corresponding dark, flat, or both for {}'.format(filename))
                #return None

            newfile_save = os.path.join(img_savepath, ('reduced_' + filename))
            header_out = science[0].header
            header = self.head.header_calculations(fits_in=path, coords=coords, proper_motion=prop_mot)

            header['HISTORY'] = 'Dark corrected with mdark_{}.fits'.format(science[0].header['EXPTIME'])
            header['HISTORY'] = 'Flat corrected with mflat_{}.fits'.format(science[0].header['FILTER'])
            header['HISTORY'] = 'Reduced using data_reduction.py'

            compressed = CompImageHDU(data=reduced, header=header)
            compressed.writeto(newfile_save)
