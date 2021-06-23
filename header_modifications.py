from astropy.time import Time
from astropy.io import fits
from barycorrpy import utc_tdb
from astropy import units as u
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, Angle, FK5
import numpy as np
import datetime
import os
import json



class Header_Info():
    def __init__(self, config_path, log):
        self.log = log
        with open(os.path.join(config_path, 'config.json')) as f:
            self.config = json.load(f)
            for x in self.config.keys():
                try:
                    self.config[x] = float(self.config[x])
                except:
                    continue


    def years_since_2015(self, year):
        return year - 2015

    def get_lst(self, date, long):

        j2000 = datetime.datetime(2000, 1, 1, 12, 0, 0, tzinfo=datetime.timezone.utc)
        date = datetime.datetime(date.year, date.month, date.day, date.hour, date.minute, date.second, tzinfo=datetime.timezone.utc)
        days = (date-j2000).total_seconds()/(3600*24)
        decimal_hours_of_day = (int(date.hour) + int(date.minute)/60 + int(date.second)/3600)/24

        lst = 100.46 + 0.985647*days + long + 15*decimal_hours_of_day
        while lst > 360:
            lst -= 360
        while lst < 0:
            lst += 360
        return lst


    def radec_to_altaz(self, ra, dec, latitude, longitude, lst):
        ha = (lst - ra) % 24
        if ha > 12:
            ha -= 24
        # Convert to degrees
        #test
        ha *= 15
        (dec_r, latitude_r, longitude_r, HA_r) = np.radians([dec, latitude, longitude, ha])

        alt_r = np.arcsin(np.sin(dec_r)*np.sin(latitude_r)+np.cos(dec_r)*np.cos(latitude_r)*np.cos(HA_r))
        az_r = np.arctan2(-np.cos(dec_r)*np.sin(HA_r), np.sin(dec_r)*np.cos(latitude_r)-np.sin(latitude_r)*np.cos(dec_r)*np.cos(HA_r))
        (az, alt) = np.degrees([az_r, alt_r])
        az = az % 360
        self.log.info('ALTAZ: {} {}'.format(alt, az))

    def coord_transform(self, ra, dec, year):
        coords_j2000 = SkyCoord(ra=ra, dec=dec, unit=u.deg, frame='icrs')
        coords_apparent = coords_j2000.transform_to(FK5(equinox='J{}'.format(year)))
        return coords_apparent

    def proper_motion_adjust(self, coordinates, proper_motion, years_since_2015):
        ra_mas_to_deg = float(proper_motion[0])*0.00099999995874704*0.000277778
        dec_mas_to_deg = float(proper_motion[0])*0.00099999995874704*0.000277778
        ra_out = coordinates.ra + ra_mas_to_deg*u.deg
        dec_out = coordinates.dec + dec_mas_to_deg*u.deg
        radec_deg_adjusted = SkyCoord(ra=ra_out, dec=dec_out, unit=u.deg, frame='icrs')
        radec_deg_adjusted_ = SkyCoord(ra=ra_out, dec=dec_out, unit=u.deg, frame=FK5(equinox='J2021'))
        self.log.info(radec_deg_adjusted.ra.deg)
        self.log.info(radec_deg_adjusted_.ra.deg)
        return radec_deg_adjusted

    def header_calculations(self, file=None, coords=None, proper_motion=None, header=None, plate_solved=False):
        dir_, filename = os.path.split(file)
        if plate_solved:
            ra = header['CRVAL1']
            dec = header['CRVAL2']
            coords_solved = SkyCoord(ra=ra, dec=dec, unit=u.deg, frame='icrs')
            header_out = self.calculations(header=header, coords=coords_solved)
        else:
            fits_in = fits.open(file)
            header_out = self.calculations(header=fits_in[0].header, coords=coords)
        self.log.info(f'{filename}: Header calculations complete')
        return header_out

    def calculations(self, coords=None, proper_motion=None, header=None):
        #requires fits file with header
        #See if alt and az are accurate enough

        jd = header['JD']    #Julian date at start of exposure
        date = header['DATE-OBS']
        date_obj = datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S')
        #print(date, date_obj)
        jd_helio = header['JD-HELIO']    #Heliocentric JD at start of exposure
        mid_exp = header['EXPOSURE']/2   #middle of exposure
        jd_utc_mid = Time(date_obj + datetime.timedelta(seconds=mid_exp), scale='utc').jd   #JD at mid exposure
        jd_utc_mid_frame = Time(date_obj + datetime.timedelta(seconds=mid_exp), scale='utc')   #JD at mid exposure but the Time frame

        bjd_tdb = utc_tdb.JDUTC_to_BJDTDB(jd_utc_mid, lat=self.config['Site_lat'], longi=self.config['Site_long'], alt=self.config['Site_alt'])

#       Converts julian date to barycentric julian date
        header['SITELAT'] = (self.config['Site_lat'], 'Geographic latitude of observatory')
        header['SITELONG'] = (self.config['Site_long'], 'Geographic longitude of observatory')
        header['SITEALT'] = (self.config['Site_alt'], 'Altitude of observatory (m) above MSL')
        header['JD_SOBS'] = (jd, 'JD (UTC) at the start of observation')
        header['JD_UTC'] = (jd_utc_mid, 'JD (UTC) at the midpoint of the exposure')
        header['BJD_TDB'] = (float(bjd_tdb[0]), 'Barycentric JD at the midpoint of the exposure')

        #print(f'{filename} JD SObs: {jd}')
        if coords:
            coords_2k = self.coord_transform(ra=coords.ra, dec=coords.dec, year=date_obj.year)
            earth = EarthLocation(lon=self.config['Site_long']*u.deg, lat=self.config['Site_lat']*u.deg, height=self.config['Site_alt']*u.m)
            altaz = coords_2k.transform_to(frame=AltAz(location=earth, obstime=date_obj))
            #print(f'{filename} Alt Az: {altaz.alt}, {altaz.az}')
            zenith_dist = 90 - np.absolute(float(altaz.alt.deg))
            #print(f'{filename} Zenith distance: {zenith_dist}'.format(zenith_dist))
            ha_mean = jd_utc_mid_frame.sidereal_time('mean', longitude=self.config['Site_long']*u.deg) - coords.ra
            ha_apparent = jd_utc_mid_frame.sidereal_time('apparent', longitude=self.config['Site_long']*u.deg) - coords.ra
            airmass = 1/(np.cos(np.pi/2 - altaz.alt.deg*np.pi/180))

            #lst = get_lst(date=date_obj, long=self.config['Site_long'])
            #radec_to_altaz(ra=coords.ra.deg, dec=coords.dec.deg, latitude=self.config['Site_lat'], longitude=self.config['Site_long'], lst=lst)




            header['ALT_OBJ'] = (altaz.alt.deg, 'Target Altitude at mid exposure')
            header['AZ_OBJ'] = (altaz.az.deg, 'Target Azimuth at mid exposure')
            header['HA_AP'] = (float(ha_mean.to_string(unit=u.hour, decimal=True)), 'Mean Hourangle at mid exposure')
            header['HA_MEAN'] = (float(ha_apparent.to_string(unit=u.hour, decimal=True)), 'Apparent Hourangle at mid exposure')
            header['ZD_OBJ'] = (zenith_dist, 'Distance from zenith at mid exposure')
            header['AIRMASS'] = (airmass, 'Target airmass at mid exposure')
            header['RA_OBJ'] = (coords.ra.hour, 'Target Right ascension (hours)')
            header['DEC_OBJ'] = (coords.dec.deg, 'Target declination (degrees)')
            header['RAOBJ2K'] = (coords_2k.ra.hour, 'Target J2000 Right ascension (hours)')
            header['DECOBJ2K'] = (coords_2k.dec.deg, 'Target J2000 declination (degrees)')

        return header
