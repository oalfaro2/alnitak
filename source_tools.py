from astropy.table import Table
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import find_peaks
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

saturation = 30000
def hms_to_degree(angle):
    string = angle.split(':')
    deg_angle = 15*float(string[0]) + 15*float(string[1])/60 + 15*float(string[2])/3600
    return deg_angle


class Star_Tools():
    def __init__(self, log):
        self.log = log
        pass


    def find_peaks(self, data):
        star_table = find_peaks(data=data, threshold=1.1*np.median(data), box_size=50, border_width=25)
        star_table.sort('peak_value', reverse=True)
        return star_table


    def bad_pix(self, source_list, image):
        mean, median, stdev = sigma_clipped_stats(image, sigma=3)
        n = 0
        bad_pix_num = 0
        stars = []
        peaks = []
        good_x = []
        good_y = []
        for _ in source_list:
            bad_pixel = False
            x_cent = source_list['x_peak'][n]
            y_cent = source_list['y_peak'][n]
            peak = image[y_cent, x_cent]
            n += 1
            if peak >= (saturation * 2) ** 2:
                bad_pixel = True
            pixels = [(y_cent, x_cent + 1), (y_cent, x_cent - 1), (y_cent + 1, x_cent), (y_cent - 1, x_cent)]
            for value in pixels:
                if image[value[0], value[1]] < 1.1 * median:
                    bad_pixel = True
            if bad_pixel:
                bad_pix_num += 1
                continue
            star = (x_cent, y_cent)
            good_x.append(x_cent)
            good_y.append(y_cent)
            stars.append(star)
            peaks.append(peak)
        star_table = Table([good_x, good_y], names = ('starx', 'stary'))

        self.log.debug(f'{bad_pix_num} sources removed from list of {len(source_list)} potential stars.')
        return star_table


    def fwhm(self, star_list, data, name, splitdata=False, radius=30):
        a = 0
        fwhm_list = []
        good_stars_x = []
        good_stars_y = []
        # a = 0
        for star in star_list:
            x_cent = star[0]
            y_cent = star[1]
            star = data[int(y_cent-radius):int(y_cent+radius), int(x_cent-radius):int(x_cent+radius)]
            starx, stary = np.indices(star.shape)
            r = np.sqrt((stary - radius)**2 + (starx - radius)**2)
            r = r.astype(np.int)

            tbin = np.bincount(r.ravel(), star.ravel())
            nr = np.bincount(r.ravel())
            radialprofile = tbin / nr

            if len(radialprofile) != 0:
                maximum = (max(radialprofile)/np.median(radialprofile)) - 1
                if maximum == 0:
                    continue
                else:
                    radialprofile = (radialprofile / np.median(radialprofile)) - 1#maximum
                f = np.linspace(0, len(radialprofile)-1, len(radialprofile))
                mean = np.mean(radialprofile)
                sigma = np.std(radialprofile)
                '''
                try:
                    popt, pcov = curve_fit(gaussianfit, f, radialprofile, p0=[1 / (np.sqrt(2 * np.pi)), mean, sigma])
                    g = np.linspace(0, len(radialprofile)-1, 10*len(radialprofile))
                    function = gaussianfit(g, *popt)
                    for x in range(len(function)):
                        if function[x] <= (1/2)*maximum:
                            fwhm = 2*x
                            if fwhm >=3 and fwhm <=30:
                                fwhm_list.append(fwhm)
                                good_stars_x.append(x_cent)
                                good_stars_y.append(y_cent)

                                plt.plot(f, radialprofile, 'b+:', label='data')
                                plt.plot(f, gaussianfit(f, *popt), 'ro:', label='fit')
                                plt.plot([0, fwhm/2], [(1/2)*maximum, (1/2)*maximum], 'g-.')
                                plt.plot([fwhm/2, fwhm/2], [0, (1/2)*maximum], 'g-.', label='HWHM')
                                plt.legend()
                                plt.xlabel('x position, HWHM = {}'.format(fwhm/2))
                                plt.ylabel('normalized counts')
                                plt.grid()
                                plt.savefig(r'/Users/owenalfaro/Desktop/testplots/GaussianPlot{}-{}.png'.format(name,a))
                                plt.close()
                                a += 1

                                break
                except RuntimeError:
                    #print("Could not find a Gaussian Fit...using whole pixel values to estimate fwhm")
                    '''
                for x in range(len(radialprofile)):
                    if radialprofile[x] <= (1/2)*maximum:
                        fwhm = 2*x
                        if fwhm >=3 and fwhm <=30:

                            plt.plot(f, radialprofile, 'b+:', label='data')
                            #plt.plot(f, gaussianfit(f, *popt), 'ro:', label='fit')
                            plt.plot([0, fwhm/2], [(1/2)*maximum, (1/2)*maximum], 'g-.')
                            plt.plot([fwhm/2, fwhm/2], [0, (1/2)*maximum], 'g-.', label='HWHM')
                            plt.legend()
                            plt.xlabel('x position, HWHM = {}'.format(fwhm/2))
                            plt.ylabel('normalized counts')
                            plt.grid()
                            plt.savefig(r'/Users/owenalfaro/Desktop/testplots/GaussianPlot{}-{}.png'.format(name,a))
                            a += 1
                            plt.close()
                            fwhm_list.append(fwhm)
                            good_stars_x.append(x_cent)
                            good_stars_y.append(y_cent)
                            break

            else:
                #print('Radial profile has length of 0...')
                continue
            if splitdata:
                return fwhm_list, good_stars_x, good_stars_y
        star_table = QTable([good_stars_x, good_stars_y], names=('starx', 'stary'))
        print('FWHM list: {}'.format(fwhm_list))
        return star_table

    def gaussianfit(self, x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
