from pylightcurve import Lightcurve
import pyfits
import datetime
from astropy.time import Time
import math
from math import floor
import numpy as np
from numpy import nan
import pandas as pd
from pandas import DataFrame

class KeplerLightcurve(object, Lightcurve):
    """
    This class is designed to simplify the handling of Kepler data.
    """

    quarters = {
        #0 : [datetime.datetime(2009,5,2,0,55,0), datetime.datetime(2009,5,11,17,51,30)],
        1 : [datetime.datetime(2009,5,13,0,15,0), datetime.datetime(2009,6,15,11,32,57)],
        2 : [datetime.datetime(2009,6,20,0,25,12), datetime.datetime(2009,9,16,23,9,27)],
        3 : [datetime.datetime(2009,9,18,17,19,58), datetime.datetime(2009,12,16,23,55,6)],
        4 : [datetime.datetime(2009,12,19,21,4,1), datetime.datetime(2010,3,19,16,53,28)]
        }

    def __init__(self, kic, start, end, title="Kepler Lightcurve"):
        times = [start, end]
        t = Time(times, scale='utc')
        uris = self.get_url_for_date_range(kic, start, end)
        fits_file = self._download(uris)
        data, meta = self._parse_fits(fits_file)
        start = data.index.searchsorted(t[0].datetime)
        end = data.index.searchsorted(t[1].datetime)
        self.import_data(data[start:end], meta)
        self.title = title
        #print self.data.index[0].value[:-5]
        self.kic = kic
        self.clc = np.array(self.data['pdcsap_flux_norm'])
        self.cts = (Time(self.data.index.astype(datetime.datetime), format='datetime',  scale='utc').jd - 2454833.0)*3600*24

    def identity_string(self):
        """
        Returns a string which identifies the lightcurve.

        Returns
        -------
        str
           An identifier of the light curve based on its length and KIC.
        """
        return str(self.kic) + "lc_len_" + str(len(self.clc))

    def add_data(self, curve=None, detrend=False, detrendmethod='none', nbins=101, order=3, knee=None, maxgap=1):
        """
        Add light curve data to the object.

        :: note 
           This method is on its way out! New methods have been implimented to store
           store the data better, and they will replace this one soon!

        Parameters
        ----------
        curvefile : string
           The file path file pointing to a light curve fits files.
        detrend : bool, optional, default: False
           A boolean flag which determines whether light curve should be detrended.
        detrendmethod : string, optional, default: 'none'
           The method for detrending the data. The options are 'savitzky_golay' to use
           :func:`.savitzky_golay`, 'runningmedian' to use :func:`.running_median`, or
           'highpassfilter' tp use :func:`.highpass_filter_lightcurve`
        nbins : int, optional, default: 101
           The width of the detrending window in bins of the light curve.
        order : int, optional, default: 3
           The polynomial order of the detrending filter.
        maxgap : int, optional, default+1
           The largest gap size (in bins) allowed before the light curve is deemed to contain gaps.

        Exceptions
        ----------
        NameError
           This needs to be replaced with an exception specific to the package!
           Error is raised if there is an I/O error while accessing a light curve file.
        """
        if curve == None:
            raise NameError("[Error] No light curve file given")

        try:
            dcurve = pyfits.open(curve)
        except IOError:
            raise NameError("[Error] An IO error occured when trying to access "+curve )
            mis_file = open('ioerror-files.log', 'a')
            mis_file.write(curve+'\n')
            mis_file.close()
            return
            
        if len(self.clc) == 0:
            self.id = dcurve[0].header['KEPLERID']
        elif dcurve[0].header['KEPLERID'] != self.id:
            raise NameError("Tried to add data from KIC"+str(dcurve[0].header['KEPLERID'])+" to KIC"+str(self.id))

        if dcurve[0].header['OBSMODE'] == 'long cadence':
            self.cadence = 'long'
        else:
            self.cadence = 'short'

        self.quarter = str(self.quarter)+str(dcurve[0].header['QUARTER'])

        # Assemble the new data into the class
        self.clc = np.append(self.clc, copy.deepcopy(dcurve[1].data['PDCSAP_FLUX']))
        self.cts = np.append(self.cts, copy.deepcopy(dcurve[1].data['TIME']*24*3600))
        self.cle = np.append(self.cle, copy.deepcopy(dcurve[1].data['PDCSAP_FLUX_ERR']))

        from astropy.time import Time
        times = dcurve[1].data['time']+2454833.0
        nans, x= self._nan_helper(times)
        times[nans]= np.interp(x(nans), x(~nans), times[~nans])
        times = Time(times, format='jd', scale='tcb')
        flux = dcurve[1].data['PDCSAP_FLUX']
        flux = flux / np.median(flux)
        flux = flux.byteswap().newbyteorder()
        file_data = DataFrame({'pdcsap_flux_norm': flux}, index=times.datetime)
        for each in fits[0].header:
            meta['pdcsap_flux_norm'][each] = fits[0].header[each]

        self.import_data(file_data, meta)
        
        dcurve.close()
        self.datagap = self.gap_checker(self.clc, maxgap=maxgap)
        self.interpolate()
        self.clc = self._dcoffset(self.clc) # remove a DC offset (calculated as the median of the light curve)
        self.data['pdcsap_flux_norm'] = self._dcoffset(self.data['pdcsap_flux_norm']) # remove a DC offset (calculated as the median of the light curve)
        if detrend:
            self.detrend(nbins, order)

        del dcurve
        
    def _parse_fits(self, filepaths):
        """
        Parses a Kepler FITS file from
        http://archive.stsci.edu/pub/kepler/lightcurves
        
        """
        
        data = pd.DataFrame()
        # We need to concatenate all of the FITS files.
        
        if not filepaths:
            pass
        
        for filepath in filepaths:
            fits = pyfits.open(filepath)
            from astropy.time import Time
            times = fits[1].data['time']+2454833.0
            nans, x= self._nan_helper(times)
            times[nans]= np.interp(x(nans), x(~nans), times[~nans])
            times = Time(times, format='jd', scale='tcb')
            
            flux = fits[1].data['PDCSAP_FLUX']
            flux = flux / np.median(flux)
            flux = flux.byteswap().newbyteorder()
            file_data = DataFrame({'pdcsap_flux_norm': flux}, index=times.datetime)
            data = data.append(file_data)
            meta = {}
            meta['pdcsap_flux_norm'] = {}
            for each in fits[0].header:
                meta['pdcsap_flux_norm'][each] = fits[0].header[each]

        return data, meta

    def get_kepler_quarter(self, start, end):
        """
        Takes a datetime object, and uses it to work out which 
        Kepler quarter it's in.
        """

        # The Kepler quarters, from page 9 of the Data Characteristics Handbook
        kepler_quarter = self.quarters

        quarter_list = []
        for q_num in kepler_quarter:
            if (
                (start <= kepler_quarter[q_num][1] or end <= kepler_quarter[q_num][1]) 
                and
                (start < kepler_quarter[q_num][1] and end > kepler_quarter[q_num][0] )
                ):
                    # if true then the satellite with sat_num is available
                    quarter_list.append(q_num)

        if not quarter_list:
            # if no satellites were found then raise an exception
            raise Exception("This time range doesn't fall within a period when Kepler data is available.")
        else:
            return quarter_list

    def get_url_for_date_range(self, kic, *args):
        """Returns a URL to the GOES data for the specified date.

        Parameters
        ----------
        kic : 'str'
            The Kepler input catalogue number for the star.
        args : TimeRange, datetimes, date strings
            Date range should be specified using a TimeRange, or start
            and end dates at datetime instances or date strings.
        """
        quarters = []
        # If only one argument, assume it's a quarter
        if len(args) == 1 and isinstance(args[0], int):
            quarters = [args[0]]
        if len(args) == 1 and isinstance(args[0], list):
            quarters = args[0]
        #If there are two arguments they'll be times
        elif len(args) == 2:
            start = Time(args[0], scale='utc').datetime
            end = Time(args[1], scale='utc').datetime
            if end < start:
                raise ValueError('start time > end time')

        urls = []
        # determine which quarter the dates are in
        if not quarters:
            quarters = self.get_kepler_quarter(start, end)
        base_url = 'http://archive.stsci.edu/pub/kepler/lightcurves/'
    
        kepler_quarter = self.quarters
    
        kic = kic.zfill(9)
        for quarter in quarters:
            if (kepler_quarter[quarter][1].timetuple().tm_mon >= 2 and \
                kepler_quarter[quarter][1].timetuple().tm_mon <= 10):
                hour = (kepler_quarter[quarter][1].timetuple().tm_hour-7)%24
            else:
                hour = (kepler_quarter[quarter][1].timetuple().tm_hour-8)%24
            url = (base_url + kic[0:4] + "/" + kic + "/" + "kplr%s-%s%s%s%s%s_llc.fits") % \
            (
                 kic,
                 kepler_quarter[quarter][1].year,
                 str(kepler_quarter[quarter][1].timetuple().tm_yday).zfill(3),
                 str(hour).zfill(2),
                 str(kepler_quarter[quarter][1].timetuple().tm_min).zfill(2),
                 str(kepler_quarter[quarter][1].timetuple().tm_sec).zfill(2)
                )
  
            urls.append(url)
        return urls

    def psd(self):
        """
        This function hard-wires the values for the column to be calculated by the
        psd() method from `Lightcurve`.
        This is done purely for simplicity, but it's also a boon for backwards
        compatibility.
        """
        sk, f = super(KeplerLightcurve, self).psd(column='pdcsap_flux_norm')
        return sk, f

    def gap_checker(self, d, maxgap=1):
        """
        Check for NaN gaps in the data greater than a given value.

        Parameters
        ----------
        d : :class:`numpy.ndarray`
           The array to check for gaps in the data.

        maxgap : int, optional, default: 1
           The maximum allowed size of gaps in the data.

        Returns
        -------
        bool
           ``True`` if there is a gap of maxgap or greater exists in ``d``, otherwise ``False``.
        """

        z = np.invert(np.isnan(d))
        y = np.diff(z.nonzero()[0])
        if len(y < maxgap+1) != len(y):
          return True
        else:
          return False

    def interpolate(self):
        """
        A method for interpolating the light curves, to compensate for NaN values.

        Examples
        --------

           >>> camelot = bf.Lightcurve(curves)
           >>> camelot.interpolate()

        """

        #for a in np.arange(len(self.lc)):
        z = self.clc
        nans, za= self._nan_helper(z)
        z[nans]= np.interp(za(nans), za(~nans), z[~nans]).astype('float32')
        self.clc = z