from pylightcurve import Lightcurve
import pyfits
import datetime
import sunpy.time
from sunpy.time import parse_time, TimeRange, is_time_in_given_format
import math
from math import floor
import numpy as np
from numpy import nan
import pandas as pd
from pandas import DataFrame

class GOESLightcurve(Lightcurve):
    """
    This class is designed to simplify the handling of GOES data.
    """
    
    def __init__(self, start, end):
        uris = self.get_url_for_date_range(start, end)
        filepaths = self._download(uris)
        header, data = self._parse_fits(filepaths)
        meta = {}
        for column in data.columns.values.tolist():
            meta[column] = {}
        #    for entry in goes.header:
        #        meta[column][entry] = goes.meta[entry]

        start_i = data.index.searchsorted(parse_time(start))
        end_i = data.index.searchsorted(parse_time(end))
        self.import_data(data[start_i:end_i], meta)

    def _parse_fits(self, filepaths):
        """Parses a GOES FITS file from
        http://umbra.nascom.nasa.gov/goes/fits/"""
        
        data = pd.DataFrame()
        # We need to concatenate all of the FITS files.
        
        for filepath in filepaths:
            fits = pyfits.open(filepath)
            header = fits[0].header
            if len(fits) == 4:
                if is_time_in_given_format(fits[0].header['DATE-OBS'], '%d/%m/%Y'):
                    start_time = datetime.datetime.strptime(fits[0].header['DATE-OBS'], '%d/%m/%Y')
                elif is_time_in_given_format(fits[0].header['DATE-OBS'], '%d/%m/%y'):
                    start_time = datetime.datetime.strptime(fits[0].header['DATE-OBS'], '%d/%m/%y')
                else:
                    raise ValueError("Date not recognized")
                xrsb = fits[2].data['FLUX'][0][:, 0]
                xrsa = fits[2].data['FLUX'][0][:, 1]
                seconds_from_start = fits[2].data['TIME'][0]
            elif 1 <= len(fits) <= 3:
                start_time = parse_time(header['TIMEZERO'])
                seconds_from_start = fits[0].data[0]
                xrsb = fits[0].data[1]
                xrsa = fits[0].data[2]
            else:
                raise ValueError("Don't know how to parse this file")
    
            times = [start_time + datetime.timedelta(seconds=int(floor(s)),
                                                     microseconds=int((s - floor(s)) * 1e6)) for s in seconds_from_start]
    
            # remove bad values as defined in header comments
            
            xrsb[xrsb == -99999] = nan
            xrsa[xrsa == -99999] = nan
    
            # fix byte ordering
            newxrsa = xrsa.byteswap().newbyteorder()
            newxrsb = xrsb.byteswap().newbyteorder()
    
            file_data = DataFrame({'xrsa': newxrsa, 'xrsb': newxrsb}, index=times)
            file_data = file_data.replace(to_replace='0', value=np.nan)
            data = data.append(file_data)

        return header, data


    def get_goes_sat_num(self, start, end):
        """
        Parses the query time to determine which GOES satellite to use."""

        goes_operational = {
            2: TimeRange('1981-01-01', '1983-04-30'),
            5: TimeRange('1983-05-02', '1984-07-31'),
            6: TimeRange('1983-06-01', '1994-08-18'),
            7: TimeRange('1994-01-01', '1996-08-13'),
            8: TimeRange('1996-03-21', '2003-06-18'),
            9: TimeRange('1997-01-01', '1998-09-08'),
            10: TimeRange('1998-07-10', '2009-12-01'),
            11: TimeRange('2006-06-20', '2008-02-15'),
            12: TimeRange('2002-12-13', '2007-05-08'),
            13: TimeRange('2006-08-01', '2006-08-01'),
            14: TimeRange('2009-12-02', '2010-10-04'),
            15: TimeRange('2010-09-01', datetime.datetime.utcnow())}

        sat_list = []
        for sat_num in goes_operational:
            if ((start > goes_operational[sat_num].start() and
                 start < goes_operational[sat_num].end()) and
                (end > goes_operational[sat_num].start() and
                 end < goes_operational[sat_num].end())):
                    # if true then the satellite with sat_num is available
                    sat_list.append(sat_num)

        if not sat_list:
            # if no satellites were found then raise an exception
            raise Exception('No operational GOES satellites within time range')
        else:
            return sat_list

    def get_url_for_date_range(self, *args):
        """Returns a URL to the GOES data for the specified date.

        Parameters
        ----------
        args : TimeRange, datetimes, date strings
            Date range should be specified using a TimeRange, or start
            and end dates at datetime instances or date strings.
        satellite_number : int
            GOES satellite number (default = 15)
        data_type : string
            Data type to return for the particular GOES satellite. Supported
            types depend on the satellite number specified. (default = xrs_2s)
        """
        # TimeRange
        if len(args) == 1 and isinstance(args[0], TimeRange):
            start = args[0].start()
            end = args[0].end()
        elif len(args) == 2:
            start = parse_time(args[0])
            end = parse_time(args[1])
        if end < start:
            raise ValueError('start time > end time')

        length = (end - start)
        length = float(length.days) + float(length.seconds)/float(60*24*60) \
        + (float(length.microseconds)/float(1000*60*24))

        length = int(math.ceil(length))
        urls = []
        one_day = datetime.timedelta(days=1)
        for day in range(length):
            file_start = start + one_day * day
            file_end = start + one_day * (day+1)
            # find out which satellite and datatype to query from the query times
            sat_num = self.get_goes_sat_num(file_start, file_end)
            base_url = 'http://umbra.nascom.nasa.gov/goes/fits/'

            if start < parse_time('1999/01/15'):
                url = (base_url + "%s/go%02d%s.fits") % (file_start.strftime("%Y"),
                    sat_num[0], file_start.strftime("%y%m%d"))
            else:
                url = (base_url + "%s/go%02d%s.fits") % (file_start.strftime("%Y"),
                    sat_num[0], file_start.strftime("%Y%m%d"))

            urls.append(url)
        return urls

