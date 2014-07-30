from pylightcurve import Lightcurve
import pyfits
import datetime
import sunpy.time
from sunpy.time import parse_time, TimeRange, is_time_in_given_format
from sunpy.net import hek
import math
from math import floor
import numpy as np
from numpy import nan
import pandas as pd
from pandas import DataFrame
import dateutil.relativedelta as rd
from netCDF4 import *
from dateutil import parser
from math import floor, log10
import copy

class GOESLightcurve(Lightcurve):
    """
    This class is designed to simplify the handling of GOES data.
    """
    title="Title"
    detrended=False
    tags = {}
    #detrend_method
    def __init__(self, start, end, **kwargs):
        if "title" in kwargs:
            self.title=kwargs["title"]
        if "cadence" in kwargs:
            self.cadence = kwargs['cadence']
        else:
            self.cadence = "2sec"

        if self.cadence == "2sec":
            uris = self.get_url_for_date_range(start, end)
            filepaths = self._download(uris)
            header, data = self._parse_fits(filepaths)
        elif self.cadence == "1min":
            uris = self.get_url_for_date_range_noaa(start, end)
            filepaths = self._download(uris)
            header, data = self._parse_cdf(filepaths)
        meta = {}
        for column in data.columns.values.tolist():
            meta[column] = {}
        #    for entry in goes.header:
        #        meta[column][entry] = goes.meta[entry]

        start_i = data.index.searchsorted(parse_time(start))
        end_i = data.index.searchsorted(parse_time(end))
        self.import_data(data[start_i:end_i], meta)

        self._get_hek_flares(start, end)
        
        self.cts = self.time_seconds()
        if "default" in kwargs:
            self.default = kwargs["default"]
            self.clc=np.array(self.data[kwargs["default"]])

    def __getitem__(self, key):
        curve = copy.deepcopy(self)
        curve.data = curve.data[key]
        curve.highlight = curve.highlight[key]
        curve.cts = curve.time_seconds()
        if self.default:
            curve.clc=np.array(curve.data[self.default])
        return curve
            
    def _parse_cdf(self, filepaths):
        """
        Parse a GOES CDF file from
        http://satdat.ngdc.noaa.gov/sem/goes/data/new_avg
        """
        data = pd.DataFrame()
        # Concatenate the CDF Files.
        # We need to work on each of the filepaths.
        header = {}
        for filepath in filepaths:
            f = Dataset(filepath, 'r', format='NETCDF4')

            aavg = f.variables['A_AVG']
            xrsa = np.array(aavg[:])
            
            bavg = f.variables['B_AVG']
            xrsb = np.array(bavg[:])

            xrsb[xrsb == -99999] = nan
            xrsa[xrsa == -99999] = nan

            newxrsa = xrsa.byteswap().newbyteorder()
            newxrsb = xrsb.byteswap().newbyteorder()

            seconds_from_start = f.variables['time_tag']
            seconds_from_start = np.array(seconds_from_start[:])

            times = [datetime.utcfromtimestamp(seconds_from_start[s]/1000) 
                     for s in range(len(seconds_from_start))]
            
            file_data = pd.DataFrame({'xrsa': newxrsa,
                          'xrsb': newxrsb}, index=times)
            file_data = file_data.replace(to_replace=aavg.missing_value, value=np.nan)

            # Parse the header information for each of the columns

            header['xrsa'] = {}
            for i in range(len(aavg.ncattrs())):
                name = aavg.ncattrs()[i]
                header['xrsa'][name] = aavg.getncattr(name)

            header['xrsb'] = {}
            for i in range(len(bavg.ncattrs())):
                name = bavg.ncattrs()[i]
                header['xrsb'][name] = bavg.getncattr(name)

            data = data.append(file_data)

        return header, data

    def _get_hek_flares(self, start, end):
        """
        A function to add SSW Latest Events tags to the light curve.
        """
        client = hek.HEKClient()
        event_type = 'FL'
        result = client.query(hek.attrs.Time(start,end),
                              hek.attrs.EventType(event_type))
        tag_list = []
        if len(result) != 0:
            for i in range(len(result)):
                start = parser.parse(result[i]['event_starttime'])
                end = parser.parse(result[i]['event_endtime'])
                peak = parser.parse(result[i]['event_peaktime'])
                classific = result[i]['fl_goescls']
                method =  result[i]['frm_name']
                tag_list.append({'start_time':start, 
                                 'end_time':end, 
                                 'peak_time':peak,
                                 'method':method,
                                 'class':classific  })

            tags = pd.DataFrame(tag_list)
            self.tags['SSW_Latest_Events'] = tags[tags.method=="SSW Latest Events"].reset_index(drop=True)

        return None

    def import_tags(self, tags, name):
        """
        This method takes a DataFrame with tags, for example flare locations, and
        adds them to the tag store for the lightcurve.

        Parameters
        ----------
        tags : Pandas DataFrame
           A Pandas DataFrame containing the tags.

        name : The name of the tag category.
        """

        self.tags[name] = tags
        return self


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
                    start_time = datetime.strptime(fits[0].header['DATE-OBS'], '%d/%m/%Y')
                elif is_time_in_given_format(fits[0].header['DATE-OBS'], '%d/%m/%y'):
                    start_time = datetime.strptime(fits[0].header['DATE-OBS'], '%d/%m/%y')
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
    
            times = [start_time + timedelta(seconds=int(floor(s)),
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
            15: TimeRange('2010-09-01', datetime.utcnow())}

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

    def get_url_for_date_range_noaa(self, *args):
        """
        Returns a URL to the GOES data for the specified date,
        from the NOAA source.

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

        length = rd.relativedelta(end, start)
        length = float(length.months) + float(length.days)/float(31) 
        length = int(math.ceil(length))

        days_back = start.day
        start = start + rd.relativedelta(days=-(days_back-1))
        urls = []
        for month in range(length):
            file_start = start + rd.relativedelta(months = month ) 
            file_end = start + rd.relativedelta(months = month+1) + rd.relativedelta(days=-1)
            
            # find out which satellite and datatype to query from the query times
            sat_num = self.get_goes_sat_num(file_start, file_end)
            base_url = 'http://satdat.ngdc.noaa.gov/sem/goes/data/new_avg/'

            url = (base_url + "%s/%s/goes%02d/netcdf/g%02d_xrs_1m_%s_%s.nc") % (file_start.strftime("%Y"),
                                                                                file_start.strftime("%m"),
                                                                                sat_num[0], sat_num[0],
                                                                                file_start.strftime("%Y%m%d"),
                                                                                file_end.strftime("%Y%m%d")
            )
            urls.append(url)
        return urls

            
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
        one_day = timedelta(days=1)
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

    def _add_line(self, axis, location, label, color='green'):
        ymin, ymax = axis.get_ybound()
        xmin, xmax = axis.get_xbound()

        axis.axvline(location, color=color, linestyle='--', linewidth=1)
        axis.text(location,ymax,
             label,
             rotation=90, backgroundcolor='white',
             color=color, size='smaller',
             verticalalignment='center',horizontalalignment='center',
        bbox=dict(facecolor='white',
                  edgecolor=color,
                  boxstyle='round'))

        return None

    def plot_tags(self, ax, column, labelcolumn):
        for tag in self.tags:
            tags = self.tags[tag]
            color = ax._get_lines.color_cycle
            ccolor=next(color)
            for i in range(len(tags)):
                flare = tags.ix[i]
                peak = flare[column]
                text = str(flare[labelcolumn])
                self._add_line(ax, peak, text, color=ccolor)

    def goes_class(self, energy):
        energye = int(floor(log10(energy)))
        energies = np.array([-9, -8, -7, -6, -5, -4, -3, -2])
        categories = ['Below A', 'A', 'B', 'C', 'M', 'X', 'Above X']
        subcat= floor(10*(energy/10**floor(log10(energy))))/10
        return categories[np.where(energies==energye)[0][0]]+str(subcat)