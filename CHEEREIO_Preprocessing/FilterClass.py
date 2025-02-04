import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import os.path
from datetime import datetime,timedelta
import copy
import gcpy
import xesmf as xe

def daily_iterator(start_date, end_date):
    current_date = start_date
    while current_date <= end_date:
        yield current_date  # Yield the current date for processing
        current_date += timedelta(days=1)

def custom_L2L_regridder(
        llres_in,
        llres_out,
        weightsdir='.',
        reuse_weights=False,
        in_extent=[-180, 180, -90, 90],
        out_extent=[-180, 180, -90, 90],
        use_method='bilinear'
    ):
    """
    Create an xESMF regridder between two lat/lon grids built into the filter object

    This is based on GCpy's make_regridder_L2L() routine, but allows for customization of the method used

    Args:
        llres_in: str
            Resolution of input grid in format 'latxlon', e.g. '4x5'
        llres_out: str
            Resolution of output grid in format 'latxlon', e.g. '4x5'

    Keyword Args (optional):
        weightsdir: str
            Directory in which to create xESMF regridder NetCDF files
            Default value: '.'
        reuse_weights: bool
            Set this flag to True to reuse existing xESMF regridder NetCDF files
            Default value: False
        in_extent: list[float, float, float, float]
            Describes minimum and maximum latitude and longitude of input grid
            in the format [minlon, maxlon, minlat, maxlat]
            Default value: [-180, 180, -90, 90]
        out_extent: list[float, float, float, float]
            Desired minimum and maximum latitude and longitude of output grid
            in the format [minlon, maxlon, minlat, maxlat]
            Default value: [-180, 180, -90, 90]
        use_method: str
            One of the acceptable ESMF regridding methods:
            bilinear, conservative, conservative_normed, patch, nearest_s2d, nearest_d2s
            See https://xesmf.readthedocs.io/en/stable/notebooks/Compare_algorithms.html for details
    """
    llgrid_in = gcpy.make_grid_LL(llres=llres_in, in_extent=in_extent, out_extent=out_extent)
    llgrid_out = gcpy.make_grid_LL(llres=llres_out, out_extent=[-180, 180, -90, 90])

    if in_extent == [-180, 180, -90, 90] and out_extent == [-180, 180, -90, 90]:
        weightsfile = os.path.join(
            weightsdir, '{}_{}_{}.nc'.format(
                use_method, llres_in, llres_out))
    else:
        in_extent_str = str(in_extent).replace(
            '[', '').replace(
            ']', '').replace(
            ', ', 'x')
        out_extent_str = str(out_extent).replace(
            '[', '').replace(
            ']', '').replace(
            ', ', 'x')
        weightsfile = os.path.join(
            weightsdir, '{}_{}_{}_{}_{}.nc'.format(
                use_method, llres_in, llres_out, in_extent_str, out_extent_str))

    try:
        regridder = xe.Regridder(llgrid_in,
                                 llgrid_out,
                                 method=use_method,
                                 filename=weightsfile,
                                 reuse_weights=False)
    except BaseException:
        regridder = xe.Regridder(llgrid_in,
                                 llgrid_out,
                                 method=use_method,
                                 filename=weightsfile,
                                 reuse_weights=False)

    return regridder
        
class MapFilter():
    def __init__(self, MaskToNaNs=False, MERRA2_CN_fl = "/home/milletd/shared/data/geos-chem/GEOS_0.5x0.625/MERRA2/2015/01/MERRA2.20150101.CN.05x0625.nc4"):
        # The init function for this class reads in the MERRA2 landmask file and removes the unnecessary 'time' dimension
        MERRA2_CN = xr.load_dataset(MERRA2_CN_fl)
        self.landmask = MERRA2_CN['FRLAND'].squeeze('time')
        if(MaskToNaNs):
            self.landmask = self.landmask.where(self.landmask != 0, np.nan)
    
    def isolate_YYYYMMDD(self, filenames):
        """
        Isolates the YYYYMMDD component from a list of filenames in the format '{sourcedir}/YYYYMMDD_{speciesname}_Column.nc'.

        Args:
            filenames [list]: A list of filenames.

        Returns:
            A list of YYYYMMDD datetimes.
          """

        YYYYMMDD_list = []
        for filename in filenames:
            # Split the filename into its components.
            parts = os.path.basename(filename).split("_")

            # Extract the YYYYMMDD component.
            #YYYYMMDD = parts[0] # - gotta convert this into a datetime.datetime
            # in order to compare its output to timeperiod
            YYYYMMDD = datetime.strptime(parts[0], '%Y%m%d')

            # Append the YYYYMMDD to the list.
            YYYYMMDD_list.append(YYYYMMDD)

        return YYYYMMDD_list
    
    def isolate_YYYYMM(self, filenames):
        """
        Isolates the YYYYMM component from a list of filenames in the format '{sourcedir}/YYYYMM_{speciesname}_screened.nc'.
        Relevant only for Kelley's original isoprene files

        Args:
            filenames [list]: A list of filenames.

        Returns:
            A list of YYYYMM datetimes.
          """

        YYYYMM_list = []
        for filename in filenames:
            # Split the filename into its components.
            parts = os.path.basename(filename).split("_")

            # Extract the YYYYMMDD component.
            #YYYYMMDD = parts[0] # - gotta convert this into a datetime.datetime
            # in order to compare its output to timeperiod
            YYYYMM = datetime.strptime(parts[0], '%Y%m')

            # Append the YYYYMMDD to the list.
            YYYYMM_list.append(YYYYMM)

        return YYYYMM_list
        
     
    def return_window(self, date_targeted, window_width=0):
        """
        For a given date, this method identifies a widow around that date.
        
        Args:
        date_targeted (datetime): The target date.
        window_width (int): Half the width of the window in days.

        Returns:
            A list of two datetimes corresponding to the start and end of the relevant window.
        """
        start_date = date_targeted - timedelta(days=window_width)
        end_date = date_targeted + timedelta(days=(window_width+1)) # Fencepost problem fix!
        
        return [start_date, end_date]
        
    
    def globObs(self, longspeciesname, timeperiod, sourcedir="/users/5/brewe222/CrIS_Column_TATM_fix/"):
        """
        For a given species and set of dates, this function returns all relevant CrIS column files.
        By default, it looks in Kelley's CrIS_Column_TATM_fix directory, but that can be changed.
        
        Args:
            longspeciesname (str): the full name of the species in question (e.g. "Ethane", "Isoprene").
                See Kelley's CrIS_Column_TATM_fix directory for examples
            timeperiod (list of datetimes): a list of two datetime.datetimes to select for opening with xarray
            sourcedir (string): where to look for the filenames
        Returns:
            A list of relevant files.
          """
        
        # If we're using Kelley's original Isoprene files, the naming conventions are different
        # And the files are monthly with daily data contained within, instead of daily
        if longspeciesname == "Isoprene":
            # First convert the dates to YYYYMM datetimes
            # Have to modify the window because isoprene files are formatted differently
            modified_window = copy.deepcopy(window)
            modified_window[1] = window[1] - timedelta(days = 1)
            tp_yyyymm = [datetime(date.year, date.month, 1) for date in modified_window]
            
            obs_list = glob(f'{sourcedir}/*_daily_CrIS_Isoprene_screened.nc')
            obs_list.sort()
            filtered_obs = [obs for obs, obsdate in zip(obs_list, self.isolate_YYYYMM(self, obs_list))
                            if obsdate >= tp_yyyymm[0] and obsdate <= tp_yyyymm[1]]
            
        else:
            obs_list = glob(f'{sourcedir}/*_{longspeciesname}_Column.nc')
            obs_list.sort()
            filtered_obs = [obs for obs, obsdate in zip(obs_list, self.isolate_YYYYMMDD(self, obs_list))
                            if obsdate >= timeperiod[0] and obsdate < timeperiod[1]]
        return filtered_obs
    
    
    def output_column(self, relevant_files, use_temp_mask=True):
        """
        Read in all relevant CrIS column files. If multiple files are specified, average the column output together.
        
        Args:
            relevant_files (list): a list of filenames to read in. Intended to be the output of GlobObs.

        Returns:
            A DataArray of Column Values from CrIS, with dimensions (P90,lat,lon,ANN).
          """
        if(len(relevant_files)>1):
            print("====================================================================================================")
            print("Multiple files supplied. Columns represent {} days worth of data".format(len(relevant_files)))
            print("====================================================================================================")
        elif(len(relevant_files)==1):
            print("====================================================================================================")
            print("Single files supplied. Columns represent {}.".format(relevant_files[0]))
            print("====================================================================================================")
        else:
            print("====================================================================================================")
            print("ERROR: NO RELEVANT FILES SUPPLIED, TRY AGAIN")
            print("====================================================================================================")
        # Open all files as a combined dataset
        datasets = [xr.open_dataset(file) for file in relevant_files]
        combined = xr.concat(datasets, dim="time").rename({"time": "day"})  # Concatenate along a new "day" dimension
        # return only the required information         
        if(use_temp_mask):
            column_array = combined['Column']
            TempMask = xr.where(combined['Tsfc'] >273.15, 1, np.nan).expand_dims(P90=column_array.coords['P90'],
                                                                                 ANN=column_array.coords['ANN'])
            combined['Column'] = column_array*TempMask

        return(combined)
    
    def output_isop_column(self, relevant_files, window):
        """
        Because Kelley's original ROCRv1 Isoprene files are formatted differently from the ROCRv2 files,
        This has to work extremely differently from output_column.
        
        Isoprene files are monthly, and contain daily data, but have no 'datetime' object associated with them.
        Thus, this formulation needs to be rather more convoluted than 'output_column'.
        
        In particular, there are very different processes required if the window in question goes over a month
        barrier.
        
        Args:
            relevant_files (list): a list of filenames to read in. Intended to be the output of GlobObs.
            window (list): a list of datetimes to average data in between.

        Returns:
        A DataArray of Column Values from CrIS, with dimensions (lat,lon).
        """
        print("====================================================================================================")
        print("This is an isoprene case and requires special handling. See Code for details.")
        print("====================================================================================================")
     
        
        # Have to modify the window because isoprene files are formatted differently
        # and because xarray's slice is inclusive, unlike everything else in python
        modified_window = copy.deepcopy(window)
        modified_window[1] = window[1] - timedelta(days = 1)
        tp_yyyymm = [datetime(date.year, date.month, 1) for date in modified_window]
        days = [pd.to_numeric(date.strftime('%d')) for date in modified_window]            
        
        if(len(relevant_files) == 1):
            # This is the simple case - just take the slice across those days and average it
            # *NOTE THAT SLICE IS INCLUSIVE AND days ARE 1-INDEXED BECAUSE THE REAL WORLD ISN'T PYTHON*
            print("====================================================================================================")
            print("Single files supplied. Columns represent {}.".format(relevant_files[0]))
            print("====================================================================================================")
            isopdata = xr.open_dataset(relevant_files[0], engine='netcdf4')
            # Need a special case handling here, too, to handle window sizes of 1
            if days[0] == days[1]:
                output = isopdata.sel(day=days[0]).expand_dims('day', axis=0)
            else:
                output = isopdata.sel(day=slice(days[0]-1, days[1]-1))#.mean(dim="day")
            
        elif(len(relevant_files)==2):
            # This is the complex case, where values go across months
            print("====================================================================================================")
            print("Two files supplied. Will combine across months")
            print("====================================================================================================")
            print("Days: {}".format(days))
            datasets = [xr.open_dataset(file) for file in relevant_files]
            # Select relevant days from each month
            data_month1 = datasets[0].sel(day=slice(days[0]-1, None))
            data_month2 = datasets[1].sel(day=slice(0, days[1]))
            
            # Concatenate the selected data along the day dimension
            combined_data = xr.concat([data_month1, data_month2], dim="day")

            # Compute the average over the day dimension
            output = combined_data#.mean(dim="day")

        elif(len(relevant_files) > 2):
            # I suppose if you made a big enough window, this could be possible, but it's not why I'm making this.
            print("====================================================================================================")
            raise Exception("More than two files supplied. A window this wide hasn't been contemplated for this program")
            print("====================================================================================================")
    
            
        else:
            print("====================================================================================================")
            raise Exception("ERROR: NO RELEVANT FILES SUPPLIED, TRY AGAIN")
            print("====================================================================================================")

        return output['Isop']
        
    def run_averaging(self, combined_data, window_width):
        """
        This is where the actual averaging gets performed.
        If we decide we want to make averaging more complicated, this is where we change it.
        
        Args:
            combined_data (DataArray): A data array of column values that certainly has dimensions (lat, lon, day) 
            and may or may not include (P90 and ANN).
            window_width (int): the number of days on either side 

        Returns:
            A DataArray of Column Values from CrIS averaged across the date range, preserving all other dims.
        """

        output = copy.deepcopy(combined_data)
        
        Averaged_Column = combined_data['Column'].mean(dim="day", skipna=True)
        output['Column']=Averaged_Column
        # for the non-Column data, just keep the data from the center of the window
        # Because of 0-indexing, this is just the value of window_width
        output=output.isel(day=window_width).drop_vars("Date")
            
        return output
    
    def get_column(self, longspeciesname, date_targeted, window_width=0,
                   sourcedir="/users/5/brewe222/CrIS_Column_TATM_fix/", run_averaging=True):
        """
        This is a master method to run an extraction
        
        Args:
        longspeciesname (str): the full name of the species in question (e.g. "Ethane", "Isoprene").
                See Kelley's CrIS_Column_TATM_fix directory for examples
        date_targeted (datetime): The target date.
        window_width (int): Half the width of the window in days.

        Returns:
            A Column map with the correct dimensions (P90, lat, lon, ANN)
            If window_width is specified, then the data will be averaged across the timewindow.
        """
        # First identify the window to run over
        window = self.return_window(self, date_targeted, window_width)
        
        # Next, return the relevant files for that window
        relevant_files = self.globObs(self, longspeciesname, window, sourcedir)
        
        if longspeciesname == "Isoprene":
            Column = self.output_isop_column(self, relevant_files, window)

        else:
            # Now open the relevant files
            Column = self.output_column(self, relevant_files)
        
        # This can be run either with or without the averaging applied
        # If run in this mode with window_width>=0, average the columns together
        if run_averaging:
            Averaged_Column = self.run_averaging(self, Column, window_width)
            return Averaged_Column
        else:
            return Column
        
    
    def run_Landfilter(self, column_array):
        """
        Given a DataArray of column values with P90 and ANN dimensions
        this routine expands the fracland map it read in in init to those dimensions
        and then multiplies the two together
        
        Args:
        column_array (DataArray): a data array containing a set of column values 

        Returns:
            A list of two datetimes corresponding to the start and end of the relevant window.
        """
        aligned_landmask = self.landmask.copy()
        
        aligned_landmask = aligned_landmask.assign_coords(
            lat=column_array.coords['lat']
        )
        
        # Handle the Isoprene no-P90 case and the regular P90 case
        if('P90' in column_array.dims):
            expanded_fracland = aligned_landmask.expand_dims(P90=column_array.coords['P90'],
                                                          ANN=column_array.coords['ANN'])
        else:
            expanded_fracland = aligned_landmask.expand_dims(ANN=column_array.coords['ANN'])
        
        # Apply the landmask
        masked_column = column_array * expanded_fracland
        
        # return masked_column
        return masked_column

    def Regrid_CrIS(self, regridder, dataset_to_regrid, targetdatetime, 
                    llres_out, extent=[-180, 180, -90, 90]):
        """
        This routine is a wrapper to run the regridding for the relevant CrIS files.
        It requires an already created regridder, so we don't rebuild that every time.

        The last lines fix the dimensions so that the output file represents the desired day.

        Args:
            regridder: regridder object
                A regridder built with MapFilter.custom_L2L_regridder
                Is preset for a given resolution and regridding algorithm
                (bilinear by default)
            dataset_to_regrid: xarray DataSet
                A CrIS DataSet in Xarray, post averaging and landmasking.
            targetdatetime: datetime
                A datetime object for the day the averaged and regridded file should represent
                the output of datetime.strptime(targetdate, "%Y%m%d")
                This is necessary here because the regridding adds it back in. 
            llres_out: str
                Resolution of output grid in format 'latxlon', e.g. '4x5'
        
        Keyword Args (optional):
            out_extent: list[float, float, float, float]
                Desired minimum and maximum latitude and longitude of output grid
                in the format [minlon, maxlon, minlat, maxlat]
                Default value: [-180, 180, -90, 90]
        """


        targetdate_np = np.datetime64(targetdatetime.strftime("%Y-%m-%d"), "s")

        output = xr.apply_ufunc(
            regridder, 
            dataset_to_regrid,  
            input_core_dims=[["lat", "lon"]],  # Specify core dimensions to apply regridding
            output_core_dims=[["lat_new", "lon_new"]],  # Output core dims (91x144)
            exclude_dims=set(("lat", "lon")),  # These dimensions will change in the output
            vectorize=True,  # Apply function across all non-core dimensions
        ).rename({"lat_new": "lat", "lon_new": "lon"})

        OutGrid = gcpy.make_grid_LL(llres=llres_out, out_extent=extent)
        
        output = output.assign_coords(lat=OutGrid['lat'], lon=OutGrid['lon'])
        output = output.assign_coords(time=targetdate_np)

        return output
