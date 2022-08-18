
################################################################################
################################################################################
from distutils.command.config import config
import cdsapi
import xarray as xr
import os
import numpy as np
from datetime import datetime
import pdb

class get_CAMS_data():
    def __init__(self,config_vars, sat_data):
        self.outfile_name = self.create_outfile_name(config_vars)
        self.model_start_date = max(sat_data.time).strftime('%Y-%m-%d')
        self.model_end_date = min(sat_data.time).strftime('%Y-%m-%d')

        self.time_list = [self.find_closest_model_time(time_stamp) for time_stamp in sat_data.time]
        self.unique_times = sorted(set(self.time_list))
        self.species = ['ozone', 'temperature']

        lat_max, lat_min, lon_max, lon_min = self.get_satellite_coord_extent(sat_data)
        self.north_lat = lat_max
        self.south_lat = lat_min
        self.west_lon = lon_min
        self.east_lon = lon_max

        self.download_model_data(config_vars)

        self.ozone_data = self.read_data(config_vars)

    def create_outfile_name(self,config_vars):
        """
        """
        # If the download path doesn't exist then make it
        if not os.path.exists(config_vars.model_download_path):
            os.mkdir(config_vars.model_download_path)

        # Create an outfile model name based on the time it was
        # created (ie make a timestamp)
        model_file_name =int(datetime.now().timestamp())

        config_vars.model_download_filename = os.path.join(
            config_vars.model_download_path,'{}_CAMS_model.grib'.format(str(model_file_name)))

    def download_model_data(self,config_vars):
        """
        """
        c = cdsapi.Client()

        outfile = config_vars.model_download_filename

        pressure_levs = ['1', '2', '3','5', '7', '10',
                    '20', '30', '50','70', '100', '150',
                    '200', '250', '300','400', '500', '600',
                    '700', '800', '850','900', '925', '950',
                    '1000']
        
        date_var = '{}/{}'.format(self.model_start_date,self.model_end_date)

        area_extent = [int(self.north_lat), int(self.west_lon), 
            int(self.south_lat), int(self.east_lon)]

        c.retrieve(
            'cams-global-reanalysis-eac4',
            {
                'date': date_var,
                'format': 'grib',
                'variable': self.species,
                'pressure_level': pressure_levs,
                'time': self.unique_times,
                'area': area_extent,
            },
            outfile)

    def find_closest_model_time(self,time_stamp):
        """
        Find the closest satellite time to the model times
        Model times are three hourly
        """
        model_time_list = np.arange(0,24,3)
        # Need to check about UTC vs time zones etc

        # Maybe add an exception here for time in wrong format
        decimal_min = time_stamp.minute / 60
        hour_decimal  = time_stamp.hour +  decimal_min
        nearest_model_ind = (np.abs(model_time_list - hour_decimal)).argmin()
        nearest_model_hour = model_time_list[nearest_model_ind]
        model_hour = str(nearest_model_hour).zfill(2) + ':00'

        return model_hour

    def get_satellite_coord_extent(self,sat_data):
        """
        """
        lat_max = max(sat_data.latitude)
        lat_min = min(sat_data.latitude)
        lon_max = max(sat_data.longitude)
        lon_min = min(sat_data.longitude)

        # Round up or down to nearest integer for latitude
        lat_max = np.ceil(lat_max)
        lat_min = np.floor(lat_min)

        # Do the same for longitude - will currently be inefficient when crossing
        # the international date line as it will download the whole globe.
        # This can be addressed later

        lon_max = np.ceil(lon_max)
        lon_min = np.floor(lon_min)

        return lat_max, lat_min, lon_max, lon_min


    def read_data(self,config_vars):
        """
        """
        fname = config_vars.model_download_filename
        data = xr.open_dataset(fname, engine='cfgrib')
 
        return data


    def delete_file(self,config_vars):
        """
            Delete the CAMS model file at given file path
        """
        if config_vars.verbose:
            print("Deleting model file {}".format(config_vars.model_download_filename))
        os.remove(config_vars.model_download_filename)
        
