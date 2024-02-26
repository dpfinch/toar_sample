################################################################################
################################################################################
import numpy as np
import pdb
import xarray as xr
import netCDF4 as nc4
from datetime import datetime, timedelta

class get_TROPESS_data():
    def __init__(self, config_vars,sat_data) -> None:
        # self.outfile_name = self.create_outfile_name(config_vars)
        self.sat_start_date = min(sat_data.time)
        self.sat_end_date = max(sat_data.time)

        self.time_list = [self.find_closest_model_time(time_stamp, config_vars) for time_stamp in sat_data.time]
        self.unique_times = sorted(set(self.time_list))
        self.species = 'o3'

        lat_max, lat_min, lon_max, lon_min = self.get_satellite_coord_extent(sat_data)
        self.north_lat = lat_max
        self.south_lat = lat_min
        self.west_lon = lon_min
        self.east_lon = lon_max
        self.ozone_data = self.read_data(config_vars)

    def read_data(self,config_vars):
        
        opendap_base_url = 'https://tropess.gesdisc.eosdis.nasa.gov/opendap/TCR2_6HR_VERTCONCS/TRPSCRO36H3D.1/'
        year_to_get = self.sat_start_date.year

        filename = 'TROPESS_reanalysis_6hr_o3_{}.nc'.format(year_to_get)

        full_url = '{}{}'.format(opendap_base_url,filename)
        dataset = nc4.Dataset(full_url)

        # o3 dataset (1464, 27,160,360) (time, level, lat, lon)
        model_lat = dataset.variables['lat'][:].data
        model_lon = dataset.variables['lon'][:].data
        model_time = dataset.variables['time'][:].data
        #NASA model file says time is in units of seconds since 2005-01-01 00:00:00
        # Seems not to be true - its just seconds since the begining of that year
        model_times = [datetime(year_to_get,1,1) + timedelta(seconds = int(x)) for x in model_time]

        model_lat_min_ind = np.nanargmin((np.abs(model_lat - self.south_lat)))
        model_lat_max_ind = np.nanargmin((np.abs(model_lat - self.north_lat)))
        model_lon_min_ind = np.nanargmin((np.abs(model_lon - self.west_lon)))
        model_lon_max_ind = np.nanargmin((np.abs(model_lon - self.east_lon)))
        
        sat_start_time_seconds =  self.sat_start_date - datetime(year_to_get,1,1) 
        sat_start_time_seconds = sat_start_time_seconds.total_seconds()
        sat_end_time_seconds = self.sat_end_date -  datetime(year_to_get,1,1) 
        sat_end_time_seconds = sat_end_time_seconds.total_seconds()
        model_time_min_ind = np.nanargmin((np.abs(model_time - sat_start_time_seconds)))
        model_time_max_ind = np.nanargmin((np.abs(model_time - sat_end_time_seconds))) 
       
        model_data = dataset.variables[self.species][model_time_min_ind:model_time_max_ind,:,model_lat_min_ind:model_lat_max_ind,model_lon_min_ind:model_lon_max_ind]
        # idx = (np.abs(array - value)).nanargmin()

        return model_data

    def find_closest_model_time(self,time_stamp, config_vars):
        """
        Find the closest satellite time to the model times
        Model times are three hourly for CAMS, six for NASA
        """

        if config_vars.model_type == 'NASA_model':
            model_time_list = np.arange(0,24,6)
        else:
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


    
if __name__ == '__main__':
    url = 'https://tropess.gesdisc.eosdis.nasa.gov/opendap/TCR2_6HR_VERTCONCS/TRPSCRO36H3D.1/TROPESS_reanalysis_6hr_o3_2005.nc'
    dataset = nc4.Dataset(url)
    
    pdb.set_trace()

    