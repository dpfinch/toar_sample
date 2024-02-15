################################################################################
################################################################################
import numpy as np
import pdb
import xarray as xr
import netCDF4 as nc4


class get_TROPESS_data():
    def __init__(self, config_vars,sat_data) -> None:
        # self.outfile_name = self.create_outfile_name(config_vars)
        self.model_start_date = min(sat_data.time)
        self.model_end_date = max(sat_data.time)

        self.time_list = [self.find_closest_model_time(time_stamp) for time_stamp in sat_data.time]
        self.unique_times = sorted(set(self.time_list))
        self.species = ['ozone']

        lat_max, lat_min, lon_max, lon_min = self.get_satellite_coord_extent(sat_data)
        self.north_lat = lat_max
        self.south_lat = lat_min
        self.west_lon = lon_min
        self.east_lon = lon_max
        self.ozone_data = self.read_data(config_vars)

    def read_data(self,config_vars):
        
        opendap_base_url = 'https://tropess.gesdisc.eosdis.nasa.gov/opendap/TCR2_6HR_VERTCONCS/TRPSCRO36H3D.1/'
        year_to_get = self.model_start_date.year
        year_to_get = 2008
        filename = 'TROPESS_reanalysis_6hr_o3_{}.nc'.format(year_to_get)

        full_url = '{}{}'.format(opendap_base_url,filename)
        dataset = nc4.Dataset(full_url)
        # model_data = xr.open_dataset(full_url)
        # idx = (np.abs(array - value)).nanargmin()

        pdb.set_trace()

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
    
if __name__ == '__main__':
    url = 'https://tropess.gesdisc.eosdis.nasa.gov/opendap/TCR2_6HR_VERTCONCS/TRPSCRO36H3D.1/TROPESS_reanalysis_6hr_o3_2005.nc'
    dataset = nc4.Dataset(url)
    pdb.set_trace()