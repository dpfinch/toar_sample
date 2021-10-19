################################################################################
################################################################################
from create_config_vars import config_vars
import get_satellite_data
from file_writer import output_to_file
import utils
import netCDF4 as nc
import sys
import os
from urllib.parse import urljoin
from requests import get
from requests.exceptions import HTTPError
import pdb
import pandas as pd
import numpy as np
import pathlib
from scipy.interpolate import interp1d
#  This prevents a userwarning about converting masked arrays to nans
import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)

class get_model_data():
    def __init__(self, config_vars, filename_format, sat_year):

        self.filepath_format = os.path.join(config_vars.model_download_path,
                                            filename_format)

        model_out_data = self.extract_model_data(sat_year)
        self.latitude = model_out_data[0]
        self.longitude = model_out_data[1]
        self.levels = model_out_data[2]
        self.datetime = model_out_data[3]
        self.model_o3 = model_out_data[4]
        self.model_t = model_out_data[5]
        self.model_ps = model_out_data[6]

    def extract_model_data(self, sat_year):
        """
            Info here ...
        """
        model_o3_filename = self.filepath_format.format('o3', sat_year)
        model_t_filename = self.filepath_format.format('t', sat_year)
        model_ps_filename = self.filepath_format.format('ps', sat_year)

        for file_test in [model_o3_filename, model_t_filename, model_ps_filename]:
            if not os.path.isfile(file_test):
                print("--> Unable to find model file {}".format(file_test))
                return None, None, None, None, None, None, None

        model_o3_dataset = nc.Dataset(model_o3_filename,'r')
        model_latitude = model_o3_dataset.variables['lat'][:]
        model_longitude = model_o3_dataset.variables['lon'][:]
        # Check if the model is 0 - 360 and if so, make it -180 - 180
        if max(model_longitude) > 180:
            model_longitude = model_longitude - 180
        if config_vars.verbose:
            print("--> Adding 0.01 hPa to top of model levels as TOA")
        model_levels = list(model_o3_dataset.variables['lev'][:])
        model_levels.append(0.01)
        model_time = model_o3_dataset.variables['time']
        model_time_unit = model_time.units
        model_o3 = model_o3_dataset.variables['o3'][:] # in ppb
        # Convert to mixing ratio
        model_o3 = model_o3 * 1e-9

        model_dt = utils.days_since_to_dt(model_time, model_time_unit)

        model_o3_dataset.close()
        # Get the temperature data
        model_t_dataset = nc.Dataset(model_t_filename,'r')
        model_t = model_t_dataset.variables['t'][:]
        model_t_dataset.close()
        # Get the surface pressure data
        model_ps_dataset = nc.Dataset(model_ps_filename,'r')
        model_ps = model_ps_dataset.variables['ps'][:]
        model_ps_dataset.close()

        return (model_latitude, model_longitude, model_levels, model_dt,
                model_o3, model_t, model_ps)


def download_model_montly_data(config_vars,model_filepath):
    '''
        info here...
    '''
    nc_filename =  pathlib.PurePath(model_filepath).name
    full_url = urljoin(config_vars.model_repo_url, nc_filename)

    # If the download path doesn't exist then make it
    if not os.path.exists(config_vars.model_download_path):
        os.mkdir(config_vars.model_download_path)

    out_fname = os.path.join(config_vars.model_download_path, nc_filename)

    if config_vars.verbose:
        print("--> Downloading {} to {}".format(nc_filename,out_fname))

    # Test that the model repo url works
    try:
        response = get(full_url)
        response.raise_for_status()
    except HTTPError as http_err:
        print('*** HTTP error occurred while trying to access model repo:')
        print(f'{http_err}')
        sys.exit()
    except Exception as err:
        print('An error occurred while trying to access model repo:')
        print(f'{err}')
        sys.exit()
    else:
        # Download the model file at this url
        with open(out_fname,'wb') as dwnld_file:
            dwnld_file.write(response.content)

        if config_vars.verbose:
            print("* Download complete *")

def get_relevant_model_data(config_vars, sat_year):
    # Current catch for resolution
    if config_vars.model_temporal_res != '1 monthly':
        print("--> No access to anything but one monthly data currently.")
        print("--> Quiting program.")
        exit()

    if config_vars.verbose:
        print("--> Finding matching model file")

    filename_format = "mon_{}_{}.nc"
    # Check if the relevant o3, temperature and surface pressure files exists
    for variable in ['o3','t','ps']:
        nc_filename = filename_format.format(variable, sat_year)
        model_filepath = os.path.join(config_vars.model_download_path,
                                        nc_filename)

        if not os.path.isfile(model_filepath):
            download_model_montly_data(config_vars,model_filepath)

    model_data = get_model_data(config_vars,filename_format, sat_year)
    return model_data

def sample_model(config_vars, satellite_info):

    model_year = config_vars.start_date.year
    model_data = get_relevant_model_data(config_vars,
                                                model_year)

    if config_vars.verbose:
        print("--> Sampling model at satellite observations coordinates.")

    # Loop through the satellite files
    for sat_file in satellite_info.file_list:
        sat_data = get_satellite_data.extract_data(config_vars,sat_file)

        satellite_df = pd.DataFrame({'sat_lat':sat_data.latitude,
                                    'sat_lon':sat_data.longitude},
                                    index = sat_data.time)
        satellite_df.index.name = 'Date_Time'
        satellite_df.reset_index(inplace = True)

        # Create empty array to fill with sampled model
        # Currently only set up for RAL files
        # TODO: make this more flexible
        # Have the same number of levels as the model data
        model_o3_profile = np.zeros(sat_data.o3.shape)

        for sample_index, sample_row in satellite_df.iterrows():

            # Find the lat/lon index of the model that matches best with the satellite
            lat_index = utils.get_coord_value(sample_row.sat_lat, model_data.latitude)
            lon_index = utils.get_coord_value(sample_row.sat_lon, model_data.longitude)

            if sample_row.Date_Time.year != model_year:
                model_data = get_relevant_model_data(config_vars,model_year)
                model_year = sample_row.Date_Time.year


            sample_month_index = sample_row.Date_Time.month - 1
            # Get the full column at given coordinates
            sampled_o3 = model_data.model_o3[sample_month_index,:,
                                            lat_index,
                                            lon_index]
            sampled_t = model_data.model_t[sample_month_index,:,
                                            lat_index,
                                            lon_index]

            model_molec_cm2 = utils.model_vmr_to_molec_cm2(config_vars,sampled_t,
                                                            model_data.levels,
                                                            sampled_o3,molec_weight = 48)

            if type(model_molec_cm2) == np.ma.core.MaskedArray:
                # Extract data from masked array
                model_molec_cm2 = model_molec_cm2.filled()
            if type(sat_data.levels) == np.ma.core.MaskedArray:
                sat_data.levels = sat_data.levels.filled()
            if type(sat_data.ak_levels) == np.ma.core.MaskedArray:
                sat_data.ak_levels = sat_data.ak_levels.filled()


            interped_model_o3 = utils.regrid_column(model_molec_cm2,
                                            model_data.levels,
                                            sat_data.ak_levels)

            model_w_aks = apply_aks_to_model(config_vars,
                                            interped_model_o3,
                                            sat_data, sample_index)

            model_w_prior = model_w_aks + sat_data.prior[sample_index]


            model_o3_profile[sample_index,:] = model_w_prior

        sat_data.model_o3 = model_o3_profile

        # Send sampled model to netcdf file
        output_to_file(config_vars, sat_data)

    return sat_data


def apply_aks_to_model(config_vars,model_column,sat_data, sample_index):
    """
        Apply the averaging kernals from the satellite files to the sampled model
        Averaging kernels applied by the dot product of model and aks from sat data

        *** not needed *** => x' = xa + A(xcomp - xa)
        Where xa is the prior state, A is the averaging kernel, xcomp is the model
    """
    # sat_prior = sat_data.prior[sample_index]
    sat_ak = sat_data.aks[sample_index]

    # sat_model_diff = model_column - sat_prior
    # applied_aks = sat_prior  + np.dot(sat_model_diff, sat_ak)
    applied_aks = np.dot(model_column, sat_ak)

    return applied_aks

if __name__ == "__main__":
    # Read in the satellite data and find what files are needed
    satellite_info = get_satellite_data.meta_data(config_vars)
    utils.satellite_date_check(config_vars, satellite_info)

    # Sample the model at the satellite path coordinates
    sat_data = sample_model(config_vars, satellite_info)

    pdb.set_trace()



################################################################################
### END OF PROGRAM
################################################################################
