################################################################################
################################################################################
import utils
import netCDF4 as nc
import sys
import os
import pdb
import numpy as np
#  This prevents a userwarning about converting masked arrays to nans
import warnings
warnings.simplefilter("ignore", UserWarning)
warnings.simplefilter("ignore", RuntimeWarning)


def output_to_file(config_vars, sat_data):
    """
        Output the sampled model to netcdf file.
    """

    first_time_str = sat_data.time[0].strftime('%Y%m%d%H%M')
    last_time_str = sat_data.time[-1].strftime('%Y%m%d%H%M')

    outfilename = 'Sampled_Model_{}_{}.nc'.format(first_time_str,last_time_str)

    outfile_full_path = (os.path.join(config_vars.output_dir, outfilename))

    if config_vars.verbose:
        print("--> Creating file: {}".format(outfilename))
        print("    In directory: {}".format(config_vars.output_dir))

    if not os.path.isdir(config_vars.output_dir):
        print("Cannot find: {}. Creating directory.".format(config_vars.output_dir))
        os.makedirs(config_vars.output_dir)

    outfile = nc.Dataset(outfile_full_path,'w')

    lat_array = sat_data.latitude
    lon_array = sat_data.longitude
    lev_array = sat_data.sat_level_mids
    time_array = utils.dt_to_seconds_since(sat_data.time)

    out_lat = outfile.createDimension('Latitude',len(lat_array))
    out_lon = outfile.createDimension('Longitude',len(lon_array))
    out_time = outfile.createDimension('Time',len(time_array))
    out_profiles = outfile.createDimension('n_profiles',len(time_array))
    out_nlevels = outfile.createDimension('n_levels',lev_array.shape[1])
    # out_levels = outfile.createDimension('Levels',len(lev_array))
    lats = outfile.createVariable('Latitude','f8',('Latitude',))
    lons = outfile.createVariable('Longitude','f8',('Longitude',))
    levs = outfile.createVariable('Levels','f8',('n_profiles','n_levels',))
    times = outfile.createVariable('Time','f8',('Time',))

    lats.units = 'Degrees North'
    lats[:] = lat_array
    lons.units = 'Degrees East'
    lons[:] = lon_array
    times.units = 'Seconds since 1990-01-01 00:00:00'
    times[:] = time_array
    levs.units = 'Pressure (hPa)'
    levs[:] = lev_array

    # Find if there are duplicated in the time array.
    # If there are then it cannot be used for a dimension
    unique_times = set(time_array)
    if len(unique_times) != time_array:
        if config_vars.verbose:
            print("--> Time array has duplicates. Cannot be used as ozone dimension.")
            print("--> Setting dimension to 'n_profiles' instead.")
        n_profiles_dim = 'n_profiles'
    else:
        n_profiles_dim = 'Time'

    orig_o3 = outfile.createVariable('Satellite_O3',
                                     'f8',('n_profiles','n_levels',))
    orig_o3.units = 'molecules/cm2'
    orig_o3[:] = sat_data.o3

    orig_o3 = outfile.createVariable('Sampled_model_O3',
                                     'f8',('n_profiles','n_levels',))
    orig_o3.units = sat_data.o3_units
    orig_o3[:] = sat_data.model_o3

    outfile.close()
