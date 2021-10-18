################################################################################
#################### Configuration file for TOAR sampling code #################

output_dir = ''
output_filename_prefix = ''

satellite_file_dir = '/Users/dfinch/Desktop/RAL_O3/'
multiple_satellite_files = True

# Set start and end date of satellite records as list of integers: YYYY,M,D
start_date = 2019,8,1
end_date = 2020,10,31

satellite_product = {
    'RAL_O3': True,
    'OMI-Aura': False,

}

#### SATELLITE VARIABLES ##### - replace as needed ####
# Set the file suffix for the satellite data. (e.g. nc for netcdf)
satellite_file_suffix = 'nc'
# Latitude variable name in satellite files
latitude_var_name = 'lat'
# Longitude variable name in satellite files
longitude_var_name = 'lon'
# Level or pressure variable name in satellite files
level_var_name = 'model_levs'
# Time variable name in satellite files
time_var_name  = 'time'
# Ozone variable name in satellite file
o3_var_name = 'o3_sub_col'
# o3_var_name = 'x' # AIRS-Aqua names
# Averaging kernal variable name in satellite files
ak_var_name = 'ak_rsc_tsc'
# ak_var_name = 'observation_ops/averaging_kernel'
# aPriori contribution to esimated sub column variable name in satellite files
prior_var_name = 'imak_apr_sc'
# prior_var_name = 'observation_ops/xa'
# Apply quality control to the satellite retrieval. This currently only work for RAL data
apply_quality_control = True

## The satellite ozone data will need to be converted to g/m2.
## To do this its likely that pressure & altitude are needed - enter variables below
# Does the satellite file also include pressure of each satellite observation layer?
include_pressure = False
pressure_var_name = 'pressure'
# Does the satellite file also include altitude of each satellite observation layer?
# This is assumed to be in meters
include_altitude = False #
altitude_var_name = '/observation_ops/altitude'

##### MODEL VARIABLES ##### https://tes.jpl.nasa.gov/tes/chemical-reanalysis/products/monthly-mean
model_repo_url = 'https://tes.jpl.nasa.gov/raw-data/tcr-2_files/monthly-mean_emissions/'
# TODO: Set up function to delete files
keep_model_downloads = False

verbose = True

# Set the temporal resoltion of the model (n monthly, n daily or n hourly)
model_temporal_res = '1 monthly'


################################################################################
######################### CREATING VARIABLE CLASS BELOW ########################
###### No need to edit below unless making major changes to code ###############

import os
from datetime import datetime

# Create objection containing all the variables set above
class create_config_vars:
    def __init__(self,**kwargs):

        self.output_dir = kwargs.get('output_dir','')
        self.output_filename_prefix = kwargs.get('output_filename_prefix','')
        self.satellite_file_dir = kwargs.get('satellite_file_dir','')
        self.multiple_satellite_files = kwargs.get('multiple_satellite_files',True)
        self.satellite_file_suffix = kwargs.get('satellite_file_suffix','nc')
        self.apply_quality_control = kwargs.get('apply_quality_control',False)
        self.lat_var_name = kwargs.get('lat_var_name','')
        self.lon_var_name = kwargs.get('lon_var_name','')
        self.lev_var_name = kwargs.get('lev_var_name','')
        self.time_var_name = kwargs.get('time_var_name','')
        self.ak_var_name = kwargs.get('ak_var_name','')
        self.prior_var_name = kwargs.get('prior_var_name','')
        self.o3_var_name = kwargs.get('o3_var_name','')
        self.include_altitude = kwargs.get('include_altitude',False)
        self.altitude_var_name = kwargs.get('altitude_var_name','')
        self.include_pressure = kwargs.get('included_pressure', False)
        self.pressure_var_name = kwargs.get('pressure_var_name','')
        self.model_repo_url = kwargs.get('model_repo_url','')
        self.keep_model_downloads = kwargs.get('keep_model_downloads',True)
        self.verbose = kwargs.get('verbose',True)
        self.model_temporal_res = kwargs.get('model_temporal_res', '')
        self.start_date = datetime(*kwargs.get('start_date',(2000,1,1)))
        self.end_date = datetime(*kwargs.get('end_date',(2020,1,1)))

        model_download_path = os.path.join(self.output_dir, 'model_download')
        self.model_download_path = model_download_path

# Make sure the model url ends with a slash otherwise the URL join doesn't work
if model_repo_url[-1] != '/':
    model_repo_url = model_repo_url + '/'
# Make sure the model temporal resolution makes sense
try:
    num_res = int(model_temporal_res.split()[0])
except (NameError, SyntaxError) as e:
    print("Model temporal resoltion ({}) not a valid input./n Format should be N hourly, N daily or N monthly".format(
        model_temporal_res
    ))
    exit()
if model_temporal_res.split()[1].lower() not in ['hourly','daily','monthly']:
    print("Model temporal resoltion ({}) not a valid input./n Format should be N hourly, N daily or N monthly".format(
        model_temporal_res
    ))
    exit()

config_vars = create_config_vars(output_dir = output_dir,
                                output_filename_prefix = output_filename_prefix,
                                satellite_file_dir =  satellite_file_dir,
                                multiple_satellite_files = multiple_satellite_files,
                                start_date = start_date,
                                end_date = end_date,
                                satellite_file_suffix = satellite_file_suffix,
                                apply_quality_control = apply_quality_control,
                                lat_var_name = latitude_var_name,
                                lon_var_name = longitude_var_name,
                                lev_var_name = level_var_name,
                                time_var_name = time_var_name,
                                o3_var_name = o3_var_name,
                                ak_var_name = ak_var_name,
                                prior_var_name = prior_var_name,
                                include_altitude = include_altitude,
                                altitude_var_name = altitude_var_name,
                                include_pressure = include_pressure,
                                pressure_var_name = pressure_var_name,
                                model_repo_url = model_repo_url,
                                keep_model_downloads = keep_model_downloads,
                                verbose = verbose,
                                model_temporal_res = model_temporal_res)
