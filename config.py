################################################################################
#################### Configuration file for TOAR sampling code #################

output_dir = '/Volumes/DougsHD/'
output_filename_prefix = 'RAL_'

satellite_file_dir = '/Volumes/DougsHD/RAL_O3/'
multiple_satellite_files = True

# Set start and end date of satellite records as list of integers: YYYY,M,D
start_date = 2019,8,1
end_date = 2020,10,31

# Set which satellite product you are using. Make sure all others are set to False
satellite_product = {
    'RAL_O3': True,
    'OMI-Aura': False,
    'IASI': False,
}

#### SATELLITE VARIABLES ##### - replace as needed ####
# Set the file suffix for the satellite data. (e.g. nc for netcdf)
satellite_file_suffix = 'nc'
# Latitude variable name in satellite files
latitude_var_name = 'lat'
# Longitude variable name in satellite files
longitude_var_name = 'lon'
# Level or pressure variable name in satellite files
level_var_name = 'levs'
# Averaging kernels may be on a different level grid - if so put the name here
# If they are on the same levels as ozone then either write the same name or write None
ak_level_var_name = 'model_levs'
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
#########################   END OF VARIABLE ASSIGNMENT  ########################
