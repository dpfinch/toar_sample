################################################################################
#################### Configuration file for TOAR sampling code #################

output_dir = '/'

# Set which satellite product you are using. Make sure all others are set to False
satellite_product = {
    'OMI-RAL': True,
    'OMI-Aura': False,
    'IASI': False,
    'OMPS-NPP': False,
}

satellite_file_dir = '/geos/d21/dfinch/'

satellite_file_dir = '/'

# Set the file suffix for the satellite data. (e.g. nc for netcdf)
satellite_file_suffix = 'h5'

# Set start and end date of satellite records as list of integers: YYYY,M,D
start_date = 2019,8,1
end_date = 2019,8,31

#### SATELLITE VARIABLES ##### - replace as needed ####
# Latitude variable name in satellite files
latitude_var_name = 'lat'
# Longitude variable name in satellite files
longitude_var_name = 'lon'
# Level or pressure variable name in satellite files - can be either 1D or 2D
level_var_name = 'levs'
# Averaging kernels may be on a different level grid - if so put the name here
# If they are on the same levels as ozone then either write the same name or None
ak_level_var_name = None
# Time variable name in satellite files
time_var_name  = ['GeolocationData/Year','GeolocationData/DayOfYear','GeolocationData/SecondsInDay']
# Ozone variable name in satellite file
o3_var_name = 'ScienceData/O3MixingRatio'
# Averaging kernal variable name in satellite files
ak_var_name = 'ScienceData/AveragingKernels'
# aPriori contribution to esimated sub column variable name in satellite files
prior_var_name = 'ScienceData/ProfileO3FirstGuess'
# Apply quality control to the satellite retrieval. This currently only work for RAL data
apply_quality_control = True # Currently doesn't do anything. Will be fixed soon.

##### MODEL VARIABLES #####

# Choose which model to use (CAMS advised)

CAMS_model = True
NASA_TES_model = False

keep_model_downloads = False

verbose = True

# Set the temporal resoltion of the model (n monthly, n daily or n hourly)
model_temporal_res = '1 monthly'

################################################################################
#########################   END OF VARIABLE ASSIGNMENT  ########################
