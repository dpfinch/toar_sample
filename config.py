################################################################################
#################### Configuration file for TOAR sampling code #################

output_dir = '/'

# Set which satellite product you are using. Make sure all others are set to False
satellite_product = {
    'OMI-RAL': False,
    'OMI-Aura': False,
    'IASI': False,
    'OMPS-NPP': True,
}

satellite_file_dir = '/'
# Set the file suffix for the satellite data. (e.g. nc for netcdf)
satellite_file_suffix = 'h5'

# Set start and end date of satellite records as list of integers: YYYY,M,D
start_date = 2019,8,1
end_date = 2020,10,31

#### SATELLITE VARIABLES ##### - replace as needed ####
# Latitude variable name in satellite files
latitude_var_name = 'GeolocationData/Latitude'
# Longitude variable name in satellite files
longitude_var_name = 'GeolocationData/Longitude'
# Level or pressure variable name in satellite files - can be either 1D or 2D
level_var_name = 'DimPressureLevel'
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
model_repo_url = 'https://tes.jpl.nasa.gov/raw-data/tcr-2_files/monthly-mean_emissions/'
# TODO: Set up function to delete files
keep_model_downloads = True

verbose = True

# Set the temporal resoltion of the model (n monthly, n daily or n hourly)
model_temporal_res = '1 monthly'

################################################################################
#########################   END OF VARIABLE ASSIGNMENT  ########################
