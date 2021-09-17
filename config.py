################################################################################
#################### Configuration file for TOAR sampling code #################

output_dir = ''
output_filename_prefix = ''

satellite_file_dir = '/Volumes/ScotSat_NO2/RAL_O3/'
multiple_satellite_files = True

# Set start and end date of satellite records as list of integers: YYYY,M,D
start_date = 2019,8,1
end_date = 2020,10,31

satellite_product = {
    'RAL_O3': True,
    'OMI-Aura': False,

}

# Satellite variable names - replace as needed
# Latitude variable name in satellite files
latitude_var_name = 'lat'
# Longitude variable name in satellite files
longitude_var_name = 'lon'
# Level or pressure variable name in satellite files
level_var_name = 'model_levs'
# Time variable name in satellite files
time_var_name  = 'time'
# Ozone variable name in satellite file
o3_var_name = 'o3_ap_sub_col_model'
# Averaging kernal variable name in satellite files
ak_var_name = 'ak_rsc_tsc'
# Prior variable name in satellite files
prior_var_name = 'o3_ap_sub_col_model'

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
        self.lat_var_name = kwargs.get('lat_var_name','')
        self.lon_var_name = kwargs.get('lon_var_name','')
        self.lev_var_name = kwargs.get('lev_var_name','')
        self.time_var_name = kwargs.get('time_var_name','')
        self.ak_var_name = kwargs.get('ak_var_name','')
        self.prior_var_name = kwargs.get('prior_var_name','')
        self.o3_var_name = kwargs.get('o3_var_name','')
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
                                lat_var_name = latitude_var_name,
                                lon_var_name = longitude_var_name,
                                lev_var_name = level_var_name,
                                time_var_name = time_var_name,
                                o3_var_name = o3_var_name,
                                ak_var_name = ak_var_name,
                                prior_var_name = prior_var_name,
                                model_repo_url = model_repo_url,
                                keep_model_downloads = keep_model_downloads,
                                verbose = verbose,
                                model_temporal_res = model_temporal_res)
