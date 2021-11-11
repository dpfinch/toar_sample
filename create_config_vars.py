
################################################################################
######################### CREATING VARIABLE CLASS BELOW ########################
###### No need to edit below unless making major changes to code ###############
from config import *
import os
import sys
from datetime import datetime
import numpy as np

# Create objection containing all the variables set above
class create_config_vars:
    def __init__(self,**kwargs):

        self.output_dir = kwargs.get('output_dir','')
        self.satellite_product = kwargs.get('satellite_product','')
        self.satellite_file_dir = kwargs.get('satellite_file_dir','')
        self.satellite_file_suffix = kwargs.get('satellite_file_suffix','nc')
        self.apply_quality_control = kwargs.get('apply_quality_control',False)
        self.lat_var_name = kwargs.get('lat_var_name','')
        self.lon_var_name = kwargs.get('lon_var_name','')
        self.lev_var_name = kwargs.get('lev_var_name','')
        self.ak_lev_var_name = kwargs.get('ak_lev_var_name',self.lev_var_name)
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
    print("** Model temporal resoltion ({}) not a valid input./n Format should be N hourly, N daily or N monthly".format(
        model_temporal_res
    ))
    sys.exit()
if model_temporal_res.split()[1].lower() not in ['hourly','daily','monthly']:
    print("** Model temporal resoltion ({}) not a valid input./n Format should be N hourly, N daily or N monthly".format(
        model_temporal_res
    ))
    sys.exit()

if np.sum(list(satellite_product.values())) > 1:
    print("** More than one satellite product has been selected in config.py. Only choose one.")
    sys.exit()

elif np.sum(list(satellite_product.values())) < 1:
    print("** No satellite product has been selected. Choose one.")
    sys.exit()

else:
    satellite_product_name = [k for k,v in satellite_product.items() if v][0]

print("### SAMPLING MODEL FOR {} ###".format(satellite_product_name))

if ak_level_var_name == None:
    ak_level_var_name = level_var_name

config_vars = create_config_vars(output_dir = output_dir,
                                satellite_product = satellite_product_name,
                                satellite_file_dir =  satellite_file_dir,
                                start_date = start_date,
                                end_date = end_date,
                                satellite_file_suffix = satellite_file_suffix,
                                apply_quality_control = apply_quality_control,
                                lat_var_name = latitude_var_name,
                                lon_var_name = longitude_var_name,
                                lev_var_name = level_var_name,
                                ak_lev_var_name = ak_level_var_name,
                                time_var_name = time_var_name,
                                o3_var_name = o3_var_name,
                                ak_var_name = ak_var_name,
                                prior_var_name = prior_var_name,
                                model_repo_url = model_repo_url,
                                keep_model_downloads = keep_model_downloads,
                                verbose = verbose,
                                model_temporal_res = model_temporal_res)
