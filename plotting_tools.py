################################################################################
################################################################################
from create_config_vars import config_vars
import netCDF4 as nc
import sys
import os
from glob import glob
import utils
import pandas as pd
import numpy as np
import pdb

class plot_data():
    def __init__(self, config_vars):
        self.sampled_outdir = config_vars.output_dir
        
    # Setting Troposphere requirements as to HEGIFTOM:
    # From ground to 150 hPa in the tropics (within 15° )
    # From ground to 200 hPa in the subtropics (15°-30°)
    # From ground to 300 hPa in the midlatitudes (30°-60°)
    # From ground to 400 hPa in the polar regions (> 60°)
    def difference_timeseries(self):
        pass


################################################################################
### END OF PROGRAM
################################################################################
