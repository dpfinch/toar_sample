################################################################################
################################################################################
import netCDF4 as nc
import sys
import os
from glob import glob
import utils
import pandas as pd
import numpy as np
import pdb


class meta_data():
    def __init__(self,config_vars):
        self.satellite_file_dir = config_vars.satellite_file_dir
        self.satellite_file_suffix = config_vars.satellite_file_suffix
        # Remove the full stop infront of the suffix if there is one.
        self.satellite_file_suffix = self.satellite_file_suffix.replace('.','')
        self.file_list = self.get_file_list(config_vars)
        self.num_files = len(self.file_list)
        self.start_date, self.end_date = self.get_satellite_date_range(config_vars)

    def get_file_list(self, config_vars):
        '''
            This needs to be updated to be more general
        '''
        files = glob('{}/**/*.{}'.format(self.satellite_file_dir,
                                        self.satellite_file_suffix), recursive = True)
        files.sort()

        if len(files) == 0:
            print("*****************************************************************")
            print("*** No satellite files found. Make sure config.py is correct. ***")
            print("*** Looking for files here:                                   ***")
            print("   ",self.satellite_file_dir)
            print("*****************************************************************")
            sys.exit()

        if config_vars.verbose:
            if len(files) > 1:
                suffix = 's'
            else:
                suffix = ''
            print("--> Found {} satellite file{}".format(len(files),suffix))
        return files

    def get_satellite_date_range(self, config_vars):
        """
            Find the start and end dates of satellite files. Assumes that ordering
            the satellite files by file name will be chronological.

        """
        # Open the first satellite file to get start date
        nc_dataset = nc.Dataset(self.file_list[0],'r')
        pass_time = nc_dataset.variables[config_vars.time_var_name]
        if config_vars.satellite_product != 'IASI':
            pass_time_unit = pass_time.units
            pass_dt = utils.days_since_to_dt(pass_time, pass_time_unit)
        else:
            pass_dt = utils.deal_with_IASI_time(config_vars,nc_dataset, only_start = True)
        if type(pass_dt) == list:
            start_date = pass_dt[0]
        else:
            # TODO: Probably need to some error catching here
            start_date = pass_dt
        nc_dataset.close()
        # Open last satellite file to get end date
        nc_dataset = nc.Dataset(self.file_list[-1],'r')
        pass_time = nc_dataset.variables[config_vars.time_var_name]
        if config_vars.satellite_product != 'IASI':
            pass_time_unit = pass_time.units
            pass_dt = utils.days_since_to_dt(pass_time, pass_time_unit)
        else:
            pass_dt = utils.deal_with_IASI_time(config_vars,nc_dataset, only_start = True)
        if type(pass_dt) == list:
            end_date = pass_dt[-1]
        else:
            # TODO: Probably need to some error catching here
            end_date = pass_dt
        nc_dataset.close()

        return start_date,end_date

class extract_data():
    def __init__(self,config_vars,sat_file):

        # TODO: Change to a 'with' statement
        # TODO: Need to make this configurable
        nc_dataset = nc.Dataset(sat_file,'r')
        # Try to get coord variables either as full names or shortened
        lat_var = config_vars.lat_var_name
        lon_var = config_vars.lon_var_name
        lev_var = config_vars.lev_var_name
        ak_lev_var = config_vars.ak_lev_var_name
        o3_var = config_vars.o3_var_name
        ak_var = config_vars.ak_var_name
        prior_var = config_vars.prior_var_name
        time_var = config_vars.time_var_name

        self.latitude = nc_dataset.variables[lat_var][:]
        self.longitude = nc_dataset.variables[lon_var][:]
        self.levels = nc_dataset.variables[lev_var][:]
        self.ak_levels = nc_dataset.variables[ak_lev_var][:]
        self.aks = nc_dataset.variables[ak_var][:]
        self.o3 = nc_dataset.variables[o3_var][:]

        # Create variable to contain the mid levels of the array
        self.sat_level_mids = np.full(self.o3.shape,-999)
        self.ak_level_mids= np.full(self.aks.shape[:2],-999)
        # Convert time variable into datetime object
        pass_time = nc_dataset.variables[time_var]
        if config_vars.satellite_product != 'IASI':
            pass_time_unit = pass_time.units
            self.time = utils.days_since_to_dt(pass_time, pass_time_unit)
        else:
            self.time = utils.deal_with_IASI_time(config_vars,nc_dataset)

        # Get units for ozone
        try:
            self.o3_units = nc_dataset.variables[o3_var].units
        except AttributeError:
            print("--> Cannot find units for ozone. Assuming v/v.")
            self.o3_units = 'v/v'
        self.prior = nc_dataset.variables[prior_var][:]
        try:
            self.prior_units = nc_dataset.variables[prior_var].units
        except AttributeError:
            print("--> Cannot find units for prior variable. Assumiung same as ozone.")
            self.prior_units = self.o3_units

        # Get units for the prior
        if self.prior_units != self.o3_units:
            print("** Units for prior do not match units for ozone. Quitting.")
            sys.exit()

        if self.o3_units in ['1','',None]:
            if config_vars.verbose:
                print("--> Ozone has units of '{}'. Assuming this is v/v".format(self.o3_units))
            self.o3_units = 'v/v'
        if self.prior_units in ['1','',None]:
            if config_vars.verbose:
                print("--> Ozone prior has units of '{}'. Assuming this is v/v".format(self.prior_units))
            self.prior_units = 'v/v'

        if self.o3_units != 'molec cm-2':
            if config_vars.verbose:
                print("--> Converting ozone from units of {} to molec/cm2".format(self.o3_units))
            if self.o3_units == 'DU':
                self.o3 = utils.DU_to_molec_cm2(self.o3)
            elif self.o3_units == 'ppb':
                self.o3 = self.o3 * 1e9
            elif self.o3_units == 'mol m-2':
                self.o3 = utils.molm2_to_molec_cm2(self.o3)
            else:
                print("** Unable to converts o3 units of {}. Qutting.".format(self.o3_units))
                sys.exit()
            self.o3_units = 'molec cm-2'

        if self.prior_units != 'molec cm-2':
            if config_vars.verbose:
                print("--> Converting prior from units of {} to molec/cm2".format(self.prior_units))
            if self.prior_units == 'DU':
                self.prior = utils.DU_to_molec_cm2(self.prior)
            elif self.prior_units == 'mol m-2':
                self.prior = utils.molm2_to_molec_cm2(self.prior)
            elif self.prior_units == 'ppb':
                self.prior = self.prior * 1e9
            else:
                print("** Unable to converts o3 prior units of {}. Qutting.".format(self.prior_units))
                sys.exit()
            self.prior_units = 'molec cm-2'

        # Get units for levels
        try:
            self.level_units = nc_dataset.variables[lev_var].units
        except AttributeError:
            print("--> Cannot find units for level. Assuming hPa.")
            self.level_units = 'hPa'
        if self.level_units.lower() not in ['hpa','hectopascal']:
            if self.level_units.lower() in ['pa','pascal']:
                self.levels = self.levels / 100
                self.level_units = 'hPa'
            else:
                print("** Unable to convert units of {} for the levels at this time".format(self.level_units))
                print("** Qutting.")
                sys.exit()
        else:
            self.level_units = 'hPa' # Just to makes sure its corrrectly written.

        nc_dataset.close()
        #  Check if the averaging kernals have some matching dimension to the ozone levels.
        # Otherwise we cannot map the kernels to the profile.
        # if len(self.levels) not in self.aks.shape:
        #     print("** Number of levels in {} does not match any dimension in the averaging kernels.".format(lev_var))
        #     print("** Number of levels = {}".format(len(self.ak_levels)))
        #     print("** Averaging kernels dimensions = {}".format(self.aks.shape))
        #     print("** They need to match. Quitting.")
        #     sys.exit()

################################################################################
### END OF PROGRAM
################################################################################
