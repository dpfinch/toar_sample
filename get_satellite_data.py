################################################################################
################################################################################
from config import config_vars
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
        self.multiple_files = config_vars.multiple_satellite_files
        self.file_list = self.get_file_list()
        self.num_files = len(self.file_list)
        self.start_date, self.end_date = self.get_satellite_date_range(config_vars)

    def get_file_list(self):
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
            print("--> Found {} satellite files".format(len(files)))
        return files

    def get_satellite_date_range(self, config_vars):
        """
            Find the start and end dates of satellite files. Assumes that ordering
            the satellite files by file name will be chronological.

        """
        # Open the first satellite file to get start date
        nc_dataset = nc.Dataset(self.file_list[0],'r')
        pass_time = nc_dataset.variables[config_vars.time_var_name]
        pass_time_unit = pass_time.units
        pass_dt = utils.days_since_to_dt(pass_time, pass_time_unit)
        if type(pass_dt) == list:
            start_date = pass_dt[0]
        else:
            # TODO: Probably need to some error catching here
            start_date = pass_dt
        nc_dataset.close()
        # Open last satellite file to get end date
        nc_dataset = nc.Dataset(self.file_list[-1],'r')
        pass_time = nc_dataset.variables[config_vars.time_var_name]
        pass_time_unit = pass_time.units
        pass_dt = utils.days_since_to_dt(pass_time, pass_time_unit)
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
        alt_var = config_vars.altitude_var_name
        pres_var = config_vars.pressure_var_name
        time_var = config_vars.time_var_name

        self.latitude = nc_dataset.variables[lat_var][:]
        self.longitude = nc_dataset.variables[lon_var][:]
        self.levels = nc_dataset.variables[lev_var][:]
        self.ak_levels = nc_dataset.variables[ak_lev_var][:]
        self.aks = nc_dataset.variables[ak_var][:]
        self.o3 = nc_dataset.variables[o3_var][:]

        #  Apply level check to make sure dimensions match
        self.level_check()

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
        if config_vars.include_altitude:
            self.altitude = nc_dataset.variables[alt_var][:]
        else:
            self.altitude = None
        if config_vars.include_pressure:
            self.pressure = nc_dataset.variables[pres_var][:]
        else:
            self.pressure = None

        pass_time = nc_dataset.variables[time_var]
        pass_time_unit = pass_time.units

        self.time = utils.days_since_to_dt(pass_time, pass_time_unit)

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
            elif self.prior_units == 'ppb':
                self.prior = self.prior * 1e9
            else:
                print("** Unable to converts o3 prior units of {}. Qutting.".format(self.prior_units))
                sys.exit()
            self.prior_units = 'molec cm-2'

        nc_dataset.close()
        #  Check if the averaging kernals have some matching dimension to the ozone levels.
        # Otherwise we cannot map the kernels to the profile.
        if len(self.levels) not in self.aks.shape:
            print("** Number of levels in {} does not match any dimension in the averaging kernels.".format(lev_var))
            print("** Number of levels = {}".format(len(self.ak_levels)))
            print("** Averaging kernels dimensions = {}".format(self.aks.shape))
            print("** They need to match. Quitting.")
            sys.exit()

    def level_check(self):
        """
            Check the levels match the profiles of ozone and the averaging kernels.
            Its assumed that if the levels are one greater than the profiles that they
            represent the level edges and therefore we can find the midpoint of the level
            and use that.
        """
        if self.levels.ndim == 1:
            if len(self.levels) not in self.o3.shape:
                print("--> Levels provided do not match O3 profile.")
            if len(self.levels) - 1 in self.o3.shape:
                print("--> Levels are one greater than O3 profile.")
                print("    Assuming they are level edges and will interpolate between them.")
                print("    Calculating mid point between levels. Assuming linear interpolation.")

                current_levels = self.levels
                num_levs = len(current_levels)
                level_diffs = np.diff(current_levels)

                if all([current_levels[x] > current_levels[x+1] for x in range(num_levs -1)]):
                    levels_mids = current_levels[1:] - (level_diffs/2)
                elif all([current_levels[x] < current_levels[x+1] for x in range(num_levs -1)]):
                    levels_mids = current_levels[:-1] - (level_diffs/2)
                else:
                    print("** Levels are not in order. Quitting")
                    sys.exit()
                self.levels = levels_mids
            else:
                print("** Quiting.")
                sys.exit()

        if self.ak_levels.ndim == 1:
            if len(self.ak_levels) not in self.aks.shape:
                print("--> Averaging kernel levels provided do not match AK profile.")
            if len(self.ak_levels) - 1 in self.aks.shape:
                print("--> AK levels are one greater than AK profile.")
                print("    Assuming they are level edges and will interpolate between them.")
                print("    Calculating mid point between AK levels. Assuming linear interpolation.")

                current_levels = self.ak_levels
                num_levs = len(current_levels)
                level_diffs = np.diff(current_levels)

                if all([current_levels[x] > current_levels[x+1] for x in range(num_levs -1)]):
                    levels_mids = current_levels[1:] - (level_diffs/2)
                elif all([current_levels[x] < current_levels[x+1] for x in range(num_levs -1)]):
                    levels_mids = current_levels[:-1] - (level_diffs/2)
                else:
                    print("** Averaging kernel levels are not in order. Quitting")
                    sys.exit()
                self.ak_levels = levels_mids

            else:
                print("** Quiting.")
                sys.exit()
        else:
            print("** Level array is not one dimensional and this program cannot handle that yet.")
            print("** Raise this issue on the GitHub repo and it'll get fixed.")

################################################################################
### END OF PROGRAM
################################################################################
