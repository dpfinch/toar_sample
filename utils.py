################################################################################
################################################################################
import sys
import os
from datetime import datetime as dt
from datetime import timedelta
from dateutil.relativedelta import relativedelta
import numpy as np
import pandas as pd
import pdb
from scipy.interpolate import interp1d

def days_since_to_dt(time_arr, since_datetime):
    '''
        Converts a "days since XXXX time to datetime"
    '''
    time_unit, unit_base_time = convert_time_since_units(since_datetime)

    # Need to ensure the time arr is in the same float type as the rest of the code
    time_float = [float(h) for h in time_arr[:]]
    if time_unit == 'seconds':
        dt_time_arr = [unit_base_time + relativedelta(seconds =+ i) for i in time_float]
    elif time_unit == 'hours':
        dt_time_arr = [unit_base_time + relativedelta(hours =+ i) for i in time_float]
    elif time_unit == 'months':
        dt_time_arr = [unit_base_time + relativedelta(months =+ i) for i in time_float]
    elif time_unit == 'weeks':
        dt_time_arr = [unit_base_time + relativedelta(weeks =+ i) for i in time_float]
    else:
        return None

    return dt_time_arr


def convert_time_since_units(since_datetime):
    '''
        Find what the approriate units are to parse
        Assumes format of "time since datetime"
    '''
    if type(since_datetime) == str:
        time_str_array = since_datetime.split()
        if time_str_array[0].lower() in ['seconds','hours','days','months']:
            time_unit = time_str_array[0].lower()
            datetime_str = '{} {}'.format(time_str_array[2],time_str_array[3])
            # Use pandas datetime function as it is more forgiving
            unit_base_time = pd.to_datetime(datetime_str)
            return time_unit, unit_base_time

        else:
            print('Cannot parse units {}'.format(since_datetime))
            return None, None

    else:
        print('Time units not string type')
        return None, None

def get_coord_value(coord, model_coord_range):
    """
        find coord index for a given latitude or longitude array
    """

    if coord > 360 or coord < -360:
        print("Satellite coordinate outside range -360 - 360 degrees")
        return 0

    value = min(range(len(model_coord_range)),
                key=lambda i: abs(model_coord_range[i] - coord))

    return value

def get_model_level_index(satellite_levs, model_levs):
    """
        As above but for a pressure level array
    """

    model_lev_index = []

    for sat_lev in satellite_levs:
        model_lev_index.append(min(range(len(model_levs)),
                key=lambda i: abs(model_levs[i] - sat_lev)))

    return model_lev_index

def calculate_box_height(t,p0,p, R_gas):
    """
        Calculate height between two pressure levels for a given temperature
    """

    g = 9.80665 # gravity

    z = -((R_gas * t)/g) * np.log(p/p0)

    return z

def calculate_box_air(t_array,p_array):
    """
        Takes temperate and pressure array and iterates through calulating box air amount
        Calculate amount of air in box based on assumptions
        These assumptions should be fine as we're comparing between satellite
        products so will remain consistent across all satellites
    """
    if type(t_array) in [float, int]:
        t_array = [t_array]
    if type(p_array) in [float, int]:
        p_array = [p_array]

    molec_g_air = 28.97 # molecular weight of air
    R_gas = 287.058 # Gas constant
    avo_num = 6.0021409e23

    box_air_amounts = []
    for box in range(len(p_array)-1):
        z = calculate_box_height(t_array[box], p_array[box], p_array[box + 1], R_gas)

        # Calculate air density in kg/m3
        # Need to convert pres (hPa) to Pa
        air_density_kg = (p_array[box] * 100) / (R_gas * t_array[box])
        air_density_g = air_density_kg * 1000
        # Convert to molec/m3
        molec_per_m3 = (air_density_g / molec_g_air) * avo_num
        # Convert to molec/cm3
        molec_cm3 = molec_per_m3 * 1e-6

        # get molec/cm2 for the whole cell
        air_amount = np.multiply(molec_cm3, z * 100)

        box_air_amounts.append(air_amount)

    return box_air_amounts

def g_m2_from_vmr(altitude_array, surface_temp, press_array, o3_vmr):
    """
        Calculates g/m2 from volume mixing ratio assuming dry air given a give
        altitude array (assumed to be in meters) and surface temp
    """
    molec_g_air = 28.97 # molecular weight of air
    R_gas = 287.058 # Gas constant
    avo_num = 6.0221409e23
    dry_air_lapse_rate = 9.8/1000 # degrees per m

    temp_arr = np.zeros(altitude_array.shape) # Create an empty array to input temperature
    o3_g_cm2 = np.zeros(altitude_array.shape) # Create empty array for the output of ozone

    temp_arr[:,0] = surface_temp # Start array with surface temperature
    for ob in range(temp_arr.shape[0]):
        for lev in range(temp_arr.shape[1]-1):
            start_temp = temp_arr[ob,lev]
            alt_diff = altitude_array[ob,lev + 1] - altitude_array[on,lev]
            temp_drop = alt_diff * dry_air_lapse_rate
            temp_arr[ob,lev + 1] = temp_arr[ob,lev] - temp_drop

        for lev in range(press_array.shape[0]):
            # Calculate air density in kg/m3
            # Need to convert pres (hPa) to Pa
            air_density_kg = (press_array[ob,lev] * 100) / (R_gas * temp_arr[lev])
            air_density_g = air_density_kg * 1000
            # Convert to molec/m3
            molec_per_m3 = (air_density_g / molec_g_air) * avo_num
            # Convert to molec/cm3
            molec_cm3 = molec_per_m3 * 1e-6
            # get molec/cm2 for the whole cell
            air_amount = np.multiply(molec_cm3, z * 100)
            o3_g_cm2[ob,lev] = air_amount * o3_vmr[ob,lev]

    return o3_g_cm2

def model_vmr_to_gm2(config_vars, sampled_t, model_levels, sampled_o3,
                    molec_weight = 48):
    """
        info here...
    """
    avo_num = 6.0221409e23

    # Calculate molecules of air per cm2 in each level of model
    model_cell_air = calculate_box_air(sampled_t, model_levels)
    # Convert v/v to molecules per cm2
    model_molec_o3_percm2 = np.multiply(model_cell_air, sampled_o3)
    model_molec_o3_perm2 = model_molec_o3_percm2 * 1e4 # Convert to m2

    model_g_m2 = (model_molec_o3_perm2 / avo_num) * molec_weight

    return model_g_m2

def regrid_column(o3_col,model_levels,sat_levels):

    """
    Determines the amount of substance within different layers of the
    atmosphere, using interpolation where needed
    """

    cumulative_col = log_interp_sum_to_level(model_levels,o3_col,sat_levels)
    output_col = np.zeros(len(cumulative_col))
    for z in range(len(sat_levels)):
        if z == 0:
            output_col[z] = cumulative_col[z]
        else:
            output_col[z] = cumulative_col[z] - cumulative_col[z-1]

    return output_col

def log_interp_sum_to_level(pres_col,amount_col,out_levels):
    """
        linearly interpolates cumulative of amount_col at pres_points based
        on the logarithm of pres_col
    """

    out_levels_log = np.log(out_levels)

    #pres_col is box bottom pressures, so append in top of atmosphere pressure.
##    pres_col = np.append(pres_col,toa_pres)
    #logarithm
    pres_col_log = np.log(pres_col)

    #get cumulative version of amount
    amount_col_cum = np.nancumsum(amount_col)
    amount_col_cum = np.append(0,amount_col_cum)
##    pdb.set_trace()
    #move any points below lowest and above higest to edge of range
    if type(out_levels_log) != float:
        for i in range(len(out_levels_log)):
            if out_levels_log[i] < min(pres_col_log):
                out_levels_log[i] = min(pres_col_log)
            elif out_levels_log[i] > max(pres_col_log):
                out_levels_log[i] = max(pres_col_log)

    # create interpolation function
    f = interp1d(pres_col_log,amount_col_cum)

    interp_column = f(out_levels_log)

    return interp_column

def DU_to_g_m2(spc_column, molc_weight = 48):
    """
        Convert dobson units to g/m2, default is for molecular weight of ozone (48)
    """
    avo_num = 6.0221409e23
    conv_rate = 2.6867e20 # 1 DU == 2.687e20 molecules per square meter
    molec_per_m2 = spc_column * conv_rate
    mol_per_m2 = molec_per_m2 / avo_num
    g_per_m2 = mol_per_m2 * molc_weight

    return g_per_m2

def molec_to_DU(column):
    """
        Convert column of molecules per cm2 to dobson units. Needed for RAL analysis
    """
    standard_temp = 273.15 # Zero celcius in kelvin
    standard_pres = 101325.0 # Pressure (Pa)
    boltz = 1.3806581e-23 # Boltzmann constant (JK-1)

    dobson_units = 1e9 * column * boltz * (standard_temp/standard_pres)

    return dobson_units

def molec_column_calc(vmr_column):
    """
        Compute sub column amount from VMR on pressure levels. Needed for RAL analysis
    """
    avo_num = 6.0221409e23 # Avogadros number
    g = 9.80665 # Gravity (m/s2)
    molec_g_air = 28.97 # molecular weight of air

    molec_col = vmr_column * avo_num/g/molec_g_air/1e2 # Taken from RAL IDL code

    return molec_col

def remove_model_file(config_vars, out_fname):
    '''
        Delete the downloaded model files once they have been sampled to save space
    '''
    # TODO: complete this function
    if config_vars.verbose:
        print("Deleting {}".format(out_fname))

    if os.path.exists(out_fname):
        os.remove(out_fname)
    else:
        print("Cannot find file {} to delete.".format(out_fname))


def satellite_date_check(config_vars, satellite_info):
    '''
        Check the satellite files date and change the config_vars to match
    '''
    if satellite_info.start_date > config_vars.start_date:
        print("--> Satellite records start from {}. \n Changing start date to match".format(
            satellite_info.start_date.date()))
        config_vars.start_date = satellite_info.start_date

    if satellite_info.end_date < config_vars.end_date:
        print("--> Satellite records end at {}. \n Changing end date to match".format(
            satellite_info.end_date.date()))
        config_vars.end_date = satellite_info.end_date
################################################################################
### END OF PROGRAM
################################################################################
