
################################################################################
################################################################################
import cdsapi
import xarray as xr

def download_model_data():
    c = cdsapi.Client()

    outfile = 'download.grib'

    pressure_levs = ['1', '2', '3',
                '5', '7', '10',
                '20', '30', '50',
                '70', '100', '150',
                '200', '250', '300',
                '400', '500', '600',
                '700', '800', '850',
                '900', '925', '950',
                '1000']

    date_var = '2003-01-01/2003-01-02'
   
    time_list = ['00:00', '03:00', '06:00',
                '09:00', '12:00', '15:00',
                '18:00', '21:00']

    species = 'ozone'

    north_lat = 60
    south_lat = 50
    west_lon = -10
    east_lon = 5

    area_extent = [north_lat, west_lon, south_lat, east_lon]

    c.retrieve(
        'cams-global-reanalysis-eac4',
        {
            'date': date_var,
            'format': 'grib',
            'variable': species,
            'pressure_level': pressure_levs,
            'time': time_list,
            'area': area_extent,
        },
        outfile)

    return outfile

def open_file(outfile_name):
    data = xr.open_dataset(outfile_name, engine='cfgrib')
    # Could do multi dataset with xr.open_mfdataset(files)
    pass

def delete_file(outfile_name):
    pass

if __name__ == '__main__':
    download_model_data()