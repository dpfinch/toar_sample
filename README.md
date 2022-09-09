# TOAR Satellite Sample

### Overview & Purpose

This software has been created as part of the second phase in the **T**ropospheric **O**zone **A**ssement **R**eport (TOAR) project. It is designed to assist in understanding satellite observations of tropospheric ozone and quantify the difference between satellite products. 

Different satellite products are inconsistent when reporting the tropospheric ozone burden in the Earths atmosphere. We need to find out whether these differences are due to the differences in sampling and data processing of each satellite or if there is something other reason for the differences. To do this we use a consistent simulated atmosphere, sampled in the same way as each satellite. As we know the observed atmosphere is consistent we should be able to replecate the differences in the satellites via this method, helping us determine where the differences are occuring from. 

### How it Works

This program reads in satellite files and extracts the latitude, longitude and time of each observation along with any other necessary information, such as averaging kernals, layer profiles and quality flags. Once the program has deteremined where and when the satellite is sampling, it download the relevant chemical reanalysis data from [CAMS](https://confluence.ecmwf.int/display/CKB/CAMS%3A+Reanalysis+data+documentation#heading-Levellistings).

This model data is a 3-hourly from 00-UTC 3D products at a roughyl 80 km resolution available between 2005-2019. High spatial and temporal resolution will be availble in the future from [JPL](https://tes.jpl.nasa.gov/tes/chemical-reanalysis/products/monthly-mean/) and the code will be updated to allow both model datasets to be used.

The program will the iterate through each observation and extract the model profile at the same coordinate and time as the observation. It will then interpolate the model levels to match the satellite levels and apply the satellite averaging kernels if necessary and produce a comparable product to the satellite input.

The sampled model data will then be written to a file in the same format as the satellite input (or as close to as possible) to be analysed at a later date.

### Downloading the Code

To download the code type the following command into your terminal:

`git clone https://github.com/dpfinch/toar_sample.git`

### Requirements

This program runs on Python and has been developed using Python 3.7. It should work with any version 3 Python but has not been tested with any Python 2 versions. There are no plans to support Python 2 development as it is now outdated.

It is recommeneded that you set up a virtual environment to ensure the correct versions of modules are used. Instruction to do this can be found [here](https://docs.python-guide.org/dev/virtualenvs/). 

To make the code into a virtual envirionment type the following into your terminal:

`virtualenv --python=python3 toar_sample`

(You may need to change `python3` to either a specific path where you python is installed or a specific version e.g. 3.7)

All necessary modules can be foind in the requirements.txt file. Run `pip install -r requirements.txt` within your virtual environment to install all required modules.

It is possible to install these one by one and as needed by running the program and errors occur when a module is needed and not been installed. This method is not advised as you may end up with the wrong version of the module and there is no guarantee it will work.

This program uses the Copernicus CDS API to pull the model data. You will need a (free) account to access the data and for the sampling software to work. Information on this can be found [here](https://ads.atmosphere.copernicus.eu/api-how-to#install-the-cds-api-key)

### Running the Program

##### config.py

The first step is to edit the `config.py`  file with the relevant information. This will be fully documents at somepoint the future but below is a starting point.

`output_dir` - Set where the output files will be written. The model data will also be downloaded here.

`satellite_product` - This is a python dictionary where you need to set one (and only one) satellite product to `True` and the rest to `False`. This determines how to processes different styles of satellite file. This is an incomplete list and if there is a satellite product you want to use that is not on the list then raise the issue on the GitHub page.

`satellite_file_dir` - Set the parents directory of the satellite data to be analysed, this program will look in subdirectories so no need to include them (e.g. /YYYY/MM/satellite_data.nc)

`satellite_file_suffix` - Set the file extension of the satellite files being used (currently only works with netcdf/nc and h5 files).

`start_date/end_date` -Set the start and end dates for you satellite data series. The program automatically checks if they are suitable for the available model data

`x_var_name` - The following variables are the names of the variables within the satellite that the program needs (e.g. latitude, longitude etc). These should be self explanitory. If the satellite file being used does not have that product then write `None`. if the variable is within a group in the netcdf file then write it as `group_name/var_name`. 

If the time variable is split across multiple variables (e.g. year, day of year & second in day) then write these as a Python list (in square brackets, separated by commas) e.g. `['year','day_of_year','second_in_day']` and the program will try and interpret this into one date & time. 

`apply_quality_control` - Choose whether to apply the same quality control that the satellite undergo to the model. The alternative is to use all satellite observations. This may not work currently.

`model_repo_url` - This points to where the model data is being downloaded from. There is no need to change this for now

`keep_model_downloads` - Choose whether to keep the downloaded model data. This will speed up subsequent runs of the program. The model data does not take up much space currently. Higher resolution versions in the future will mean this needs changing.

`verbose` - Choose whether to print progress report of the program. This is useful to debugging but also seeing what it is doing as it can be a slow program to run depending on the satellite data.

`model_temporal_res` - This is a dummy variable currently that doesn't do anything but will in the future. Leave as is.

##### Executing the script 

The main driving script of the program is `sample_model.py`. To run this via terminal, change into the code directory and type the command `python sample_model.py`. Or run within an IDE (e.g. Spyder or Visual Studio Code). 

If the verbose keyword is set in the cofig file, you should see print outs of the progress of the program. The speed of the this program varies dramatically due to the different satellite files. The larger the satellite file the longer the script will take to run. There are plans to increase the speed of the program in future developments.

### Reporting Errors

This prorgam is still very much in the development stage and there will no doubt be errors users will come across. It is designed to be as broad as possible to account for the different types of satlellite files and this comes with difficulties in being able to make the code run smoothly.

Please report any errors or bugs on the GitHub issues page so they are logged properly and can be addressed. The code will be updated regularly so check for new updates if something isn't working.

### Future Developments

As mentioned above, there are plenty of developments and fixes needed. This will be mainly driven by the needs of the community. Please suggests any changes you want or need to be incorporated into the code.

