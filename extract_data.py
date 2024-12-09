import os
import numpy as np
from netCDF4 import Dataset

def read_variables(directory):
    v10_list = []
    s10_list = []

    for root, _, files in os.walk(directory):
        for file in files:
            if file == 'cm1out.nc':
                file_path = os.path.join(root, file)
                with Dataset(file_path, 'r') as nc_file:
                    v10 = np.squeeze(nc_file.variables['v10'][:])
                    s10 = np.squeeze(nc_file.variables['s10'][:])
                v10_list.append(v10)
                s10_list.append(s10)

    return v10_list, s10_list

def read_identical_variables(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            if file == 'cm1out.nc':
                file_path = os.path.join(root, file)
                with Dataset(file_path, 'r') as nc_file:
                    r = np.squeeze(nc_file.variables['xh'][:])
                    time = np.squeeze(nc_file.variables['time'][:])
                return r, time

def save_output(v10_list, s10_list, r, time, output_file):
    with Dataset(output_file, 'w', format='NETCDF4') as nc_out:
        # Create dimensions
        nc_out.createDimension('time', len(time))
        nc_out.createDimension('r', len(r))
        nc_out.createDimension('samples', len(v10_list))

        # Create variables
        times = nc_out.createVariable('time', 'f4', ('time',))
        rs = nc_out.createVariable('r', 'f4', ('r',))
        v10s = nc_out.createVariable('v10_list', 'f4', ('samples', 'time', 'r'))
        s10s = nc_out.createVariable('s10_list', 'f4', ('samples', 'time', 'r'))

        # Assign data to variables
        times[:] = time
        rs[:] = r
        for i in range(len(v10_list)):
            v10s[i, :, :] = v10_list[i]
            s10s[i, :, :] = s10_list[i]

def main():
    directory = '/opt/data/default/kkchustu02/projects/tc_size_intensity/exps'
    output_file = 'ensemble_output.nc'

    v10_list, s10_list = read_variables(directory)
    r, time = read_identical_variables(directory)

    save_output(v10_list, s10_list, r, time, output_file)

if __name__ == "__main__":
    main()