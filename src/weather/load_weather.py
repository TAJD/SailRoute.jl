
ad weather files.

Reduce the amount of interpolating associated with large weather files.

Thomas Dickson
thomas.dickson@soton.ac.uk
25/05/2018
"""

import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe
import glob, os,  time
from datetime import datetime


def look_in_netcdf(path):
    """Load netcdf file and return xrray."""
    with xr.open_dataset(path) as ds:
        print(ds.keys())


def load_dataset(path_nc, var):
    """Load netcdf file and return a specific variable."""
    with xr.open_dataset(path_nc) as ds:
        ds.coords['lat'] = ('latitude', ds['latitude'].values)
        ds.coords['lon'] = ('longitude', ds['longitude'].values)
        ds.swap_dims({'longitude': 'lon', 'latitude': 'lat'})
        return ds[var]


def load_dataset_ensemble(path_nc, var, ens):
    """Load dataset for a specific ensemble scenario from ERA5 netcdf file."""
    with xr.open_dataset(path_nc) as ds:
        ds.coords['lat'] = ('latitude', ds['latitude'].values)
        ds.coords['lon'] = ('longitude', ds['longitude'].values)
        ds.swap_dims({'longitude': 'lon', 'latitude': 'lat'})
        ds = ds.isel(number=ens)
        return ds[var]


def return_time_index(ds):
    """Return the time index from a Dataset in unix timestamp."""
    time_vals = ds['time'].values
    times = [(t - np.datetime64('1970-01-01T00:00:00Z'))/np.timedelta64(1, 's') for t in time_vals]
    return np.array(times)


def regrid_data(ds, longs, lats):
    """Regrid dataset to new longs and lats."""
    ds_out = xr.Dataset({'lat': (['lat_b'], lats),
                         'lon': (['lon_b'], longs)})
    regridder = xe.Regridder(ds, ds_out, 'patch', reuse_weights=True)
    ds0 = regridder(ds)
    ds0.coords['lat_b'] = ('lat_b', ds0['lat'].values)
    ds0.coords['lon_b'] = ('lon_b', ds0['lon'].values)
    return ds0


def return_data(ds):
    """Return the values, longs, lats and timestamps of dataset."""
    ds = ds.bfill('longitude', limit=None)
    values = ds.values
    lons = ds['lon'].values
    lats = ds['lat'].values
    time_vals = ds['time'].values
    return values, lons, lats


def load_cluster(path_nc, longs, lats, var):
    """Load cluster data."""
    with xr.open_dataset(path_nc) as ds:
            ds.coords['lat'] = ('latitude', ds['latitude'].values)
            ds.coords['lon'] = ('longitude', ds['longitude'].values)
            ds.swap_dims({'longitude': 'lon', 'latitude': 'lat'})
    x = xr.DataArray(longs, dims='lon')
    y = xr.DataArray(lats, dims='lat')
    ds = ds.to_array(var)
    return ds.interp(longitude=x, latitude=y)


def process_wind(path_nc, longs, lats):
    """
    Return wind speed and direction data.

    Data is regridded to the location of each node.
    """
    ds_u10 = load_dataset(path_nc, 'u10')
    regrid_ds_u10 = regrid_data(ds_u10, longs[:, 0], lats[0, :])
    ds_v10 = load_dataset(path_nc, 'v10')
    regrid_ds_v10 = regrid_data(ds_v10, longs[:, 0], lats[0, :])
    ws = 1.943844 * (regrid_ds_u10**2 + regrid_ds_v10**2)**0.5
    wind_dir = np.rad2deg(np.arctan2(regrid_ds_u10, regrid_ds_v10)) + 180.0
    return ws, wind_dir


def process_waves(path_nc, longs, lats):
    """Return wave data."""
    wh = load_dataset(path_nc, 'swh')
    wd = load_dataset(path_nc, 'mwd')
    wp = load_dataset(path_nc, 'mwp')
    regrid_wh = regrid_data(wh, longs[:, 0], lats[0, :])
    regrid_wd = regrid_data(wd, longs[:, 0], lats[0, :])
    regrid_wp = regrid_data(wp, longs[:, 0], lats[0, :])
    return regrid_wh, regrid_wd, regrid_wp


def process_era5_weather(path_nc, longs, lats):
    """Return era5 weather data."""
    wisp = load_dataset(path_nc, 'wind')
    widi = load_dataset(path_nc, 'dwi')
    wh = load_dataset(path_nc, 'swh')
    wd = load_dataset(path_nc, 'mdts')
    wp = load_dataset(path_nc, 'mpts')
    rg_wisp = regrid_data(wisp, longs[:, 0], lats[0, :])
    rg_widi = regrid_data(widi, longs[:, 0], lats[0, :])
    rg_wh = regrid_data(wh, longs[:, 0], lats[0, :])
    rg_wd = regrid_data(wd, longs[:, 0], lats[0, :])
    rg_wp = regrid_data(wp, longs[:, 0], lats[0, :])
    wisp = None
    widi = None
    wh = None
    wd = None
    wp = None
    return rg_wisp, rg_widi, rg_wh, rg_wd, rg_wp


def retrieve_era5_ensemble(path_nc, ens=0):
    """Retrieve the weather data for a specific ensemble within an ERA5 dataset."""
    ds_u10 = load_dataset_ensemble(path_nc, 'u10', ens)
    ds_v10 = load_dataset_ensemble(path_nc, 'v10', ens)
    # print(np.count_nonzero(~np.isnan(ds_u10.data)))
    wisp = 1.943844 * (ds_u10**2 + ds_v10**2)**0.5
    # print(np.count_nonzero(~np.isnan(wisp.data))) # same number of nans before and after weather unit conversion
    widi = np.rad2deg(np.arctan2(ds_u10, ds_v10)) + 180.0
    wh = load_dataset_ensemble(path_nc, 'swh', ens)
    wd = load_dataset_ensemble(path_nc, 'mwd', ens)
    wp = load_dataset_ensemble(path_nc, 'mwp', ens)
    time = return_time_index(wh)
    return wisp, widi, wh, wd, wp, time


def test_load_ensemble_scenario():
    """Test the loading of an ensemble scenario.
    1. Load weather data
    2. Select variable
    3. Select ensemble
    4. Return ensemble data
    5. Convert u10 and v10 to wind speed and wind direction.
    """
    path = "/scratch/td7g11/era5/polynesia_2005_01.nc"
    look_in_netcdf(path)
    ds = load_dataset_ensemble(path, "u10", 7)
    print(ds)
    wisp, widi, wh, wd, wp, time = retrieve_era5_ensemble(path, 7)
    # print(wisp)
    print(np.count_nonzero(~np.isnan(wisp)))
    print(np.count_nonzero(~np.isnan(widi)))
    print(np.count_nonzero(~np.isnan(wh)))
    print(np.count_nonzero(~np.isnan(wd)))
    print(np.count_nonzero(~np.isnan(wp)))


def retrieve_era20_weather(path_nc):
    wisp = load_dataset(path_nc, 'wind')
    widi = load_dataset(path_nc, 'dwi')
    wh = load_dataset(path_nc, 'swh')
    wd = load_dataset(path_nc, 'mwd')
    wp = load_dataset(path_nc, 'mwp')
    time = return_time_index(wisp)
    return wisp, widi, wh, wd, wp, time


def change_area_values(array, value, lon1, lat1, lon2, lat2):
    """
    Change the weather values in a given rectangular area.

    array is an xarray DataArray
    value is the new value
    lon1 and lat1 are the coordinates of the bottom left corner of the area
    lon2 and lat2 are the coordinates of the top right of the area
    """
    lc = array.coords['lon']
    la = array.coords['lat']
    array.loc[dict(lon_b=lc[(lc > lon1) & (lc < lon2)],
                   lat_b=la[(la > lat1) & (la < lat2)])] = value
    return array


def sample_weather_scenario():
    """
    Generate a weather scenario with known values for the wind condition.
    """
    times = pd.date_range('1/1/2000', periods=72, freq='6H')
    latitude = np.linspace(0, 10, 11)
    longitude = np.linspace(0, 10, 11)
    wsp_vals = np.full((72, 11, 11), 10.0)
    wdi_vals = np.full((72, 11, 11), 0.0)
    cusp_vals = np.full((72, 11, 11), 0.0)
    cudi_vals = np.full((72, 11, 11), 0.0)
    wadi_vals = np.full((72, 11, 11), 0.0)
    wahi_vals = np.full((72, 11, 11), 0.0)
    wisp = xr.DataArray(wsp_vals, dims=['time', 'lon_b', 'lat_b'],
                        coords={'time': times,
                                'lon_b': longitude,
                                'lat_b': latitude})
    widi = xr.DataArray(wdi_vals, dims=['time', 'lon_b', 'lat_b'],
                        coords={'time': times,
                                'lon_b': longitude,
                                'lat_b': latitude})
    cusp = xr.DataArray(cusp_vals, dims=['time', 'lon_b', 'lat_b'],
                        coords={'time': times,
                                'lon_b': longitude,
                                'lat_b': latitude})
    cudi = xr.DataArray(cudi_vals, dims=['time', 'lon_b', 'lat_b'],
                        coords={'time': times,
                                'lon_b': longitude,
                                'lat_b': latitude})
    wahi = xr.DataArray(cusp_vals, dims=['time', 'lon_b', 'lat_b'],
                        coords={'time': times,
                                'lon_b': longitude,
                                'lat_b': latitude})
    wadi = xr.DataArray(cudi_vals, dims=['time', 'lon_b', 'lat_b'],
                        coords={'time': times,
                                'lon_b': longitude,
                                'lat_b': latitude})
    return wisp, widi, cusp, cudi, wahi, wadi


def get_weather_files(dir_path):
    return glob.glob(dir_path+"/*.nc")


def concatenate_weather_files(dir_path):
    """Concatenate all .nc files found in the directory set by path."""
    # import all the files as datasets
    fnames = get_weather_files(dir_path)
    ds_list = []
    for f in fnames:
        with xr.open_dataset(f, engine='netcdf4') as ds:
            ds_list.append(ds)
    ds_main = xr.concat(ds_list, dim='time')
    groups = ds_main.groupby('time')
    return groups


def aggregate_weather_files():
    # weather_path = os.path.dirname(os.path.realpath(__file__))
    path = "/mainfs/home/td7g11/weather_data/polynesia_weather/low/1976"
    # path = weather_path + "/weather_data/polynesia_weather/low/1976/"
    # get a list of all the files in a directory
    ds_whole = concatenate_weather_files(path)
    print(ds_whole.last())
    ds_whole.last().to_netcdf(path+"1976_polynesia.nc")


def aggregate_era5_files_example():
    path = "/scratch/td7g11/era5/2005_q2"
    ds_whole = concatenate_weather_files(path)
    print(ds_whole.last())
    ds_whole.last().to_netcdf("/scratch/td7g11/era5/2005_q2/"+"polynesia_2005_q2.nc")


def aggregate_era5_files(name):
    """Function to aggregate individual era5 files."""
    path = "/scratch/td7g11/era5/polynesia_"+str(name)+"/"
    ds_whole = concatenate_weather_files(path)
    ds_whole.last().to_netcdf("/scratch/td7g11/era5/polynesia_"+str(name)+"/polynesia_"+str(name)+".nc")
    print(str(name)+" concatenated")


if __name__ == '__main__':
    # test_load_ensemble_scenario()
    years = ["1997", "1998"]
    names = [year+"_q"+str(q) for year in years for q in range(1,5)]
#     aggregate_era5_files("1995_q4")
#     aggregate_era5_files("1996_q4")
#     aggregate_era5_files("2004_q4")
#     aggregate_era5_files("2005_q4")
    for n in names:
        aggregate_era5_files(n)
