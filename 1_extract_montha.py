import xarray as xr
import re

# --- 1. Load dataset
file_path = "/IBTrACS.last3years.v04r01.nc"
ds = xr.open_dataset(file_path)

# --- 2. Extract only Montha (storm index 401)
montha_ds = ds.isel(storm=[401])

# --- 3. Clean attributes to remove problematic Unicode
def clean_attrs(ds):
    for var in ds.variables:
        for attr in list(ds[var].attrs):
            val = ds[var].attrs[attr]
            if isinstance(val, str):
                ds[var].attrs[attr] = re.sub(r'[^\x00-\x7F]+', '', val)
    for attr in list(ds.attrs):
        val = ds.attrs[attr]
        if isinstance(val, str):
            ds.attrs[attr] = re.sub(r'[^\x00-\x7F]+', '', val)
    return ds

montha_ds = clean_attrs(montha_ds)

# --- 4. Save as NetCDF using SciPy engine (offline safe)
montha_ds.to_netcdf("montha_track.nc", engine="scipy")


print("✅ Saved cleaned NetCDF: montha_track.nc")
