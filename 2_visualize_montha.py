import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from matplotlib.animation import FuncAnimation
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

# --------------------------
# Load data
# --------------------------
track_ds = xr.open_dataset("data/sample_ibtracs_montha_subset.nc")
track_lat = track_ds['lat'].values[0]
track_lon = track_ds['lon'].values[0]
track_time = track_ds['time'].values[0]
valid_idx = ~np.isnan(track_lat) & ~np.isnan(track_lon)
track_lat, track_lon, track_time = track_lat[valid_idx], track_lon[valid_idx], track_time[valid_idx]

era_ds = xr.open_dataset("data/sample_era5_instant.nc")
precip_ds = xr.open_dataset("data/sample_era5_accum.nc")
u10, v10, msl, tp = era_ds['u10'], era_ds['v10'], era_ds['msl'], precip_ds['tp']

# --------------------------
# Smooth the track
# --------------------------
n_interp = len(track_lat) * 4
f_lat = interp1d(np.arange(len(track_lat)), track_lat, kind='cubic')
f_lon = interp1d(np.arange(len(track_lon)), track_lon, kind='cubic')
track_lat_smooth = f_lat(np.linspace(0, len(track_lat)-1, n_interp))
track_lon_smooth = f_lon(np.linspace(0, len(track_lon)-1, n_interp))
total_duration = track_time[-1] - track_time[0]
step = total_duration / (n_interp - 1)
track_time_smooth = np.array([track_time[0] + i * step for i in range(n_interp)], dtype='datetime64[s]')

# --------------------------
# Custom colormaps
# --------------------------


precip_colors = ["#e0f3f8", "#91bfdb", "#4575b4", "#fee090", "#fc8d59", "#d73027"]
precip_cmap = LinearSegmentedColormap.from_list("wmo_rain", precip_colors)
precip_norm = PowerNorm(gamma=0.6, vmin=0, vmax=25)

wind_cmap = plt.cm.get_cmap("turbo")

# --------------------------
# Map setup
# --------------------------
def setup_map(ax):
    ax.set_extent([70, 100, 5, 30], crs=ccrs.PlateCarree())
    india_shp = "data/India_State_Boundary.shp"
    coast_shp = "data/ne_50m_coastline.shp"
    india_reader = shpreader.Reader(india_shp)
    coast_reader = shpreader.Reader(coast_shp)
    ax.add_geometries(india_reader.geometries(), ccrs.PlateCarree(),
                      edgecolor='black', facecolor='none', linewidth=0.8)
    ax.add_geometries(coast_reader.geometries(), ccrs.PlateCarree(),
                      edgecolor='gray', facecolor='none', linewidth=0.7)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, color='gray', alpha=0.4, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    ax.text(78, 20, 'India', fontsize=13, fontweight='bold')
    ax.text(87, 18, 'Bay of Bengal', fontsize=12, fontstyle='italic')
    ax.text(93, 13, 'Andaman Islands', fontsize=11, fontstyle='italic')

# --------------------------
# Animation function
# --------------------------
def animate_cyclone(track_lat, track_lon, track_time, u10, v10, msl, tp):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9),
                                   subplot_kw={'projection': ccrs.PlateCarree()})
    setup_map(ax1)
    setup_map(ax2)
    ax1.set_title("Rainfall + Pressure Contours", fontsize=15, y=1.06)
    ax2.set_title("Wind Magnitude + Wind Vectors", fontsize=15, y=1.06)

    times_era = u10['valid_time'].values
    lat_grid = u10['latitude'].values
    lon_grid = u10['longitude'].values

    pressure_plot = precip_plot = wind_mag_plot = quiver_plot = None
    track_line, = ax2.plot([], [], 'r--', lw=1.5, transform=ccrs.PlateCarree())
    point, = ax2.plot([], [], 'ro', markersize=5, transform=ccrs.PlateCarree())

    # Colorbar scaffolds
    first_idx = np.argmin(np.abs(times_era - track_time[0]))
    p0 = msl.isel(valid_time=first_idx).values / 100
    pr0 = tp.isel(valid_time=first_idx).values * 1000
    ws0 = np.sqrt(u10.isel(valid_time=first_idx)**2 + v10.isel(valid_time=first_idx)**2).values

    # Colorbars
    sm_rain = plt.cm.ScalarMappable(cmap=precip_cmap)
    sm_rain.set_clim(vmin=0, vmax=19)
    fig.colorbar(sm_rain, ax=ax1, orientation='horizontal', fraction=0.05,
                 pad=0.07, label='Rainfall (mm / 3hr)')

    sm_wind = plt.cm.ScalarMappable(cmap=wind_cmap)
    sm_wind.set_clim(vmin=0, vmax=25)
    fig.colorbar(sm_wind, ax=ax2, orientation='horizontal', fraction=0.05,
                 pad=0.07, label='Wind Speed (m/s)')

    # --------------------------
    # Frame update
    # --------------------------
    def update(frame):
        nonlocal pressure_plot, precip_plot, wind_mag_plot, quiver_plot
        t = track_time[frame]
        era_idx = np.argmin(np.abs(times_era - t))

        u = u10.isel(valid_time=era_idx).values
        v = v10.isel(valid_time=era_idx).values
        p = msl.isel(valid_time=era_idx).values / 100
        pr = tp.isel(valid_time=era_idx).values * 1000
        ws = np.sqrt(u**2 + v**2)

        # ---- PANEL 1 ----
        if precip_plot is not None:
            precip_plot.remove()
        if pressure_plot is not None:
            for c in pressure_plot.collections:
                c.remove()
        if hasattr(ax1, "_label_artists"):
            for lbl in ax1._label_artists:
                lbl.remove()

        for artist in getattr(ax1, "_frame_artists", []):
            artist.remove()
        ax1._frame_artists = []
        ax1._label_artists = []  # reset label tracker

        # Rainfall shading
        precip_plot = ax1.pcolormesh(
            lon_grid, lat_grid, pr,
            cmap=precip_cmap, norm=precip_norm,
            alpha=0.9, transform=ccrs.PlateCarree()
        )

        

        # Smooth pressure slightly
        p_smooth = gaussian_filter(p, sigma=1.2)

        # Pressure contours
        pressure_plot = ax1.contour(
            lon_grid, lat_grid, p_smooth,
            levels=np.arange(990, 1036, 4),
            colors='red', linewidths=1.0, alpha=0.85,
            transform=ccrs.PlateCarree()
        )
        # Store labels for cleanup
        labels = ax1.clabel(
            pressure_plot, inline=True, fontsize=8, fmt='%1.0f',
            colors='darkred', inline_spacing=4
        )
        ax1._label_artists.extend(labels)

        # Cyclone track + red cyclone center
        track_obj = ax1.plot(track_lon[:frame], track_lat[:frame],
                             'r--', lw=1.2, alpha=0.6, transform=ccrs.PlateCarree())[0]
        dot_obj = ax1.plot(track_lon[frame], track_lat[frame],
                             'ro', markersize=6, markeredgecolor='black',transform=ccrs.PlateCarree(), zorder=10)[0]
        ax1._frame_artists.extend([track_obj, dot_obj])

        # ---- PANEL 2 ----
        if wind_mag_plot is not None:
            wind_mag_plot.remove()
        if quiver_plot is not None:
            quiver_plot.remove()

        wind_mag_plot = ax2.pcolormesh(
            lon_grid, lat_grid, ws,
            cmap=wind_cmap, vmin=0, vmax=25,
            transform=ccrs.PlateCarree()
        )
        quiver_plot = ax2.quiver(
            lon_grid[::7], lat_grid[::7],
            u[::7, ::7], v[::7, ::7],
            color='black', scale=500, width=0.002, alpha=0.6,
            transform=ccrs.PlateCarree()
        )

        track_line.set_data(track_lon[:frame], track_lat[:frame])
        point.set_data([track_lon[frame]], [track_lat[frame]])

        fig.suptitle(f"Cyclone Montha | {np.datetime_as_string(t, unit='h')}",
                     fontsize=20, fontweight='bold', y=0.93)

        return precip_plot, pressure_plot, wind_mag_plot, quiver_plot, track_line, point




    anim = FuncAnimation(fig, update, frames=n_interp, blit=False)
    anim.save("montha.mp4", fps=8, dpi=150)
    plt.close()

# --------------------------
# Run
# --------------------------

animate_cyclone(track_lat_smooth, track_lon_smooth, track_time_smooth, u10, v10, msl, tp)

