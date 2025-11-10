#  Cyclone Montha — Bay of Bengal (November 2025)

A Python-based visualization of Cyclone **Montha (2025)** using reanalysis and best-track data.  
This project integrates **IBTrACS** cyclone metadata with **ERA5** atmospheric variables to simulate the system’s evolution in the Bay of Bengal.

---

##  Objective
To develop hands-on experience in Python-based scientific data analysis and visualize the evolution of a tropical cyclone using multi-source datasets.

---

##  Workflow

### 1. Data Extraction
- Extracted Cyclone *Montha* (storm index 401) from **IBTrACS (Last 3 Years)** dataset.
- Cleaned attributes to remove encoding errors.
- Saved subset as `montha_track.nc`.

*A preprocessing script (`1_prepare_ibtracs_montha.py`) is included to demonstrate how the Montha track was originally extracted and cleaned from the full IBTrACS dataset — running it is optional, as the processed data is already provided in `/data/`.*

### 2. Data Integration
- Loaded ERA5 variables: `u10`, `v10`, `msl`, `tp`.
- Matched ERA5 timestamps with cyclone track.
- Smoothed cyclone track using cubic interpolation for cleaner animation.

### 3. Visualization
- Dual-panel dynamic animation:
  - **Left Panel:** Rainfall (tp) + Pressure (msl) contours  
  - **Right Panel:** Wind magnitude and direction (`u10`, `v10`)
- Visualization built using **xarray**, **matplotlib**, and **cartopy**.




