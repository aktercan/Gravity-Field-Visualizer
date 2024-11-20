# Gravity-Field-Visualizer

**Gravity-Field-Visualizer** is a Python-based tool that processes gravitational potential data, analyzes spherical harmonic coefficients, and visualizes gravity anomalies. This project is tailored for geophysicists and researchers studying Earth's gravity field and its temporal variations.

---

## Features

- **Read Gravitational Potential Coefficients:**
  - Extracts spherical harmonic coefficients and metadata from `.gfc` files.
  - Handles both standard and error-included `.gfc` files.
- **Gaussian Filtering:**
  - Smooths spherical harmonic coefficients using a Gaussian filter.
  - User-defined filter radius to balance detail and noise suppression.
- **Gravity Anomaly Calculation:**
  - Computes gravity disturbances using spherical harmonics and Legendre polynomials.
  - Supports user-defined maximum degree for harmonic analysis.
- **Global Visualization:**
  - Displays gravity anomalies on a world map.
  - Optionally overlays coastlines for geographical context.
- **Customizable:**
  - Easily adapt the code for various `.gfc` datasets and visualization preferences.

## Configure Parameters

The script uses these default parameters:
    •    First File: EGSIEM_COMB_90_NEQ_2007_03.gfc
    •    Second File: EGSIEM_COMB_90_NEQ_2007_09.gfc
    •    Maximum Degree: 30
    •    Filter Radius: 500 km

You can modify these directly in the main() function in main.py.

## View Results

    •    The computed anomalies will be displayed as a contour map.
    •    If coastline data (coast.dat) is available, it will be overlayed.
    
## Detailed Functionality

1. Reading Gravitational Coefficients

    •    Function: read_potential_coefficients(filename)
    •    Reads spherical harmonic coefficients (Cnm, Snm) and metadata such as:
    •    Maximum degree
    •    Earth’s gravity constant (GM)
    •    Reference radius (R)
    •    Handles error data if present in the .gfc file.

2. Gaussian Filtering

    •    Function: filter_coefficients_gaussian(radius, max_degree)
    •    Applies a Gaussian filter to suppress high-frequency noise.
    •    Customizable filter radius:
    •    Small radius → Retains high-frequency details.
    •    Large radius → More smoothing.

3. Gravity Anomaly Calculation

    •    Function: gravity_disturbance(Lambda, Theta, cnm, snm, GM, R, max_degree)
    •    Computes gravity disturbances based on spherical harmonics and Legendre polynomials.
    •    Outputs a 2D grid of gravity anomalies.

4. Visualization

    •    Function: show_grid(lambda_vals, theta_vals, grid, titletext)
    •    Displays gravity anomalies as a global contour map.
    •    Features:
    •    Longitude (x-axis) and latitude (y-axis).
    •    Adjustable resolution and color mapping.
    •    Optional coastline overlay from coast.dat.

## Example Data

Sample .gfc files provided in the DATA/ directory:
    •    EGSIEM_COMB_90_NEQ_2007_03.gfc (March 2007)
    •    EGSIEM_COMB_90_NEQ_2007_09.gfc (September 2007)

These datasets allow comparison of gravity field changes over six months.

## Required Libraries

    •    numpy - For numerical calculations.
    •    matplotlib - For visualization.
    

