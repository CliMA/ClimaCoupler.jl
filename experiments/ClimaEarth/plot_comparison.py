#!/usr/bin/env python3
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import random
from pathlib import Path
from matplotlib.colors import TwoSlopeNorm
import matplotlib.animation as animation
from PIL import Image
import glob

def get_variable_file(output_dir, variable_name, file_modifier="_10m_inst.nc"):
    """
    Find the file containing the specified variable in the output directory.

    Args:
        output_dir (str): Path to the output directory
        variable_name (str): Name of the variable to find
        file_modifier (str): Suffix to append to variable name when looking for files (default: "_10m_inst.nc")

    Returns:
        str: Path to the file containing the variable, or None if not found
    """
    # Check if the variable file exists directly with the specified modifier
    var_file = os.path.join(output_dir, f"{variable_name}{file_modifier}")
    if os.path.exists(var_file):
        print(f"Found {variable_name} in file: {var_file}")
        return var_file

    print(f"Warning: Could not find {var_file}, searching through all netCDF files...")

    # If not found, search through all netCDF files
    for file in os.listdir(output_dir):
        if file.endswith(".nc"):
            file_path = os.path.join(output_dir, file)
            try:
                with xr.open_dataset(file_path) as ds:
                    if variable_name in ds.variables:
                        print(f"Found {variable_name} in alternative file: {file_path}")
                        return file_path
            except Exception as e:
                print(f"Error checking {file}: {e}")

    print(f"Warning: Variable {variable_name} not found in any netCDF file in {output_dir}")
    return None

def get_variable_ranges(output_dir, variables, time_indices, predefined_ranges={}, file_modifier="_10m_inst.nc"):
    """
    Get the min/max ranges and colormaps ONLY for variables specified in predefined_ranges.

    Args:
        output_dir (str): Path to the output directory
        variables (list): List of variable names (used to check against predefined_ranges)
        time_indices (list): List of time indices (not used for calculation here)
        predefined_ranges (dict): Dictionary mapping variable names to (vmin, vmax) tuples.
        file_modifier (str): Suffix to append to variable name when looking for files (default: "_10m_inst.nc")

    Returns:
        dict: Dictionary with variable names from predefined_ranges as keys and (vmin, vmax, colormap) tuples as values
    """
    ranges = {}

    # Only process variables that are explicitly in predefined_ranges
    for var_name in predefined_ranges.keys():
        if var_name not in variables: # Check if the predefined variable is in the list we are processing
             print(f"Info: Variable '{var_name}' from predefined_ranges is not in the current list of variables to plot. Skipping range calculation for it.")
             continue
        try:
            # We still need to find the file to potentially determine the best colormap,
            # even though vmin/vmax are predefined.
            var_file = get_variable_file(output_dir, var_name, file_modifier)
            if var_file is None:
                print(f"Warning: Could not find file for '{var_name}' listed in predefined_ranges. Cannot determine colormap.")
                continue # Skip if file not found

            # Determine colormap based on variable name convention
            # (We don't need to load data just for cmap determination if we rely on names)
            if var_name in ['waup', 'nh_pressure']:
                cmap = 'RdBu_r'
            elif var_name in ['thetaa', 'hus', 'ta', 'tke', 'taup', 'husup', 'clwvi', 'lwp', 'pr', 'hussfc']:
                cmap = 'gist_ncar'
            else:  # Assume positive-only variables like 'arup', 'clw', 'cli', 'clwup', 'cliup'
                cmap = 'viridis'

            # Use predefined ranges
            if len(predefined_ranges[var_name]) == 3:
                # If colormap is explicitly specified: (vmin, vmax, cmap)
                vmin, vmax, custom_cmap = predefined_ranges[var_name]
                cmap = custom_cmap  # Use the explicitly specified colormap
            else:
                # If only range is specified: (vmin, vmax)
                vmin, vmax = predefined_ranges[var_name]
                # Use the automatically determined colormap

            print(f"Using predefined range for {var_name}: ({vmin}, {vmax}) with cmap '{cmap}'")

            ranges[var_name] = (vmin, vmax, cmap)

        except Exception as e:
            print(f"Error processing predefined range for {var_name}: {e}")
            continue

    # Print variables that will use automatic ranging
    auto_ranged_vars = [v for v in variables if v not in predefined_ranges]
    if auto_ranged_vars:
        print(f"Variables using automatic per-level ranging: {', '.join(auto_ranged_vars)}")

    return ranges

def calculate_global_variable_ranges(output_dir, variables, time_indices, z_indices, file_modifier="_10m_inst.nc", custom_ranges=None):
    """
    Calculate global min/max ranges for all variables across all time indices and z-levels.
    This ensures consistent colorbar ranges across all frames in movies.

    Args:
        output_dir (str): Path to the output directory
        variables (list): List of variable names to process
        time_indices (list): List of time indices to include
        z_indices (list): List of z-indices to include (for 3D variables)
        file_modifier (str): Suffix to append to variable name when looking for files
        custom_ranges (dict): Dictionary of custom ranges and colormaps to override automatic determination

    Returns:
        dict: Dictionary with variable names as keys and (vmin, vmax, cmap) tuples as values
    """
    global_ranges = {}

    print("\nCalculating global variable ranges for consistent movie colorbars...")

    for variable_name in variables:
        try:
            print(f"Processing {variable_name}...")

            # Find the variable file
            var_file = get_variable_file(output_dir, variable_name, file_modifier)
            if var_file is None:
                print(f"  Warning: Could not find file for {variable_name}. Skipping.")
                continue

            # Open the dataset
            with xr.open_dataset(var_file) as ds:
                if variable_name not in ds.variables:
                    print(f"  Warning: Variable {variable_name} not found in file. Skipping.")
                    continue

                var_data = ds[variable_name]

                # Check for z dimension
                z_dim = None
                if 'z' in var_data.dims:
                    z_dim = 'z'
                elif 'z_reference' in var_data.dims:
                    z_dim = 'z_reference'

                # Initialize min/max for this variable
                var_min = float('inf')
                var_max = float('-inf')

                if z_dim is None:
                    # 2D variable - process all time indices
                    print(f"  {variable_name} is 2D. Processing {len(time_indices)} time indices...")

                    for time_idx in time_indices:
                        try:
                            if time_idx < len(ds.time):
                                data_slice = var_data.isel(time=time_idx)
                                slice_min = float(data_slice.min())
                                slice_max = float(data_slice.max())

                                var_min = min(var_min, slice_min)
                                var_max = max(var_max, slice_max)
                        except Exception as e:
                            print(f"    Warning: Could not process time index {time_idx}: {e}")
                            continue

                else:
                    # 3D variable - process all time indices and z-levels
                    print(f"  {variable_name} is 3D. Processing {len(time_indices)} time indices and {len(z_indices)} z-levels...")

                    for time_idx in time_indices:
                        for z_idx in z_indices:
                            try:
                                if time_idx < len(ds.time) and z_idx < len(ds[z_dim]):
                                    data_slice = var_data.isel(time=time_idx, **{z_dim: z_idx})
                                    slice_min = float(data_slice.min())
                                    slice_max = float(data_slice.max())

                                    var_min = min(var_min, slice_min)
                                    var_max = max(var_max, slice_max)
                            except Exception as e:
                                print(f"    Warning: Could not process time={time_idx}, z={z_idx}: {e}")
                                continue

                # Determine colormap - check custom ranges first, then use variable name convention
                if custom_ranges and variable_name in custom_ranges:
                    if len(custom_ranges[variable_name]) == 3:
                        # If colormap is explicitly specified: (vmin, vmax, cmap)
                        _, _, cmap = custom_ranges[variable_name]
                    else:
                        # If only range is specified: (vmin, vmax), use automatic determination
                        if variable_name in ['waup', 'nh_pressure']:
                            cmap = 'RdBu_r'
                        elif variable_name in ['thetaa', 'hus', 'ta', 'tke', 'taup', 'husup', 'clwvi', 'lwp', 'pr', 'hussfc']:
                            cmap = 'gist_ncar'
                        else:  # Assume positive-only variables
                            cmap = 'viridis'
                else:
                    # Use automatic colormap determination based on variable name convention
                    if variable_name in ['waup', 'nh_pressure']:
                        cmap = 'RdBu_r'
                    elif variable_name in ['thetaa', 'hus', 'ta', 'tke', 'taup', 'husup', 'clwvi', 'lwp', 'pr', 'hussfc']:
                        cmap = 'gist_ncar'
                    else:  # Assume positive-only variables
                        cmap = 'viridis'

                # Store the global range
                if var_min != float('inf') and var_max != float('-inf'):
                    global_ranges[variable_name] = (var_min, var_max, cmap)
                    print(f"  {variable_name}: global range ({var_min:.6f}, {var_max:.6f}), cmap: {cmap}")
                else:
                    print(f"  Warning: Could not determine range for {variable_name}")

        except Exception as e:
            print(f"Error processing {variable_name}: {e}")
            continue

    print(f"Global ranges calculated for {len(global_ranges)} variables.")
    return global_ranges

def plot_variable_levels(output_dir, variable_name, base_output_folder, time_index, var_ranges, z_indices=[0, 1, 2, 3, 5, 10, -10, -1], file_modifier="_10m_inst.nc"):
    """
    Plot lat/lon maps for specified z-levels of a given variable.
    Uses predefined ranges if available, otherwise calculates range per level.
    Handles both 2D and 3D variables.
    Args:
        output_dir (str): Path to the output directory
        variable_name (str): Name of the variable to plot
        base_output_folder (str): Base path to save the output plots
        time_index (int): Time index to plot (default: 0)
        var_ranges (dict): Dictionary of predefined variable ranges (vmin, vmax, cmap)
        z_indices (list): List of z-indices to plot (only used for 3D variables)
        file_modifier (str): Suffix to append to variable name when looking for files
    """
    # Find and open the file
    var_file = get_variable_file(output_dir, variable_name, file_modifier)
    if var_file is None:
        print(f"Error: Could not find file containing variable '{variable_name}'")
        return

    # Check if predefined range exists for this variable
    use_predefined_range = variable_name in var_ranges
    if use_predefined_range:
        predefined_vmin, predefined_vmax, predefined_cmap = var_ranges[variable_name]
        print(f"Plotting {variable_name} using predefined range: ({predefined_vmin}, {predefined_vmax}), cmap: {predefined_cmap}")
    else:
        print(f"Plotting {variable_name} using automatic range per level.")

    print(f"Loading data from {var_file}")
    ds = xr.open_dataset(var_file)

    # Get the variable data
    if variable_name not in ds.variables:
        print(f"Error: Could not find variable {variable_name}")
        ds.close()
        return

    var_data = ds[variable_name]

    # Check for z dimension (try 'z' first, then 'z_reference' as fallback)
    z_dim = None
    if 'z' in var_data.dims:
        z_dim = 'z'
    elif 'z_reference' in var_data.dims:
        z_dim = 'z_reference'

    # Get metadata
    var_units = getattr(var_data, 'units', '')
    var_long_name = getattr(var_data, 'long_name', variable_name)

    # Handle 2D vs 3D variables
    if z_dim is None:
        # This is a 2D variable (e.g., surface variables like pr, tas)
        print(f"Variable {variable_name} is 2D (no vertical dimension). Plotting as single 2D field.")

        # Create output folder for 2D variable
        output_folder = os.path.join(base_output_folder, variable_name, "2D_surface")
        os.makedirs(output_folder, exist_ok=True)

        # Create figure
        fig, ax = plt.subplots(figsize=(12, 8))

        # Get data slice (just time)
        data_slice = var_data.isel(time=time_index)

        # Calculate actual min and max of the slice for text annotation
        data_min = float(data_slice.min())
        data_max = float(data_slice.max())

        # Determine plot range and colormap - use global ranges if available
        norm = None
        if use_predefined_range:
            plot_vmin = predefined_vmin
            plot_vmax = predefined_vmax
            plot_cmap = predefined_cmap
        else:
            # Use global ranges for consistent colorbars across frames
            plot_vmin = var_ranges[variable_name][0] if variable_name in var_ranges else data_min
            plot_vmax = var_ranges[variable_name][1] if variable_name in var_ranges else data_max
            plot_cmap = var_ranges[variable_name][2] if variable_name in var_ranges else 'viridis'

            # Handle potential divergence around zero for RdBu_r colormap
            if plot_cmap == 'RdBu_r' and plot_vmin < 0 < plot_vmax:
                norm = TwoSlopeNorm(vmin=plot_vmin, vcenter=0, vmax=plot_vmax)

        # Plot the data slice
        im = ax.imshow(
            data_slice,
            cmap=plot_cmap,
            norm=norm,
            vmin=plot_vmin if norm is None else None,
            vmax=plot_vmax if norm is None else None,
            extent=[
                float(data_slice.lon.min()),
                float(data_slice.lon.max()),
                float(data_slice.lat.min()),
                float(data_slice.lat.max())
            ],
            origin='lower',
            aspect='auto'
        )

        # Add labels and colorbar
        plt.colorbar(im, ax=ax, label=f"{var_long_name} [{var_units}]")
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
        ax.set_title(f"{var_long_name} (2D Surface Variable)\nTime: {ds.time.values[time_index]}")

        # Determine formatting for min/max values
        if abs(data_max) < 0.01 and data_max != 0:
            max_text = f'Max: {data_max:.2e}'
        else:
            max_text = f'Max: {data_max:.2f}'

        if abs(data_min) < 0.01 and data_min != 0:
            min_text = f'Min: {data_min:.2e}'
        else:
            min_text = f'Min: {data_min:.2f}'

        # Add min/max text in the upper right corner
        ax.text(0.98, 0.98, f'{max_text}\n{min_text}',
                transform=ax.transAxes,
                verticalalignment='top',
                horizontalalignment='right',
                fontweight='bold',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='lightgray'))

        # Save and close
        output_file = os.path.join(output_folder, f"time_{time_index}.png")
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        print(f"Saved 2D plot to {output_file}")
        plt.close(fig)

    else:
        # This is a 3D variable - plot each specified z-level
        print(f"Variable {variable_name} is 3D with {z_dim} dimension. Plotting {len(z_indices)} levels.")

        for z_index in z_indices:
            try: # Add try-except for individual level processing
                z_height = float(ds[z_dim][z_index])
                output_folder = os.path.join(base_output_folder, variable_name, f"z{z_index}_{z_height:.1f}m")
                os.makedirs(output_folder, exist_ok=True)

                # Create figure
                fig, ax = plt.subplots(figsize=(12, 8))

                # Get data slice
                data_slice = var_data.isel({z_dim: z_index, 'time': time_index})

                # Calculate actual min and max of the slice for text annotation
                data_min = float(data_slice.min())
                data_max = float(data_slice.max())

                # Determine plot range and colormap - use global ranges if available
                norm = None
                if use_predefined_range:
                    plot_vmin = predefined_vmin
                    plot_vmax = predefined_vmax
                    plot_cmap = predefined_cmap
                    # Handle potential divergence around zero for RdBu_r colormap with predefined range
                    if plot_cmap == 'RdBu_r' and plot_vmin < 0 < plot_vmax:
                        norm = TwoSlopeNorm(vmin=plot_vmin, vcenter=0, vmax=plot_vmax)
                else:
                    # Use global ranges for consistent colorbars across frames
                    plot_vmin = var_ranges[variable_name][0] if variable_name in var_ranges else data_min
                    plot_vmax = var_ranges[variable_name][1] if variable_name in var_ranges else data_max
                    plot_cmap = var_ranges[variable_name][2] if variable_name in var_ranges else 'viridis'

                    # Handle potential divergence around zero for RdBu_r colormap
                    if plot_cmap == 'RdBu_r' and plot_vmin < 0 < plot_vmax:
                        norm = TwoSlopeNorm(vmin=plot_vmin, vcenter=0, vmax=plot_vmax)

                # Plot the data slice
                im = ax.imshow(
                    data_slice,
                    cmap=plot_cmap,
                    norm=norm, # Use norm for diverging colormaps if needed
                    vmin=plot_vmin if norm is None else None, # Pass vmin/vmax only if norm is not used
                    vmax=plot_vmax if norm is None else None,
                    extent=[
                        float(data_slice.lon.min()),
                        float(data_slice.lon.max()),
                        float(data_slice.lat.min()),
                        float(data_slice.lat.max())
                    ],
                    origin='lower',
                    aspect='auto'
                )

                # Add labels and colorbar
                plt.colorbar(im, ax=ax, label=f"{var_long_name} [{var_units}]")
                ax.set_xlabel('Longitude')
                ax.set_ylabel('Latitude')
                ax.set_title(f"{var_long_name} at z={z_height:.1f}m\nTime: {ds.time.values[time_index]}")

                # Determine formatting for min/max values (using slice min/max)
                if abs(data_max) < 0.01 and data_max != 0:
                    max_text = f'Max: {data_max:.2e}'
                else:
                    max_text = f'Max: {data_max:.2f}'

                if abs(data_min) < 0.01 and data_min != 0:
                    min_text = f'Min: {data_min:.2e}'
                else:
                    min_text = f'Min: {data_min:.2f}'

                # Add min/max text in the upper right corner
                ax.text(0.98, 0.98, f'{max_text}\n{min_text}',
                        transform=ax.transAxes,
                        verticalalignment='top',
                        horizontalalignment='right',
                        fontweight='bold',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='lightgray'))

                # Save and close
                output_file = os.path.join(output_folder, f"time_{time_index}.png")
                plt.savefig(output_file, bbox_inches='tight', dpi=300)
                print(f"Saved plot to {output_file}")
                plt.close(fig)

            except IndexError:
                 print(f"Warning: z-index {z_index} is out of bounds for variable {variable_name}. Skipping this level.")
                 plt.close(fig) # Ensure figure is closed if error occurs before saving
            except Exception as e:
                 print(f"Error plotting level z={z_index} for {variable_name}: {e}")
                 plt.close(fig) # Ensure figure is closed on other errors

    ds.close()

def plot_vertical_profiles(output_dir, variables, time_index, output_folder, var_ranges, file_modifier="_10m_inst.nc", max_height=10000, num_profiles=50, z_profile_max_m=15000):
    """
    Plot vertical profiles of all variables from randomly sampled lat/lon points.
    Generates three versions: log z, linear z, linear z limited.

    Args:
        output_dir (str): Path to the output directory
        variables (list): List of variable names to plot
        time_index (int): Time index to plot
        output_folder (str): Base directory to save the plots (subdirs will be created)
        var_ranges (dict): Dictionary of variable ranges
        file_modifier (str): Suffix to append to variable name when looking for files
        max_height (int): Maximum height to consider from data in meters (for initial range finding)
        num_profiles (int): Number of random profiles to sample and plot per variable
        z_profile_max_m (int): Maximum height for the limited linear z plot.
    """
    # Base output folder is now used for subdirectories
    # os.makedirs(output_folder, exist_ok=True) # Moved directory creation lower

    # --- Determine grid dimensions and sample locations ---
    nlat, nlon = None, None
    first_var_file = None
    # Find the first valid variable file to get dimensions
    for var_name in variables:
        var_file = get_variable_file(output_dir, var_name, file_modifier)
        if var_file:
            try:
                with xr.open_dataset(var_file) as ds:
                    if 'lat' in ds.dims and 'lon' in ds.dims:
                        nlat = ds.dims['lat']
                        nlon = ds.dims['lon']
                        first_var_file = var_file # Store for getting time value later
                        print(f"Determined grid dimensions from {var_name}: {nlat} lat x {nlon} lon")
                        break # Stop after finding the first valid file
            except Exception as e:
                print(f"Could not read dimensions from {var_file}: {e}")
                continue

    if nlat is None or nlon is None:
        print("Error: Could not determine grid dimensions (lat, lon) from any variable file.")
        return

    # Generate random unique lat/lon index pairs
    all_indices = [(lat, lon) for lat in range(nlat) for lon in range(nlon)]
    if num_profiles > len(all_indices):
        print(f"Warning: Requested {num_profiles} profiles, but only {len(all_indices)} grid points exist. Plotting all points.")
        sampled_indices = all_indices
        num_profiles = len(all_indices) # Adjust num_profiles to actual number
    else:
        sampled_indices = random.sample(all_indices, num_profiles)
    print(f"Sampling {num_profiles} random locations.")
    # --- End sampling ---

    # Set up the subplot grid
    n_vars = len(variables)
    n_cols = 3
    n_rows = (n_vars + n_cols - 1) // n_cols  # Ceiling division

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows), squeeze=False) # Ensure axes is always 2D
    axes = axes.flatten()  # Flatten to make indexing easier

    # Store z values to ensure consistent y-axis
    all_z_values = []
    z_coords = None # To store the actual z coordinate array

    # First pass to get z values
    for var_name in variables:
        try:
            var_file = get_variable_file(output_dir, var_name, file_modifier)
            if var_file is None:
                continue

            with xr.open_dataset(var_file) as ds:
                if var_name not in ds.variables:
                    continue

                var_data = ds[var_name]
                # Check for z dimension (try 'z' first, then 'z_reference' as fallback)
                z_dim = None
                if 'z' in var_data.dims:
                    z_dim = 'z'
                elif 'z_reference' in var_data.dims:
                    z_dim = 'z_reference'
                else:
                    continue

                all_z_values.append(ds[z_dim].values)
                if z_coords is None: # Store the first valid z coordinate array
                    z_coords = ds[z_dim].values

        except Exception as e:
            print(f"Error reading z-values for {var_name}: {e}")
            continue

    # Get common z range if we found any z values
    z_min_actual, z_max_actual = 0, max_height # Default values
    if all_z_values:
        # Use the stored z_coords for plotting if available, otherwise calculate min/max
        if z_coords is not None:
             z_min_actual = np.min(z_coords)
             z_max_actual = np.max(z_coords)
        else: # Fallback if z_coords wasn't captured (shouldn't happen if all_z_values is populated)
             z_min_actual = min(np.min(z) for z in all_z_values)
             z_max_actual = max(np.max(z) for z in all_z_values)
        # z_min = max(1, z_min_actual) # Avoid zero in log scale - Apply this later per scale type
        # z_max = min(z_max_actual, max_height) # Apply max_height limit later per scale type
    elif z_coords is None:
        print("Warning: Could not determine z coordinates for plotting.")
        # Attempt to proceed without z, plotting might fail or look incorrect
        z_coords = np.arange(10) # Placeholder
        z_min_actual, z_max_actual = 0, 10


    # Plot each variable
    plot_index = 0
    processed_vars = [] # Keep track of variables actually plotted
    for var_name in variables:
        try:
            var_file = get_variable_file(output_dir, var_name, file_modifier)
            if var_file is None:
                print(f"Skipping {var_name}: File not found.")
                continue

            with xr.open_dataset(var_file) as ds:
                if var_name not in ds.variables:
                    print(f"Skipping {var_name}: Variable not in file {var_file}.")
                    continue

                var_data = ds[var_name]
                # Check for z dimension (try 'z' first, then 'z_reference' as fallback)
                z_dim = None
                if 'z' in var_data.dims:
                    z_dim = 'z'
                elif 'z_reference' in var_data.dims:
                    z_dim = 'z_reference'
                else:
                    print(f"Skipping {var_name}: Missing z or z_reference dimension.")
                    continue

                if 'lat' not in var_data.dims or 'lon' not in var_data.dims:
                    print(f"Skipping {var_name}: Missing required dimensions (lat, lon).")
                    continue

                # Get metadata
                var_units = getattr(var_data, 'units', '')
                var_long_name = getattr(var_data, 'long_name', var_name)
                current_z_coords = ds[z_dim].values # Use z from this specific file for indexing

                # Plot profiles from all sampled locations onto the same subplot
                plotted_count = 0
                for lat_idx, lon_idx in sampled_indices:
                    # Ensure indices are valid for this specific variable's dimensions
                    if lat_idx < ds.dims['lat'] and lon_idx < ds.dims['lon']:
                        profile = var_data.isel(time=time_index, lat=lat_idx, lon=lon_idx)
                        # Use the common z_coords for the y-axis values in the plot call
                        # Ensure the profile length matches z_coords length
                        if len(profile.values) == len(current_z_coords):
                             axes[plot_index].plot(profile.values, current_z_coords, '-', linewidth=0.5, alpha=0.4)
                             plotted_count += 1
                        else:
                             print(f"Warning: Mismatch in z-dimension length for {var_name} at ({lat_idx}, {lon_idx}). Profile length {len(profile.values)}, expected {len(current_z_coords)}. Skipping this profile.")
                    else:
                         print(f"Warning: Sampled index ({lat_idx}, {lon_idx}) out of bounds for {var_name}. Skipping.")


                if plotted_count > 0:
                    axes[plot_index].set_yscale('log')
                    axes[plot_index].set_ylabel('Height (m)')
                    axes[plot_index].set_xlabel(f'{var_long_name} [{var_units}]')
                    axes[plot_index].set_title(f'{var_long_name}')
                    axes[plot_index].grid(True, which='both', linestyle='--', linewidth=0.5)

                    # Set consistent y-axis limits (will be applied later based on scale)
                    # axes[plot_index].set_ylim(z_min, z_max)

                    # Set consistent x-axis limits if we have them
                    if var_name in var_ranges:
                        vmin, vmax, _ = var_ranges[var_name]
                        axes[plot_index].set_xlim(vmin, vmax)

                    processed_vars.append(var_name) # Add to list of processed vars
                    plot_index += 1 # Move to the next subplot only if we plotted something
                else:
                    print(f"No profiles plotted for {var_name}.")


        except Exception as e:
            print(f"Error processing variable {var_name}: {e}")
            continue # Skip to the next variable

    # Remove empty subplots
    for j in range(plot_index, len(axes)):
        fig.delaxes(axes[j])

    # Add overall title with time information
    time_val = "Unknown Time"
    if first_var_file: # Use the file we stored earlier
        try:
            with xr.open_dataset(first_var_file) as ds:
                 if time_index < len(ds.time):
                     time_val = ds.time.values[time_index]
                 else:
                     time_val = f"Index {time_index} (out of bounds)"
        except Exception as e:
            print(f"Could not read time from {first_var_file}: {e}")

    fig.suptitle(f'Vertical Profiles from {num_profiles} Random Locations\nTime: {time_val}')

    # --- Save plot in three different formats ---
    plot_configs = [
        {'scale': 'log', 'subdir': 'log_z', 'z_min': max(1, z_min_actual), 'z_max': min(z_max_actual, max_height)},
        {'scale': 'linear', 'subdir': 'linear_z', 'z_min': 0, 'z_max': min(z_max_actual, max_height)},
        {'scale': 'linear', 'subdir': 'linear_z_limited', 'z_min': 0, 'z_max': min(z_max_actual, max_height, z_profile_max_m)}
    ]

    for config in plot_configs:
        current_output_folder = os.path.join(output_folder, config['subdir'])
        os.makedirs(current_output_folder, exist_ok=True)

        print(f"Applying settings for {config['subdir']}: scale={config['scale']}, ylim=({config['z_min']}, {config['z_max']})")
        for k in range(plot_index): # Apply scale and limits to all used axes
            axes[k].set_yscale(config['scale'])
            axes[k].set_ylim(config['z_min'], config['z_max'])

        # Adjust layout and save
        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust rect to make space for suptitle
        output_file = os.path.join(current_output_folder, f'vertical_profiles_sampled_time_{time_index}.pdf')
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        print(f"Saved {config['scale']} vertical profiles to {output_file}")

    plt.close(fig) # Close the figure after all saves are done

def plot_extreme_vertical_profiles(output_dir, variables, time_index, output_folder, var_ranges, target_variable, num_extreme=50, file_modifier="_10m_inst.nc", max_height=10000, z_profile_max_m=15000):
    """
    Plot vertical profiles of all variables from locations with the N largest
    and N smallest values of a target variable (across the vertical dimension).
    Generates three versions: log z, linear z, linear z limited.

    Args:
        output_dir (str): Path to the output directory
        variables (list): List of variable names to plot
        time_index (int): Time index to plot
        output_folder (str): Base directory to save the plots (subdirs will be created)
        var_ranges (dict): Dictionary of variable ranges
        target_variable (str): The variable used to find extreme locations.
        num_extreme (int): Number of largest and smallest profiles to find (total 2*N, potentially fewer if overlap).
        file_modifier (str): Suffix to append to variable name when looking for files
        max_height (int): Maximum height to consider from data in meters.
        z_profile_max_m (int): Maximum height for the limited linear z plot.
    """
    # Base output folder is now used for subdirectories
    # os.makedirs(output_folder, exist_ok=True) # Moved directory creation lower

    # --- Find extreme locations based on target_variable ---
    target_var_file = get_variable_file(output_dir, target_variable, file_modifier)
    if target_var_file is None:
        print(f"Error: Could not find file for target variable '{target_variable}'. Cannot determine extreme locations.")
        return

    print(f"Finding extreme locations based on '{target_variable}' from {target_var_file}")
    extreme_indices = set() # Use a set to automatically handle duplicates

    try:
        with xr.open_dataset(target_var_file) as ds_target:
            if target_variable not in ds_target.variables:
                print(f"Error: Target variable '{target_variable}' not found in {target_var_file}.")
                return
            # Check for z dimension (try 'z' first, then 'z_reference' as fallback)
            target_z_dim = None
            if 'z' in ds_target[target_variable].dims:
                target_z_dim = 'z'
            elif 'z_reference' in ds_target[target_variable].dims:
                target_z_dim = 'z_reference'
            else:
                print(f"Error: Target variable '{target_variable}' is missing z or z_reference dimension.")
                return

            if 'lat' not in ds_target[target_variable].dims or 'lon' not in ds_target[target_variable].dims:
                print(f"Error: Target variable '{target_variable}' is missing required dimensions (lat, lon).")
                return

            target_data = ds_target[target_variable].isel(time=time_index)
            nlat = ds_target.dims['lat']
            nlon = ds_target.dims['lon']

            # Find min/max across the vertical dimension for each lat/lon point
            # Keep coordinates to easily find indices later
            max_vals_per_col = target_data.max(dim=target_z_dim, keep_attrs=True)
            min_vals_per_col = target_data.min(dim=target_z_dim, keep_attrs=True)

            # Flatten the max/min values and get the indices of the sorted values
            # We need to handle potential NaNs which sort unpredictably
            flat_max_indices = np.argsort(np.nan_to_num(max_vals_per_col.values, nan=-np.inf).flatten())
            flat_min_indices = np.argsort(np.nan_to_num(min_vals_per_col.values, nan=np.inf).flatten())

            # Get the indices of the N largest max values
            largest_indices_flat = flat_max_indices[-num_extreme:]
            # Get the indices of the N smallest min values
            smallest_indices_flat = flat_min_indices[:num_extreme]

            # Convert flat indices back to (lat, lon) pairs
            largest_indices_2d = np.unravel_index(largest_indices_flat, (nlat, nlon))
            smallest_indices_2d = np.unravel_index(smallest_indices_flat, (nlat, nlon))

            # Add to the set of extreme indices
            for i in range(len(largest_indices_2d[0])):
                extreme_indices.add((largest_indices_2d[0][i], largest_indices_2d[1][i]))
            for i in range(len(smallest_indices_2d[0])):
                extreme_indices.add((smallest_indices_2d[0][i], smallest_indices_2d[1][i]))

            # Store time value for plot title
            time_val = "Unknown Time"
            if time_index < len(ds_target.time):
                time_val = ds_target.time.values[time_index]
            else:
                time_val = f"Index {time_index} (out of bounds)"

    except Exception as e:
        print(f"Error processing target variable {target_variable} to find extremes: {e}")
        return

    if not extreme_indices:
        print(f"No extreme indices found for {target_variable}. Skipping plotting.")
        return

    num_unique_profiles = len(extreme_indices)
    print(f"Found {num_unique_profiles} unique extreme locations based on {target_variable} (requested {2*num_extreme}).")
    sampled_indices = list(extreme_indices) # Convert set to list for iteration

    # --- Plotting setup (similar to plot_vertical_profiles) ---
    n_vars = len(variables)
    n_cols = 3
    n_rows = (n_vars + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows), squeeze=False)
    axes = axes.flatten()

    all_z_values = []
    z_coords = None

    # First pass to get z values (same logic as plot_vertical_profiles)
    for var_name in variables:
        try:
            var_file = get_variable_file(output_dir, var_name, file_modifier)
            if var_file is None: continue
            with xr.open_dataset(var_file) as ds:
                if var_name not in ds.variables: continue
                # Check for z dimension (try 'z' first, then 'z_reference' as fallback)
                z_dim = None
                if 'z' in ds[var_name].dims:
                    z_dim = 'z'
                elif 'z_reference' in ds[var_name].dims:
                    z_dim = 'z_reference'
                else:
                    continue
                all_z_values.append(ds[z_dim].values)
                if z_coords is None: z_coords = ds[z_dim].values
        except Exception as e:
            print(f"Error reading z-values for {var_name}: {e}")
            continue

    # Get common z range if we found any z values
    z_min_actual, z_max_actual = 0, max_height # Default values
    if all_z_values:
        if z_coords is not None:
             z_min_actual = np.min(z_coords)
             z_max_actual = np.max(z_coords)
        else:
             z_min_actual = min(np.min(z) for z in all_z_values)
             z_max_actual = max(np.max(z) for z in all_z_values)
        # z_min = max(1, z_min_actual) # Apply later
        # z_max = min(z_max_actual, max_height) # Apply later
    elif z_coords is None:
        print("Warning: Could not determine z coordinates for plotting.")
        z_coords = np.arange(10) # Placeholder
        z_min_actual, z_max_actual = 0, 10

    # --- Plot each variable at the extreme locations ---
    plot_index = 0
    processed_vars = []
    for var_name in variables:
        try:
            var_file = get_variable_file(output_dir, var_name, file_modifier)
            if var_file is None:
                print(f"Skipping {var_name}: File not found.")
                continue

            with xr.open_dataset(var_file) as ds:
                if var_name not in ds.variables:
                    print(f"Skipping {var_name}: Variable not in file {var_file}.")
                    continue

                var_data = ds[var_name]
                # Check for z dimension (try 'z' first, then 'z_reference' as fallback)
                z_dim = None
                if 'z' in var_data.dims:
                    z_dim = 'z'
                elif 'z_reference' in var_data.dims:
                    z_dim = 'z_reference'
                else:
                    print(f"Skipping {var_name}: Missing z or z_reference dimension.")
                    continue

                if 'lat' not in var_data.dims or 'lon' not in var_data.dims:
                    print(f"Skipping {var_name}: Missing required dimensions (lat, lon).")
                    continue

                var_units = getattr(var_data, 'units', '')
                var_long_name = getattr(var_data, 'long_name', var_name)
                current_z_coords = ds[z_dim].values

                plotted_count = 0
                for lat_idx, lon_idx in sampled_indices:
                    if lat_idx < ds.dims['lat'] and lon_idx < ds.dims['lon']:
                        profile = var_data.isel(time=time_index, lat=lat_idx, lon=lon_idx)
                        if len(profile.values) == len(current_z_coords):
                             axes[plot_index].plot(profile.values, current_z_coords, '-', linewidth=0.5, alpha=0.4)
                             plotted_count += 1
                        else:
                             print(f"Warning: Mismatch in z-dimension length for {var_name} at ({lat_idx}, {lon_idx}). Profile length {len(profile.values)}, expected {len(current_z_coords)}. Skipping this profile.")
                    else:
                         print(f"Warning: Sampled index ({lat_idx}, {lon_idx}) out of bounds for {var_name}. Skipping.")

                if plotted_count > 0:
                    axes[plot_index].set_yscale('log')
                    axes[plot_index].set_ylabel('Height (m)')
                    axes[plot_index].set_xlabel(f'{var_long_name} [{var_units}]')
                    axes[plot_index].set_title(f'{var_long_name}')
                    axes[plot_index].grid(True, which='both', linestyle='--', linewidth=0.5)
                    axes[plot_index].set_ylim(z_min_actual, z_max_actual)
                    if var_name in var_ranges:
                        vmin, vmax, _ = var_ranges[var_name]
                        axes[plot_index].set_xlim(vmin, vmax)
                    processed_vars.append(var_name)
                    plot_index += 1
                else:
                    print(f"No profiles plotted for {var_name} at extreme locations.")

        except Exception as e:
            print(f"Error processing variable {var_name} for extreme profiles: {e}")
            continue

    # --- Finalize plot ---
    for j in range(plot_index, len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle(f'Vertical Profiles from {num_unique_profiles} Locations with Extreme {target_variable.upper()} Values\nTime: {time_val}')

    # --- Save plot in three different formats ---
    plot_configs = [
        {'scale': 'log', 'subdir': 'log_z', 'z_min': max(1, z_min_actual), 'z_max': min(z_max_actual, max_height)},
        {'scale': 'linear', 'subdir': 'linear_z', 'z_min': 0, 'z_max': min(z_max_actual, max_height)},
        {'scale': 'linear', 'subdir': 'linear_z_limited', 'z_min': 0, 'z_max': min(z_max_actual, max_height, z_profile_max_m)}
    ]

    for config in plot_configs:
        current_output_folder = os.path.join(output_folder, config['subdir'])
        os.makedirs(current_output_folder, exist_ok=True)

        print(f"Applying settings for {config['subdir']}: scale={config['scale']}, ylim=({config['z_min']}, {config['z_max']})")
        for k in range(plot_index): # Apply scale and limits to all used axes
            axes[k].set_yscale(config['scale'])
            axes[k].set_ylim(config['z_min'], config['z_max'])

        # Adjust layout and save
        plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust rect to make space for suptitle
        output_file = os.path.join(current_output_folder, f'vertical_profiles_extreme_{target_variable}_time_{time_index}.pdf')
        plt.savefig(output_file, bbox_inches='tight', dpi=300)
        print(f"Saved {config['scale']} extreme vertical profiles based on {target_variable} to {output_file}")

    plt.close(fig) # Close the figure after all saves are done

def plot_negative_area_histograms(output_dir, variables, time_index=-1, output_folder=None, file_modifier="_10m_inst.nc"):
    """
    Find grid cells with negative area and plot histograms of variable values in those cells.

    Args:
        output_dir (str): Path to the output directory
        variables (list): List of variable names to plot
        time_index (int): Time index to plot
        output_folder (str): Where to save the plots
        file_modifier (str): Suffix to append to variable name when looking for files
    """
    if output_folder is None:
        output_folder = os.path.join("./plots/global/exp1", "negative_area_analysis")

    os.makedirs(output_folder, exist_ok=True)

    # Use 'arup' as the area variable
    area_file = get_variable_file(output_dir, 'arup', file_modifier)

    if area_file is None:
        print("Error: Could not find file containing 'arup' variable")
        return

    # Load the area data
    with xr.open_dataset(area_file) as ds:
        area_data = ds['arup']
        if time_index >= len(ds.time):
            time_index = -1

        # Get the time value
        time_value = ds.time.values[time_index]

        # Select the time slice and first level (where negative areas are)
        area_slice = area_data.isel(time=time_index, z=0)

        # Print some diagnostic information
        print(f"Area data range at z=0: {area_slice.min().values:.6f} to {area_slice.max().values:.6f}")

        # Find negative areas
        negative_mask = area_slice < 0
        n_negative = negative_mask.sum().item()

        if not negative_mask.any():
            print(f"No negative area values found at time index {time_index}")
            return

        # Get the indices of negative areas
        neg_indices = np.where(negative_mask)
        negative_values = area_slice.values[negative_mask]

        # Create area histogram with specified range
        plt.figure(figsize=(10, 6))
        plt.hist(negative_values, bins=50, range=(-0.03, 0))
        plt.xlabel('Area (arup)')
        plt.ylabel('Count')
        plt.title(f'Histogram of Negative Area Values at z=0\nTime: {time_value}')
        plt.grid(True)

        # Save the area histogram
        area_hist_file = os.path.join(output_folder, f'negative_area_histogram_time_{time_index}.pdf')
        plt.savefig(area_hist_file, bbox_inches='tight', dpi=300)
        plt.close()
        print(f"Saved negative area histogram to {area_hist_file}")

    # Now analyze each variable at the negative area locations
    n_vars = len(variables)
    n_cols = 3
    n_rows = (n_vars + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    axes = axes.flatten()

    # Plot histograms for each variable
    for i, var_name in enumerate(variables):
        try:
            var_file = get_variable_file(output_dir, var_name, file_modifier)
            if var_file is None:
                continue

            with xr.open_dataset(var_file) as ds:
                if var_name not in ds.variables:
                    continue

                var_data = ds[var_name]

                # Select the time slice
                if 'time' in var_data.dims:
                    var_slice = var_data.isel(time=time_index)
                else:
                    var_slice = var_data

                # Get metadata
                var_units = getattr(var_data, 'units', '')
                var_long_name = getattr(var_data, 'long_name', var_name)

                # Print diagnostic info for this variable
                print(f"\n{var_name} at negative area locations:")

                # Collect all values at negative area locations
                # Check for z dimension (try 'z' first, then 'z_reference' as fallback)
                z_dim = None
                if 'z' in var_slice.dims:
                    z_dim = 'z'
                elif 'z_reference' in var_slice.dims:
                    z_dim = 'z_reference'

                if z_dim is not None:
                    # For 3D variables, collect values from all levels
                    all_values = []
                    for z_idx in range(len(ds[z_dim])):
                        level_slice = var_slice.isel({z_dim: z_idx})
                        values = level_slice.values[neg_indices[0], neg_indices[1]]
                        all_values.extend(values.flatten())
                        print(f"  Level {z_idx}: range {values.min():.6f} to {values.max():.6f}")

                    if all_values:  # Plot if we have values
                        all_values = np.array(all_values)
                        print(f"  All levels: range {all_values.min():.6f} to {all_values.max():.6f}")
                        axes[i].hist(all_values, bins=50)
                else:
                    # For 2D variables
                    values = var_slice.values[neg_indices].flatten()
                    if len(values) > 0:
                        print(f"  2D range: {values.min():.6f} to {values.max():.6f}")
                        axes[i].hist(values, bins=50)

                axes[i].set_xlabel(f'{var_long_name} [{var_units}]')
                axes[i].set_ylabel('Count')
                axes[i].set_title(f'{var_long_name}')
                axes[i].grid(True)

        except Exception as e:
            print(f"Error plotting histogram for {var_name}: {e}")
            continue

    # Remove empty subplots
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    # Add overall title
    fig.suptitle(f'Variable Distributions in Negative Area Grid Cells\nTime: {time_value}')

    # Adjust layout and save
    plt.tight_layout()
    output_file = os.path.join(output_folder, f'negative_area_variable_histograms_time_{time_index}.pdf')
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print(f"Saved variable histograms to {output_file}")
    plt.close(fig)

def create_variable_movie(output_dir, variable_name, base_output_folder, var_ranges, z_indices, time_indices, file_modifier="_10m_inst.nc", fps=2):
    """
    Create a GIF movie showing the evolution of a variable over time at each z-level.
    Handles both 2D and 3D variables.

    Args:
        output_dir (str): Path to the output directory
        variable_name (str): Name of the variable to plot
        base_output_folder (str): Base path to save the output plots
        var_ranges (dict): Dictionary of predefined variable ranges (vmin, vmax, cmap)
        z_indices (list): List of z-indices to plot (only used for 3D variables)
        time_indices (list): List of time indices to include in the movie
        file_modifier (str): Suffix to append to variable name when looking for files
        fps (int): Frames per second for the GIF (default: 2)
    """
    # Find and open the file
    var_file = get_variable_file(output_dir, variable_name, file_modifier)
    if var_file is None:
        print(f"Error: Could not find file containing variable '{variable_name}'")
        return

    # Check if predefined range exists for this variable
    use_predefined_range = variable_name in var_ranges
    if use_predefined_range:
        predefined_vmin, predefined_vmax, predefined_cmap = var_ranges[variable_name]
        print(f"Creating movie for {variable_name} using predefined range: ({predefined_vmin}, {predefined_vmax}), cmap: {predefined_cmap}")
    else:
        print(f"Creating movie for {variable_name} using automatic range per level.")

    print(f"Loading data from {var_file}")
    ds = xr.open_dataset(var_file)

    # Get the variable data
    if variable_name not in ds.variables:
        print(f"Error: Could not find variable {variable_name}")
        ds.close()
        return

    var_data = ds[variable_name]

    # Check for z dimension (try 'z' first, then 'z_reference' as fallback)
    z_dim = None
    if 'z' in var_data.dims:
        z_dim = 'z'
    elif 'z_reference' in var_data.dims:
        z_dim = 'z_reference'

    # Get metadata
    var_units = getattr(var_data, 'units', '')
    var_long_name = getattr(var_data, 'long_name', variable_name)

    # Handle 2D vs 3D variables
    if z_dim is None:
        # This is a 2D variable - create movie from 2D_surface folder
        print(f"Variable {variable_name} is 2D. Creating movie from 2D surface plots.")

        output_folder = os.path.join(base_output_folder, variable_name, "2D_surface")

        # Check if PNG files exist for this 2D variable
        png_pattern = os.path.join(output_folder, "time_*.png")
        png_files = glob.glob(png_pattern)

        if not png_files:
            print(f"No PNG files found for 2D variable {variable_name}. Skipping movie creation.")
            ds.close()
            return

        # Sort PNG files by time index
        png_files.sort(key=lambda x: int(x.split('time_')[1].split('.')[0]))

        # Create the movie
        print(f"Creating movie for 2D variable {variable_name} with {len(png_files)} frames...")

        # Load images and create GIF
        images = []
        for png_file in png_files:
            try:
                img = Image.open(png_file)
                images.append(img)
            except Exception as e:
                print(f"Warning: Could not load {png_file}: {e}")
                continue

        if not images:
            print(f"No valid images found for 2D variable {variable_name}. Skipping movie creation.")
            ds.close()
            return

        # Save as GIF
        movie_file = os.path.join(output_folder, f"{variable_name}_2D_surface_movie.gif")

        # Calculate duration based on fps
        duration = int(1000 / fps)  # PIL uses milliseconds

        try:
            images[0].save(
                movie_file,
                save_all=True,
                append_images=images[1:],
                duration=duration,
                loop=0
            )
            print(f"Saved 2D movie to {movie_file}")
        except Exception as e:
            print(f"Error saving 2D movie for {variable_name}: {e}")

        # Close all images to free memory
        for img in images:
            img.close()

    else:
        # This is a 3D variable - create movie for each specified z-level
        print(f"Variable {variable_name} is 3D. Creating movies for {len(z_indices)} z-levels.")

        for z_index in z_indices:
            try:
                z_height = float(ds[z_dim][z_index])
                output_folder = os.path.join(base_output_folder, variable_name, f"z{z_index}_{z_height:.1f}m")

                # Check if PNG files exist for this level
                png_pattern = os.path.join(output_folder, "time_*.png")
                png_files = glob.glob(png_pattern)

                if not png_files:
                    print(f"No PNG files found for {variable_name} at z={z_height:.1f}m. Skipping movie creation.")
                    continue

                # Sort PNG files by time index
                png_files.sort(key=lambda x: int(x.split('time_')[1].split('.')[0]))

                # Create the movie
                print(f"Creating movie for {variable_name} at z={z_height:.1f}m with {len(png_files)} frames...")

                # Load images and create GIF
                images = []
                for png_file in png_files:
                    try:
                        img = Image.open(png_file)
                        images.append(img)
                    except Exception as e:
                        print(f"Warning: Could not load {png_file}: {e}")
                        continue

                if not images:
                    print(f"No valid images found for {variable_name} at z={z_height:.1f}m. Skipping movie creation.")
                    continue

                # Save as GIF
                movie_file = os.path.join(output_folder, f"{variable_name}_z{z_index}_{z_height:.1f}m_movie.gif")

                # Calculate duration based on fps
                duration = int(1000 / fps)  # PIL uses milliseconds

                try:
                    images[0].save(
                        movie_file,
                        save_all=True,
                        append_images=images[1:],
                        duration=duration,
                        loop=0
                    )
                    print(f"Saved movie to {movie_file}")
                except Exception as e:
                    print(f"Error saving movie for {variable_name} at z={z_height:.1f}m: {e}")

                # Close all images to free memory
                for img in images:
                    img.close()

            except IndexError:
                print(f"Warning: z-index {z_index} is out of bounds for variable {variable_name}. Skipping this level.")
            except Exception as e:
                print(f"Error creating movie for level z={z_index} for {variable_name}: {e}")

    ds.close()

def main():
    # Parse command line arguments
    # parser = argparse.ArgumentParser(description='Plot lat/lon maps for different vertical levels of a variable')
    # parser.add_argument('output_dir', type=str, help='Path to the output directory')
    # parser.add_argument('variable_name', type=str, help='Name of the variable to plot')
    # parser.add_argument('--output_folder', type=str, default='./plots', help='Path to save the output plots')
    # parser.add_argument('--time_index', type=int, default=0, help='Time index to plot (default: 0)')


    output_dir = "/home/jschmitt/ClimaCoupler.jl/experiments/ClimaEarth/output/wxquest_progedmf/output_0003/clima_atmos"
    # base_output_folder = "./plots/coupler_runs_v_0.7/20241201"
    base_output_folder = output_dir
    # List of variables to plot
    # variables = ['ta', 'hus', 'ua', 'va', 'waup', 'arup', 'cl', 'wa', 'clw', 'cli', 'tke', 'taup', 'edt', 'evu', 'lmix', 'pfull', 'rhoa', 'pr'] # 'nh_pressure',
    # variables = ['waup', 'wa', 'tas', 'mslp', 'pr', 'tke', 'arup', 'hussfc', 'clwvi', 'lwp', 'cli', 'ta', 'hus', 'cl', 'ua', 'va', 'waup', 'wa', 'clw', 'cli', 'pfull', 'rhoa'] # 'nh_pressure',
    variables = [ 'clwvi', 'lwp', 'cli', 'ua', 'va', 'wa',  'pfull', 'ta', 'tas', 'mslp', 'pr', 'hfss', 'hfls', 'ts', 'lwp', 'cl']

    # Define custom vmin/vmax ranges and colormaps for specific variables (optional)
    # Format: { 'variable_name': (vmin, vmax, cmap), ... } or { 'variable_name': (vmin, vmax), ... }
    custom_ranges = {
        # 'thetaa': (200, 320),   # Example: Set thetaa range from 200K to 320K
        # 'ta': (220, 320),   # Example: Set thetaa range from 200K to 320K
        # 'hus': (0, 0.02),       # Example: Set specific humidity range from 0 to 0.02
        # 'waup': (-2, 10),        # Example: Set vertical velocity range
        # 'arup': (-0.1, 0.5),    # Example: Set area range (adjust as needed)
        # 'clw': (0, 1e-3),       # Example: Set cloud liquid water range
        # 'cli': (0, 1e-4),       # Example: Set cloud ice range
        'tas': (200, 350, 'gist_ncar'),        # Set surface air temperature range from 200K to 350K with gist_ncar colormap
        'ua': (-50, 50, 'RdBu_r'),           # Set zonal wind component range from -50 to 50 m/s with diverging colormap
        'va': (-50, 50, 'RdBu_r'),           # Set meridional wind component range from -50 to 50 m/s with diverging colormap
        'wa': (-50, 50, 'RdBu_r'),           # Set vertical wind component range from -50 to 50 m/s with diverging colormap
        'lwp': (0, 1.0, 'gist_ncar'),        # Cap liquid water path at 1.0
        'clwvi': (0, 1.0, 'gist_ncar'),      # Cap column-integrated cloud water at 1.0
        'cl': (0, 100, 'gist_ncar'),
        # Add other variables here if needed
    }

    # Define the file modifier to use when searching for variables
    # file_modifier = "_10m_inst.nc"
    # file_modifier = "_1h_inst.nc"
    # file_modifier = "_10s_average.nc"
    # file_modifier = "_24h_average.nc"
    # file_modifier = "_1h_max.nc"
    file_modifier = "_1d_average.nc"

    # file_modifier = ".nc"  # Alternative simpler modifier

    z_indices=[0, 1, 2, 10, 20, 30, 40, 50, -10, -5, -2, -1]

    #time_indices = [7, 14, 21, 28, -1]
    time_indices = [-1]

    # Define max height for limited linear z profile plots
    z_profile_max_m = 150000 # meters

    # Option to create movies (GIFs)
    create_movies = False  # Set to False to skip movie creation
    movie_fps = 3  # Frames per second for the GIFs

    # Get variable ranges first - pass the custom ranges and file_modifier
    var_ranges = get_variable_ranges(output_dir, variables, time_indices, predefined_ranges=custom_ranges, file_modifier=file_modifier)

    # Calculate global ranges for consistent colorbars across all frames
    global_ranges = calculate_global_variable_ranges(output_dir, variables, time_indices, z_indices, file_modifier, custom_ranges)

    # Merge custom ranges with global ranges (custom ranges take precedence)
    for var_name, (vmin, vmax, cmap) in var_ranges.items():
        if var_name in global_ranges:
            # Update colormap from global ranges if not specified in custom ranges
            global_ranges[var_name] = (vmin, vmax, cmap)

    # Use global ranges for plotting
    final_var_ranges = global_ranges

    # ###### 2D map
    for variable_name in variables:
        for time_index in time_indices:
            try:
                plot_variable_levels(output_dir, variable_name, base_output_folder, time_index, final_var_ranges, z_indices, file_modifier)
            except Exception as e:
                print(f"Error plotting {variable_name} at time index {time_index}: {e}")

        # Create movie for this variable after all its levels are plotted
        if create_movies:
            print(f"\nCreating movie for {variable_name}...")
            try:
                create_variable_movie(
                    output_dir,
                    variable_name,
                    base_output_folder,
                    final_var_ranges,
                    z_indices,
                    time_indices,
                    file_modifier,
                    fps=movie_fps
                )
            except Exception as e:
                print(f"Error creating movie for {variable_name}: {e}")


if __name__ == '__main__':
    main()