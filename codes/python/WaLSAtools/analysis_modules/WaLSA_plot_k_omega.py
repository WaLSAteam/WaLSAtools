# -----------------------------------------------------------------------------------------------------
# WaLSAtools - Wave analysis tools
# Copyright (C) 2025 WaLSA Team - Shahin Jafarzadeh et al.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# Note: If you use WaLSAtools for research, please consider citing:
# Jafarzadeh, S., Jess, D. B., Stangalini, M. et al. 2025, Nature Reviews Methods Primers, 5, 21
# -----------------------------------------------------------------------------------------------------
# The following codes are baed on those originally written by David B. Jess and Samuel D. T. Grant
# -----------------------------------------------------------------------------------------------------

import numpy as np # type: ignore
from matplotlib.colors import ListedColormap # type: ignore
from WaLSAtools import WaLSA_histo_opt # type: ignore
from scipy.interpolate import griddata # type: ignore
import matplotlib.pyplot as plt # type: ignore
from matplotlib.ticker import FixedLocator, FixedFormatter # type: ignore
import matplotlib.patches as patches # type: ignore
import matplotlib.patheffects as path_effects # type: ignore
from matplotlib.colors import Normalize # type: ignore

def WaLSA_plot_k_omega(
    kopower, kopower_xscale, kopower_yscale, 
    xtitle='Wavenumber (pixel⁻¹)', ytitle='Frequency (mHz)',
    xlog=False, ylog=False, xrange=None, yrange=None,
    xtick_interval=0.05, ytick_interval=200,
    xtick_minor_interval=0.01, ytick_minor_interval=50,
    colorbar_label='Log₁₀(Power)', cmap='viridis', cbartab=0.12,
    figsize=(9, 4), smooth=True, dypix = 1200, dxpix = 1600, ax=None,
    k1=None, k2=None, f1=None, f2=None, colorbar_location='right'
):
    """
    Plots the k-omega diagram with corrected orientation and axis alignment.
    Automatically clips the power array and scales based on x and y ranges.

    Example usage:
    WaLSA_plot_k_omega(
        kopower=power,
        kopower_xscale=wavenumber,
        kopower_yscale=frequencies*1000.,
        xlog=False, ylog=False, 
        xrange=(0, 0.3), figsize=(8, 4), cbartab=0.18,
        xtitle='Wavenumber (pixel⁻¹)', ytitle='Frequency (mHz)',
        colorbar_label='Log₁₀(Oscillation Power)',
        f1=0.470*1000, f2=0.530*1000,
        k1=0.047, k2=0.25
    )
    """

    if xrange is None or len(xrange) != 2:
        xrange = [np.min(kopower_xscale), np.max(kopower_xscale)]

    if yrange is None or len(yrange) != 2:
        yrange = [np.min(kopower_yscale), np.max(kopower_yscale)]

    # Handle xrange and yrange adjustments
    if xrange is not None and len(xrange) == 2 and xrange[0] == 0:
        xrange = [np.min(kopower_xscale), xrange[1]]

    if yrange is not None and len(yrange) == 2 and yrange[0] == 0:
        yrange = [np.min(kopower_yscale), yrange[1]]

    # Clip kopower, kopower_xscale, and kopower_yscale based on xrange and yrange
    if xrange:
        x_min, x_max = xrange
        x_indices = np.where((kopower_xscale >= x_min) & (kopower_xscale <= x_max))[0]

        # Check if the lower bound is not included
        if kopower_xscale[x_indices[0]] > x_min and x_indices[0] > 0:
            x_indices = np.insert(x_indices, 0, x_indices[0] - 1)

        # Check if the upper bound is not included
        if kopower_xscale[x_indices[-1]] < x_max and x_indices[-1] + 1 < len(kopower_xscale):
            x_indices = np.append(x_indices, x_indices[-1] + 1)

        kopower = kopower[:, x_indices]
        kopower_xscale = kopower_xscale[x_indices]

    if yrange:
        y_min, y_max = yrange
        y_indices = np.where((kopower_yscale >= y_min) & (kopower_yscale <= y_max))[0]

        # Check if the lower bound is not included
        if kopower_yscale[y_indices[0]] > y_min and y_indices[0] > 0:
            y_indices = np.insert(y_indices, 0, y_indices[0] - 1)

        # Check if the upper bound is not included
        if kopower_yscale[y_indices[-1]] < y_max and y_indices[-1] + 1 < len(kopower_yscale):
            y_indices = np.append(y_indices, y_indices[-1] + 1)

        kopower = kopower[y_indices, :]
        kopower_yscale = kopower_yscale[y_indices]

    # Pixel coordinates in data space
    if xlog:
        nxpix = np.logspace(np.log10(xrange[0]), np.log10(xrange[1]), dxpix)
    else:
        nxpix = np.linspace(xrange[0], xrange[1], dxpix)

    if ylog:
        nypix = np.logspace(np.log10(yrange[0]), np.log10(yrange[1]), dypix)
    else:
        nypix = np.linspace(yrange[0], yrange[1], dypix)

    # Interpolation over data
    interpolator = griddata(
        points=(np.repeat(kopower_yscale, len(kopower_xscale)), np.tile(kopower_xscale, len(kopower_yscale))),
        values=kopower.ravel(),
        xi=(nypix[:, None], nxpix[None, :]),
        method='linear'
    )
    newimage = np.nan_to_num(interpolator)

    # Clip the interpolated image to the exact xrange and yrange
    x_indices = (nxpix >= xrange[0]) & (nxpix <= xrange[1])
    y_indices = (nypix >= yrange[0]) & (nypix <= yrange[1])
    newimage = newimage[np.ix_(y_indices, x_indices)]
    nxpix = nxpix[x_indices]
    nypix = nypix[y_indices]

    # Extent for plotting
    extent = [nxpix[0], nxpix[-1], nypix[0], nypix[-1]]

    # Set up the plot
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    # Keep the examples' IDL palette when it is available. An installed
    # package can also plot from a directory without the example assets.
    palette = cmap
    if cmap == 'viridis':
        try:
            rgb_values = np.loadtxt('Color_Tables/idl_colormap_13.txt') / 255.0
        except FileNotFoundError:
            pass
        else:
            palette = ListedColormap(rgb_values)

    # Compute minimum and maximum values of the data
    vmin = np.nanmin(kopower) 
    vmax = np.nanmax(kopower) 

    img = ax.imshow(
        WaLSA_histo_opt(newimage), extent=extent, origin='lower', aspect='auto', cmap=palette, norm=Normalize(vmin=vmin, vmax=vmax)
    )

    # Configure axis labels and scales
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    if xlog:
        ax.set_xscale('log')
    if ylog:
        ax.set_yscale('log')

    # Configure ticks
    ax.xaxis.set_major_locator(plt.MultipleLocator(xtick_interval))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(xtick_minor_interval))
    ax.yaxis.set_major_locator(plt.MultipleLocator(ytick_interval))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(ytick_minor_interval))

    # # Add a secondary x and y axes for perido and spatial size
    ax.tick_params(axis='both', which='both', direction='out', top=True, right=True)
    ax.tick_params(axis='both', which='major', length=7, width=1.1)  # Major ticks
    ax.tick_params(axis='both', which='minor', length=4, width=1.1)  # Minor ticks

    major_xticks = ax.get_xticks()
    major_yticks = ax.get_yticks()
    
    # Define a function to calculate the period (1000 / frequency)
    def frequency_to_period(frequency):
        return 1000. / frequency if frequency > 0 else np.inf
    # Generate labels for the secondary y-axis based on the primary y-axis ticks
    period_labels = [f"{frequency_to_period(tick):.1f}" if tick > 0 else "" for tick in major_yticks]

    # Set custom labels for the secondary y-axis
    ax_right = ax.secondary_yaxis('right')
    ax_right.set_yticks(major_yticks)
    ax_right.set_yticklabels(period_labels)
    ax_right.set_ylabel('Period (s)', labelpad=12)

    # Define a function to calculate spatial size (2π / wavenumber)
    def wavenumber_to_spatial_size(wavenumber):
        return 2 * np.pi / wavenumber if wavenumber > 0 else np.inf

    # Generate labels for the secondary x-axis based on the primary x-axis ticks
    spatial_size_labels = [f"{wavenumber_to_spatial_size(tick):.1f}" if tick > 0 else "" for tick in major_xticks]

    # Set custom labels for the secondary x-axis
    ax_top = ax.secondary_xaxis('top')
    ax_top.set_xticks(major_xticks)
    ax_top.set_xticklabels(spatial_size_labels)
    ax_top.set_xlabel('Spatial Size (pixel)', labelpad=12)

    if k1 is not None and k2 is not None and f1 is not None and f2 is not None:
        width = k2 - k1
        height = f2 - f1
        rectangle = patches.Rectangle(
            (k1, f1), width, height, zorder=10,
            linewidth=1.5, edgecolor='white', facecolor='none', linestyle='--'
        )
        rectangle.set_path_effects([
            path_effects.Stroke(linewidth=2.5, foreground='black'),
            path_effects.Normal()
        ])
        # Mark the filtered area
        ax.add_patch(rectangle)

    # Add colorbar
    tick_values = [
        vmin,
        vmin + (vmax - vmin) * 0.33,
        vmin + (vmax - vmin) * 0.67,
        vmax
    ]
    # Format the tick labels
    tick_labels = [f"{v:.1f}" for v in tick_values] 

    if colorbar_location == 'top':
        cbar = plt.colorbar(img, ax=ax, orientation='horizontal', pad=cbartab, location='top', aspect=30)
        cbar.set_label(colorbar_label)
        # Override the ticks and tick labels
        cbar.set_ticks(tick_values)       # Set custom tick locations
        cbar.set_ticklabels(tick_labels)  # Set custom tick labels
        cbar.ax.xaxis.set_major_locator(FixedLocator(tick_values))  # Fix tick locations
        cbar.ax.xaxis.set_major_formatter(FixedFormatter(tick_labels))  # Fix tick labels
        # Suppress auto ticks completely
        cbar.ax.xaxis.set_minor_locator(FixedLocator([]))  # Ensure no minor ticks appear
    else:
        cbar = plt.colorbar(img, ax=ax, orientation='vertical', pad=cbartab, aspect=22)
        cbar.set_label(colorbar_label)
        # Override the ticks and tick labels
        cbar.set_ticks(tick_values)       # Set custom tick locations
        cbar.set_ticklabels(tick_labels)  # Set custom tick labels
        cbar.ax.yaxis.set_major_locator(FixedLocator(tick_values))  # Fix tick locations
        cbar.ax.yaxis.set_major_formatter(FixedFormatter(tick_labels))  # Fix tick labels
        # Suppress auto ticks completely
        cbar.ax.yaxis.set_minor_locator(FixedLocator([]))  # Ensure no minor ticks appear

    return ax