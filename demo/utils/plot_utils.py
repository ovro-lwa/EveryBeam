# Copyright (C) 2021 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import glob
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import matplotlib

matplotlib.use("Agg")

"""
Plotting utilities that can be used to visualize parts of the OSKAR demos.
"""


def plot_fits_panel(ax, image, title, cmap):
    """
    Plot a single fits panel

    Parameters
    ----------
    ax : matplotlib.axis
        Axis to which image should be added
    image : np.ndarray
        Array containing the image data.
    title : str
        Title for plot
    cmap : matplotlib.colormap
        Colormap to use
    """

    im = ax.imshow(np.squeeze(image), cmap=cmap)
    plt.colorbar(im, format="%.2e")
    plt.tick_params(
        labelcolor="none", top="off", bottom="off", left="off", right="off"
    )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    ax.set_title(title)
    ax.axis("equal")


def make_fits_plots(title, glob_pattern, out_basename):
    """
    Generates a plot of the fits image with four panels
    (either all Stokes components, or all polarizations)

    Parameters
    ----------
    title : str
        Title for plot.
    glob_pattern : str
        Globbing pattern for (fits) filenames that are to be plotted.
    out_basename : str
        Basename for output images.
    """

    # Load FITS images matching the glob pattern.
    files = glob.glob(glob_pattern)
    files.sort()  # Must be sorted!
    images = []
    for file in files:
        images.append(fits.getdata(file))

    # Set titles and colour maps to use.
    cmap = ""
    titles = []
    if "_AUTO_POWER" in glob_pattern:
        cmap = "Blues_r"
        titles = ["Stokes I", "Stokes Q", "Stokes U", "Stokes V"]
    elif "_AMP" in glob_pattern:
        cmap = "CMRmap"
        titles = ["XX", "XY", "YX", "YY"]

    # Sub-plots, one for each Stokes parameter.
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(221, frameon=False)
    plot_fits_panel(ax, images[0], titles[0], cmap)
    ax = fig.add_subplot(222, frameon=False)
    plot_fits_panel(ax, images[1], titles[1], cmap)
    ax = fig.add_subplot(223, frameon=False)
    plot_fits_panel(ax, images[2], titles[2], cmap)
    ax = fig.add_subplot(224, frameon=False)
    plot_fits_panel(ax, images[3], titles[3], cmap)

    # Add main title.
    title = title.replace("_", " ")  # Replace underscores with spaces.
    fig.suptitle(title)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)

    # Save and close.
    plt.savefig(f"{out_basename}.png")
    plt.close("all")
