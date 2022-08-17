"""
Showcase the basic functions of the HSWFS class.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from typing import Any

import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from hswfs.plotting import plot_shifts, disable_ticks
from hswfs.sensor import HSWFS
from hswfs.shifts import generate_test_shifts
from hswfs.utils import crop_center, get_unit_disk_meshgrid
from hswfs.zernike import Wavefront, eval_cartesian


# -----------------------------------------------------------------------------
# MAIN CODE
# -----------------------------------------------------------------------------

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    # Set up virtual wavefront sensor, fit wavefront and compute PSF
    # -------------------------------------------------------------------------

    # Create a set of shifts for a HSWFS with a grid size of 8x8
    print("Generating shifts...", end=" ", flush=True)
    shifts = generate_test_shifts(test_case="defocus", grid_size=8)
    print("Done!", flush=True)

    # Set up a new HSWFS instance
    print("Setting up wavefront sensor...", end=" ", flush=True)
    sensor = HSWFS(relative_shifts=shifts)
    print("Done!", flush=True)

    # Fit the wavefront that corresponds to the shifts
    print("Fitting wavefront...", end=" ", flush=True)
    coefficients = sensor.fit_wavefront(n_zernike=9)
    wavefront = Wavefront(coefficients=coefficients)
    print("Done!", flush=True)

    # Compute the PSF that corresponds to that wavefront
    print("Computing PSF...", end=" ", flush=True)
    psf = sensor.get_psf(wavefront=wavefront, resolution=64)
    print("Done!", flush=True)

    # -------------------------------------------------------------------------
    # Create a plot of the results
    # -------------------------------------------------------------------------

    print("Plotting results...", end=" ", flush=True)

    # Set up a grid of plots for the results
    gs = gridspec.GridSpec(2, 2, width_ratios=[2.25, 1], height_ratios=[1, 1])
    ax_sensor: Any = plt.subplot(gs[:, :-1])
    ax_wavefront: Any = plt.subplot(gs[:-1, -1])
    ax_psf: Any = plt.subplot(gs[-1, -1])

    # Plot the shifts on the wavefront sensor
    plot_shifts(ax=ax_sensor, relative_shifts=shifts)
    ax_sensor.set_title("Shifts on wavefront sensor")
    ax_sensor.set_aspect("equal")

    # Evaluate the wavefront on a grid and plot it
    x_0, y_0 = get_unit_disk_meshgrid(resolution=512)
    wf_grid = eval_cartesian(wavefront.cartesian, x_0=x_0, y_0=y_0)
    limit = 1.1 * np.nanmax(np.abs(wf_grid))
    ax_wavefront.imshow(
        wf_grid,
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-limit,
        vmax=limit,
    )
    ax_wavefront.set_title("Wavefront")

    # Plot the PSF (cropped to center)
    ax_psf.imshow(
        crop_center(psf, size=(33, 33)),
        interpolation="nearest",
        norm=colors.LogNorm(vmin=0.001, vmax=1),
    )
    ax_psf.set_title("PSF (log-scale)")

    # Disable ticks on all subplots
    disable_ticks(ax_sensor)
    disable_ticks(ax_wavefront)
    disable_ticks(ax_psf)

    # Save the result
    plt.tight_layout()
    plt.savefig("example_result.png", pad=0, dpi=300, bbox_inches="tight")

    print("Done!", flush=True)
