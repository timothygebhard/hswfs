"""
This module provides functions related to plotting.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from hswfs.utils import get_subaperture_centers


# -----------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# -----------------------------------------------------------------------------

def disable_ticks(
    ax: Any,
) -> None:
    """
    Disable the ticks and labels on the given matplotlib `ax`. This is
    similar to calling `ax.axis('off')`, except that the frame around
    the plot is preserved.

    Args:
        ax: A matplotlib axis.
    """

    ax.tick_params(axis='both', which='both', top=False, bottom=False,
                   left=False, right=False, labelbottom=False, labelleft=False)


def plot_shifts(
    ax: Any,
    relative_shifts: np.ndarray,
) -> None:
    """
    Create a plot of a wavefront sensor (with the `relative_shifts` for
    each subaperture) on the provided `ax`.

    Args:
        ax: A matplotlib axis.
        relative_shifts: A numpy array of shape `(N, N, 2)`, where `N`
            is the grid size, containing the relative shift (or offset)
            for each subaperture.
    """

    # Determine the grid size from the shifts array
    grid_size = relative_shifts.shape[0]

    # Compute the positions of the centers of the subapertures
    x, y = get_subaperture_centers(grid_size=grid_size)

    # Draw the grid of subapertures
    for z in np.linspace(-1 / np.sqrt(2), 1 / np.sqrt(2), grid_size + 1):
        ax.plot((z, z), (-1 / np.sqrt(2), 1 / np.sqrt(2)), color='black')
        ax.plot((-1 / np.sqrt(2), 1 / np.sqrt(2)), (z, z), color='black')

    # Plot the centers of the subapertures
    ax.plot(x.flatten(), y.flatten(), 'x', ms=4, color='C2', alpha=0.5)

    # Add a red circle indicating the unit disk
    ax.add_artist(plt.Circle((0, 0), 1, color='red', ls='--', fill=False))
 
    # Determine a shrinkage factor to map the relative shifts into the right
    # reference frame (i.e., scale to the size of a subaperture in the plot)
    factor = np.sqrt(2) / grid_size / 2

    # Plot the observed position for each subapertures, which differs from the
    # center of the subaperture by the given shift vector
    ax.plot(x.flatten() + factor * relative_shifts[:, :, 0].flatten(),
            y.flatten() + factor * relative_shifts[:, :, 1].flatten(),
            '.', color='C0')

    # Fix the x- and y-limits of the plot
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)

    # Fix the aspect ratio of the plot
    ax.set_aspect('equal')
