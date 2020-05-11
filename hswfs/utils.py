"""
This module provides various utility functions.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from typing import Tuple

import numpy as np


# -----------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# -----------------------------------------------------------------------------

def get_unit_disk_meshgrid(
    resolution: int,
) -> Tuple[np.array, np.array]:
    """
    Get a (Cartesian) mesh grid of positions on the unit disk, that is,
    all positions with with a Euclidean distance <= 1 from (0, 0).
 
    Args:
        resolution: An integer specifying the size of the mesh grid,
            that is, the number of points in each dimensions.

    Returns:
        A mesh grid consisting of the tuple `x_0`, `y_0`, which are each
        numpy arrays of shape `(resolution, resolution)`. For positions
        that are on the unit disk, they contain the coordinates of the
        position; otherwise, they contain `np.nan`.
    """

    # Create a meshgrid of (Cartesian) positions: [-1, 1] x [-1, 1]
    x_0, y_0 = np.meshgrid(np.linspace(-1, 1, resolution),
                           np.linspace(-1, 1, resolution))

    # Create a mask for the unit disk (only select position with radius <= 1)
    unit_disk_mask = np.sqrt(x_0**2 + y_0**2) <= 1

    # Mask out all the position that are not on the unit disk
    x_0[~unit_disk_mask] = np.nan
    y_0[~unit_disk_mask] = np.nan

    return x_0, y_0