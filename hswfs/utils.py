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

def crop_center(array: np.ndarray, size: Tuple[int, ...]) -> np.ndarray:
    """
    Crop an n-dimensional array to the given size around its center.

    Args:
        array: The numpy array to be cropped.
        size: A tuple containing the size of the cropped array. To not
            crop along a specific axis, you can specify the size of
            that axis as -1.
    Returns:
        The input array, cropped to the desired size around its center.
    """

    # Ensure that the the array shape and the size variable match
    if array.ndim != len(size):
        raise RuntimeError(
            "Length of size must match number of dimensions of array!"
        )

    # Loop over the the axes of the array to create slices
    slices = list()
    for old_len, new_len in zip(array.shape, size):

        # Compute start and end position for axis
        start = old_len // 2 - new_len // 2 if new_len != -1 else None
        end = start + new_len if start is not None else None

        # Create a slice object for axis
        slices.append(slice(start, end))

    return np.asarray(array[tuple(slices)])


def get_subaperture_centers(
    grid_size: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the positions of the centers of the subapertures of the
    sensor. This assumes a simple geometry, where the sensor is taken
    to be the largest square that can fit inside the unit circle, and
    consists of a grid of `grid_size` x `grid_size` subapertures.

    Args:
        grid_size:  An integer specifying the size of the (quadratic)
            grid of subapertures in the HSWFS sensor.

    Returns:
        A mesh grid, consisting of two numpy arrays which specify the
        `x` and `y` positions of the centers of the subapertures.
    """

    x = np.linspace((1 / grid_size - 1), (1 - 1 / grid_size), grid_size)
    x = 1 / np.sqrt(2) * np.repeat(x.reshape(1, -1), grid_size, axis=0)
    y = np.linspace((1 - 1 / grid_size), (1 / grid_size - 1), grid_size)
    y = 1 / np.sqrt(2) * np.repeat(y.reshape(-1, 1), grid_size, axis=1)

    return x, y


def get_unit_disk_meshgrid(
    resolution: int,
) -> Tuple[np.ndarray, np.ndarray]:
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
    x_0, y_0 = np.meshgrid(
        np.linspace(-1, 1, resolution), np.linspace(-1, 1, resolution)
    )

    # Create a mask for the unit disk (only select position with radius <= 1)
    unit_disk_mask = np.sqrt(x_0**2 + y_0**2) <= 1

    # Mask out all the position that are not on the unit disk
    x_0[~unit_disk_mask] = np.nan
    y_0[~unit_disk_mask] = np.nan

    return x_0, y_0
