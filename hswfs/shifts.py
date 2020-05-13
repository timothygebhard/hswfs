"""
This module provides functions to generate shift grids to be used with
the HSWFS class.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from typing import Optional

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

import numpy as np


# -----------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# -----------------------------------------------------------------------------

def generate_random_shifts(
    grid_size: int = 16,
    smooth_std: Optional[float] = None,
    random_seed: int = 42,
) -> np.array:
    r"""
    Create a 2D grid of random shifts (or offsets, or delta) for the
    subapertures of the Hartmann-Shack wavefront sensor.

    By default, these shifts will independently of each other follow a
    uniform distribution for both :math:`\Delta x` and :math:`\Delta y`.
    This is almost certainly an oversimplification; however, I do not
    know the "true" distribution that these shifts should follow for
    real atmospheric perturbations.

    To make the results slightly more realistic, there is the option
    to smooth (i.e., convolve) the vector field of shifts with an 2D
    Gaussian kernel, which will introduce some local correlations
    between neighboring apertures.

    Args:
        grid_size: An integer specifying the size of the (quadratic)
            grid of subapertures in the HSWFS sensor. Default: 16.
        smooth_std: The standard deviation of the 2D Gaussian kernel
            that is used to smooth the results. If None, no smoothing
            is applied (default).
        random_seed: Seed to be used for the random number generator.

    Returns:
        A numpy array of shape `(grid_size, grid_size, 2)` which, for
        each subaperture in the wavefront sensor, contains the a 2D
        shift vector `(delta_x, delta_y)` specifying the *relative*
        offsets.
        "Relative" means that the values are in `[-1, 1]` and may still
        need to be multiplied with the (pixel) size of the subapertures.
    """

    # Create a new RNG object: this ensures reproducibility independently of
    # the state of the global numpy.random random number generator.
    rng = np.random.RandomState(seed=random_seed)

    # Start by simply generating a (grid_size, grid_size, 2)-sized grid of
    # shifts that are sampled independently and uniformly from [-1, 1].
    shifts = rng.uniform(-1, 1, (grid_size, grid_size, 2))

    # If desired, smooth the vector field of shifts with a 2D Gaussian kernel
    # to introduce some local correlations between neighboring grid cells.
    if smooth_std is not None:
        kernel = Gaussian2DKernel(x_stddev=smooth_std, y_stddev=smooth_std)
        shifts[:, :, 0] = convolve(shifts[:, :, 0], kernel)
        shifts[:, :, 1] = convolve(shifts[:, :, 1], kernel)

    # Make sure that all values are still in [-1, 1]
    shifts = np.clip(shifts, a_min=-1, a_max=1)

    return shifts


def generate_test_shifts(
    test_case: str,
    grid_size: int = 16,
) -> np.ndarray:
    """
    Create a special 2D grid of shifts (or offsets, or delta) for the
    subapertures of the Hartmann-Shack wavefront sensor, for which the
    the shape of the corresponding wavefront is known. This is useful
    for testing.

    Currently, there are three such special cases:

    1. **x_shift:**
       The shifts in all subapertures are equal and only in x-direction.
       In this case, the resulting wavefront can be described using only
       the Zernike polynomial :math:`Z^{-1}_{1}`; that is, if we try to
       fit the wavefront corresponding to this case, the coefficients of
       all Zernike polynomials should be 0, except for :math:`j=1`.
    2. **y_shift:**
       The shifts in all subapertures are equal and only in y-direction.
       In this case, the resulting wavefront can be described using only
       the Zernike polynomial :math:`Z^{1}_{1}`; that is, if we try to
       fit the wavefront corresponding to this case, the coefficients of
       all Zernike polynomials should be 0, except for :math:`j=2`.
    3. **defocus:**
       The shifts in all subapertures point in the direction that is
       given by the vector from the center of the sensor to the
       respective aperture, and the value of the shift is proportional
       to the square of the distance of the subaperture from the center
       of the sensor.
       In this case, the resulting wavefront can be described using only
       the Zernike polynomial :math:`Z^{0}_{2}`; that is, if we try to
       fit the wavefront corresponding to this case, the coefficients of
       all Zernike polynomials should be 0, except for :math:`j=4`.

    Args:
        test_case: A string containing the name of one of the three
            cases described above, that is, one of the following:
            `"x_shift"`, `"y_shift"` or `"defocus"`.
        grid_size: An integer specifying the size of the (quadratic)
            grid of subapertures in the HSWFS sensor. Default: 16.

    Returns:
        A numpy array of shape `(grid_size, grid_size, 2)` which, for
        each subaperture in the wavefront sensor, contains the a 2D
        shift vector `(delta_x, delta_y)` specifying the *relative*
        offsets for the chosen test case.
    """

    # Initialize shifts array as all zeroes
    shifts = np.zeros((grid_size, grid_size, 2))

    #
    if test_case == 'x_shift':
        shifts[:, :, 0] = 0.5

    elif test_case == 'y_shift':
        shifts[:, :, 1] = 0.5

    elif test_case == 'defocus':

        # Compute the sign of the x- and y-shift for each subapertures
        col, row = np.meshgrid(np.arange(grid_size), np.arange(grid_size))
        x_sign = np.sign((col + 0.5) - grid_size / 2)
        y_sign = np.sign(grid_size / 2 - (row + 0.5))

        # Compute the x- and y-shifts
        x_shift = np.round((col + 0.5) - grid_size / 2 + 0.001 * x_sign) ** 2
        y_shift = np.round(grid_size / 2 - (row + 0.5) + 0.001 * y_sign) ** 2

        # Combine the x- and y-shifts into a single 3D array
        shifts = np.array([x_sign * x_shift, y_sign * y_shift])
        shifts = np.moveaxis(shifts, 0, -1)

        # Normalize such that max(shifts) == 1/2
        shifts /= np.max(shifts) * 2

    # Raise an error for all other values of
    else:
        raise ValueError('test_case must be one of the following: '
                         '"x_shift", "y_shift", "defocus"!')

    return shifts
