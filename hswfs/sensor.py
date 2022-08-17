"""
This module provides the main class that is used to simulates a
Hartmann-Shack wavefront sensor (HSWFS).
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

import numpy as np
import scipy as sp

from hswfs.utils import get_unit_disk_meshgrid
from hswfs.zernike import (
    derive,
    eval_cartesian,
    j_to_mn,
    Wavefront,
    ZernikePolynomial,
)
from hswfs.fast_zernike import zernike_derivative_cartesian


# -----------------------------------------------------------------------------
# CLASS DEFINITIONS
# -----------------------------------------------------------------------------

class HSWFS:
    """
    A virtual Hartmann-Shack wavefront sensor.

    Args:
        relative_shifts: A numpy array of shape `(N, N, 2)`, where `N`
            is the grid size, containing the relative shift (or offset)
            for each subaperture. "Relative" means that the shift values
            are in `[-1, 1]`.
        aperture_size: An integer specifying the size of a subaperture
            in pixels (subapertures are assumed to be quadratic).
    """

    def __init__(
        self,
        relative_shifts: np.ndarray,
        aperture_size: int = 32,
    ) -> None:

        # Store constructor arguments
        self.relative_shifts = relative_shifts
        self.aperture_size = aperture_size

        # Store derived properties
        self.grid_size = self.relative_shifts.shape[0]

    def fit_wavefront(
        self,
        n_zernike: int = 9,
    ) -> np.ndarray:
        """
        Perform a modal (i.a. Zernike polynomial-based) least squares
        fit to find the wavefront that matches the sensor data. The
        result is a vector of Zernike coefficients.

        For more details (also about, for example, the notation), please
        see [Cubalchini_1979]_. For an alternative discussion of the
        same topic, see also, for example, section 4.3.2 of [Dai_2007]_.

        .. note::
            Note that the code in this package uses a slightly different
            notation for indexing the Zernike polynomials: when using
            only a single index :math:`j`, we start counting them at
            :math:`j = 0` (the constant polynomial which corresponds to
            :math:`Z^0_0` when using double indices :math:`m` and
            :math:`n`). [Cubalchini_1979]_, on the other hand, starts
            counting the Zernike polynomials at 1, which is clearly
            unpythonic :)

        .. note::
            For compatibility reasons, the coefficient vector that is
            return by the function will also contain a value for the
            first Zernike polynomial :math:`Z^0_0`, which is always set
            to be zero. the coefficient vector will, therefore, have
            `n_zernike + 1` entries.

        .. warning::
            Fitting the wavefront requires the Cartesian derivatives of
            the Zernike polynomials. Computing these "on the fly" is
            relatively slow, which is why there exists the
            :py:mod:`hswfs.fast_zernike` module, which contains
            pre-computed versions of these derivatives up to
            :math:`j = 135`. If `n_zernike` is chosen to be larger than
            this value, fitting the wavefront may take a relatively
            long time. However, for most practical purposes, using this
            many Zernike polynomials will probably be unnecessary.

        Args:
            n_zernike: The number of Zernike polynomials to be used in
                the fit. Note that :math:`Z^0_0` (the constant term) is
                ignored in the fit.
                The index `j_max` of the Zernike polynomial with the
                highest order will, therefore, be `n_zernike + 1`.

        Returns:
            A numpy array of length `n_zernike + 1`, where the
            :math:`j`-th entry is the coefficient that corresponds to
            the Zernike polynomial :math:`Z_j`.
        """

        # ---------------------------------------------------------------------
        # Compute the P vector
        # ---------------------------------------------------------------------

        # Shortcut for the total number of apertures in the grid
        n_ap = self.grid_size**2

        # Create p-vector from measured shifts:
        #   (x_shift_1, x_shift_2, ..., x_shift_N, y_shift_1,  ..., y_shift_n)
        p = np.concatenate(
            (
                self.relative_shifts[:, :, 0].reshape(n_ap),
                self.relative_shifts[:, :, 1].reshape(n_ap),
            )
        )

        # ---------------------------------------------------------------------
        # Compute the D matrix of derivatives of Zernike polynomials
        # ---------------------------------------------------------------------

        # Initialize D.
        # According to eq. (14) in [Cubalchini_1979], D is a matrix of shape:
        #   (n_zernike, 2 * n_subapertures),
        # with entries defined as:
        #   D_ab = \frac{\partial Z_a}{\partial x} (x, y)_b
        #   if 0 < b <= n_subapertures,
        # and
        #   D_ab = \frac{\partial Z_a}{\partial y} (x, y)_b'
        #   if n_subapertures < b <= 2 *n_apertures,
        # where b' = b - n_subapertures, and (x, y)_b denotes the relative
        # position of the center of the b-th subaperture assuming the entire
        # sensor grid is placed on the unit disk.
        d = np.full((n_zernike, 2 * n_ap), np.nan)

        # Compute evaluation position (x, y), that is, the relative positions
        # of the centers of the subapertures
        x_0 = (
            1
            / np.sqrt(2)
            * np.linspace(
                (1 / self.grid_size - 1),
                (1 - 1 / self.grid_size),
                self.grid_size,
            ).reshape(1, -1)
        )
        x_0 = np.repeat(x_0, self.grid_size, axis=0)
        y_0 = (
            1
            / np.sqrt(2)
            * np.linspace(
                (1 - 1 / self.grid_size),
                (1 / self.grid_size - 1),
                self.grid_size,
            ).reshape(-1, 1)
        )
        y_0 = np.repeat(y_0, self.grid_size, axis=1)

        # We compute the entries of D row by row, because rows correspond to
        # Zernike polynomials
        for row_idx, j in enumerate(range(1, n_zernike + 1)):

            # Map single-index j to double-indices m, n
            m, n = j_to_mn(j)

            # Compute derivatives in x- and y-direction. For "small" values of
            # j, we can use the pre-computed derivatives from the fast_zernike
            # module. For even higher orders, the derivatives first need to be
            # computed "on demand", which is a lot slower.
            if j <= 135:
                x_derivatives = zernike_derivative_cartesian(
                    m, n, x_0, y_0, "x"
                )
                y_derivatives = zernike_derivative_cartesian(
                    m, n, x_0, y_0, "y"
                )
            else:
                zernike_polynomial = ZernikePolynomial(m=m, n=n).cartesian
                x_derivatives = eval_cartesian(
                    derive(zernike_polynomial, "x"), x_0, y_0
                )
                y_derivatives = eval_cartesian(
                    derive(zernike_polynomial, "y"), x_0, y_0
                )

            # Store derivatives in D matrix
            d[row_idx] = np.concatenate(
                (x_derivatives.flatten(), y_derivatives.flatten())
            )

        # ---------------------------------------------------------------------
        # Find the Zernike coefficients by solving a linear equation system
        # ---------------------------------------------------------------------

        # Define the matrix E as per eq. (15) in [Cubalchini_1979]
        e = d @ d.transpose()

        # Finally, we can compute the coefficient vector of the wavefront.
        # According to eq. (16) in [Cubalchini_1979], we have:
        #   A = E^{-1} D P
        # Instead of computing the right-hand side directly by inverting E,
        # which is often ill-conditioned, resulting in major numerical
        # instabilities and incorrect fits, we transform this to:
        #   EA = DP,
        # and use a least squares solver to solve the  equation system for A.
        a = sp.linalg.lstsq(a=e, b=d @ p)[0]

        # Add an additional 0 in the first position for for Z^_0
        a = np.insert(a, 0, 0)

        return np.asarray(a)

    @staticmethod
    def get_psf(
        wavefront: Wavefront,
        resolution: int = 256,
    ) -> np.ndarray:
        r"""
        Compute the point spread function (PSF) that corresponds to the
        given `wavefront`. This needs to happen in discretized form,
        because in the general case, the integral that connects the PSF
        and the wavefront is not analytically tractable.

        The mathematics in this function --- using the FFT to compute
        the PSF --- is based on section (14.15), particularly on
        eq. (14.188), in [Malacara_2007]_, as well as the exemplary
        MATLAB implementation on appendix 8.C of [Dai_2007]_.

        For small abberations, there exists an analytical approximation
        based on the Fourier transformations of the Zernike polynomials.
        It can be found in eq. (8.28) in [Dai_2007]_, where it is also
        formally derived in appendix 8.B.

        .. warning::
            Please note that the point spread function will, in general,
            also depend on the wavelength :math:`\lambda`. This
            dependency is currently not correctly modeled by the
            `get_psf()` function; instead, :math:`\lambda = 1` is used.

        Args:
            wavefront:
            resolution: The size of the grid on which the wavefront is
                computed (and thus also the PSF).

        Returns:
            The PSF that corresponds to the wavefront that was fitted to
            the sensor data.
        """

        # Get a grid of the unit disk
        x_0, y_0 = get_unit_disk_meshgrid(resolution=resolution)

        # Compute the wavefront on a grid of the given resolution, and cast
        # np.nan to 0, because the FFT cannot deal with NaN
        wf_grid = eval_cartesian(
            expression=wavefront.cartesian, x_0=x_0, y_0=y_0
        )
        wf_grid = np.nan_to_num(wf_grid)

        # Compute the pupil function. In our simple case, this is simply the
        # unit disk, i.e., it is 1 where the radius is <= 1, and 0 elsewhere
        pupil = np.logical_not(np.logical_or(np.isnan(x_0), np.isnan(y_0)))
        pupil = pupil.astype(float)

        # Define a wavelength lambda
        # FIXME: Proper treatment of wavelength is not yet implemented!
        lambda_ = 1

        # Compute the corresponding point spread function (PSF)
        psf = pupil * np.exp(2 * np.pi * 1j / lambda_ * wf_grid)
        psf = np.fft.fft2(psf)
        psf = np.abs(psf) ** 2
        psf = np.fft.fftshift(psf) / np.max(psf)

        return np.asarray(psf)
