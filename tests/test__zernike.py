"""
Unit tests for functions in zernike.py
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

import sympy as sy

from hswfs.zernike import mn_to_j, j_to_mn, ZernikePolynomial


# -----------------------------------------------------------------------------
# TEST CASES
# -----------------------------------------------------------------------------

def test__mn_to_j() -> None:

    assert mn_to_j(0, 0) == 0
    assert mn_to_j(-1, 1) == 1
    assert mn_to_j(1, 1) == 2
    assert mn_to_j(-2, 2) == 3
    assert mn_to_j(0, 2) == 4
    assert mn_to_j(2, 2) == 5
    assert mn_to_j(-3, 3) == 6
    assert mn_to_j(-1, 3) == 7
    assert mn_to_j(1, 3) == 8
    assert mn_to_j(3, 3) == 9


def test__j_mn_to() -> None:

    assert j_to_mn(0) == (0, 0)
    assert j_to_mn(1) == (-1, 1)
    assert j_to_mn(2) == (1, 1)
    assert j_to_mn(3) == (-2, 2)
    assert j_to_mn(4) == (0, 2)
    assert j_to_mn(5) == (2, 2)
    assert j_to_mn(6) == (-3, 3)
    assert j_to_mn(7) == (-1, 3)
    assert j_to_mn(8) == (1, 3)
    assert j_to_mn(9) == (3, 3)


def test__radial_polynomial() -> None:

    rho = sy.symbols('rho')

    # -------------------------------------------------------------------------
    # Trivial case: n - m is odd
    # -------------------------------------------------------------------------

    zernike = ZernikePolynomial(m=0, n=1)
    expected = sy.sympify(0)
    assert sy.simplify(zernike.radial_part - expected) == 0

    # -------------------------------------------------------------------------
    # Interesting cases: n - m is even (test first six polynomials)
    # -------------------------------------------------------------------------

    zernike = ZernikePolynomial(m=0, n=0)
    expected = sy.sympify(1)
    assert sy.simplify(zernike.radial_part - expected) == 0

    zernike = ZernikePolynomial(m=1, n=1)
    expected = rho
    assert sy.simplify(zernike.radial_part - expected) == 0

    zernike = ZernikePolynomial(m=0, n=2)
    expected = 2 * rho**2 - 1
    assert sy.simplify(zernike.radial_part - expected) == 0

    zernike = ZernikePolynomial(m=2, n=2)
    expected = rho**2
    assert sy.simplify(zernike.radial_part - expected) == 0

    zernike = ZernikePolynomial(m=1, n=3)
    expected = 3 * rho**3 - 2 * rho
    assert sy.simplify(zernike.radial_part - expected) == 0

    zernike = ZernikePolynomial(m=3, n=3)
    expected = rho**3
    assert sy.simplify(zernike.radial_part - expected) == 0


def test__zernike_polynomial() -> None:

    rho, phi = sy.symbols('rho, phi')

    # -------------------------------------------------------------------------
    # Test first 10 Zernike polynomials
    # -------------------------------------------------------------------------

    zernike = ZernikePolynomial(m=0, n=0)
    expected = sy.sympify(1)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=-1, n=1)
    expected = 2 * rho * sy.sin(phi)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=1, n=1)
    expected = 2 * rho * sy.cos(phi)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=-2, n=2)
    expected = sy.sqrt(6) * rho**2 * sy.sin(2 * phi)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=0, n=2)
    expected = sy.sqrt(3) * (2 * rho**2 - 1)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=2, n=2)
    expected = sy.sqrt(6) * rho**2 * sy.cos(2 * phi)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=-3, n=3)
    expected = sy.sqrt(8) * rho**3 * sy.sin(3 * phi)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=-1, n=3)
    expected = sy.sqrt(8) * (3 * rho**3 - 2 * rho) * sy.sin(phi)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=1, n=3)
    expected = sy.sqrt(8) * (3 * rho**3 - 2 * rho) * sy.cos(phi)
    assert sy.simplify(zernike.polar - expected) == 0

    zernike = ZernikePolynomial(m=3, n=3)
    expected = sy.sqrt(8) * rho**3 * sy.cos(3 * phi)
    assert sy.simplify(zernike.polar - expected) == 0


def test__zernike_fourier_transform() -> None:

    k1, k2 = sy.symbols('k1, k2')

    zernike = ZernikePolynomial(m=0, n=0)
    expected = sy.besselj(2 * sy.pi * k1, 1) / (sy.pi * k1)
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=-1, n=1)
    expected = (-2 * sy.I * sy.besselj(2 * sy.pi * k1, 2) *
                sy.sin(k2) / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=1, n=1)
    expected = (-2 * sy.I * sy.besselj(2 * sy.pi * k1, 2) *
                sy.cos(k2) / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=-2, n=2)
    expected = (-sy.sqrt(6) * sy.besselj(2 * sy.pi * k1, 3) *
                sy.sin(2 * k2) / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=0, n=2)
    expected = (-sy.sqrt(3) * sy.besselj(2 * sy.pi * k1, 3) *
                1 / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=2, n=2)
    expected = (-sy.sqrt(6) * sy.besselj(2 * sy.pi * k1, 3) *
                sy.cos(2 * k2) / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=-3, n=3)
    expected = (sy.sqrt(8) * sy.I * sy.besselj(2 * sy.pi * k1, 4) *
                sy.sin(3 * k2) / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=-1, n=3)
    expected = (sy.sqrt(8) * sy.I * sy.besselj(2 * sy.pi * k1, 4) *
                sy.sin(1 * k2) / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=1, n=3)
    expected = (sy.sqrt(8) * sy.I * sy.besselj(2 * sy.pi * k1, 4) *
                sy.cos(1 * k2) / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0

    zernike = ZernikePolynomial(m=3, n=3)
    expected = (sy.sqrt(8) * sy.I * sy.besselj(2 * sy.pi * k1, 4) *
                sy.cos(3 * k2) / (sy.pi * k1))
    assert sy.simplify(zernike.fourier_transform - expected) == 0
