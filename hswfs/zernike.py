"""
This module provides methods related to Zernike polynomials and their
derivatives.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from copy import deepcopy
from typing import Callable, Dict, Optional, Sequence, Tuple, Union

import numpy as np
import sympy as sy


# -----------------------------------------------------------------------------
# AUXILIARY FUNCTIONS
# -----------------------------------------------------------------------------

def mn_to_j(
    m: int,
    n: int,
) -> int:
    r"""
    Map the indices :math:`m, n` from the double-indexing scheme to
    the corresponding index :math:`j` of the single-indexing scheme.
    Basically, we are just counting the Zernike polynomials in a
    well-defined fashion.

    Mathematically, the mapping is given by:

    .. math::
        j = \frac{n \cdot (n + 2) + m}{2}

    Args:
        m: Index :math:`m` of :math:`Z^m_n`.
        n: Index :math:`m` of :math:`Z^m_n`.

    Returns:
        The single index :math:`j` which corresponds to :math:`m, n`.
    """

    return int((n * (n + 2) + m) / 2)


def j_to_mn(
    j: int,
) -> Tuple[int, int]:
    r"""
    Map the index :math:`j` of the single-indexing scheme to the
    corresponding indices :math:`m, n` from the double-indexing scheme.

    Mathematically, the mapping is given by:

    .. math::
        n = \left\lceil (-3 + \sqrt{9 + 8 j})\, /\, 2 \right\rceil
        \quad \text{and} \quad
        m = 2 j - n \cdot (n + 2)

    Args:
        j: Index :math:`j` of :math:`Z_j`.

    Returns:
        The pair of indices :math:`m, n` which correspond to :math:`j`.
    """

    n = int(np.ceil((-3 + np.sqrt(9 + 8 * j)) / 2))
    m = int(2 * j - n * (n + 2))

    return m, n


def polar_to_cartesian(
    expression: sy.Expr,
) -> sy.Expr:
    r"""
    Convert a sympy expression (`sy.Expr`) from polar coordinates to
    Cartesian coordinates by substituting :math:`\rho` and :math:`\phi`
    by the appropriate functions of :math:`x` and :math:`y`.

    Mathematically, the coordinate transformation is given by the
    following substitutions:

    .. math::
        \rho = \sqrt{x^2 + y^2}
        \quad \text{and} \quad
        \phi = \arctan\left( \frac{y}{x} \right)

    Args:
        expression: A sympy expression in polar coordinates, i.e., it
            must contain two free symbols named "rho" and "phi".

    Returns:
        The original `expression`, converted to Cartesian coordinates.
    """

    # Define symbols for polar and cartesian coordinates
    rho, phi = sy.symbols('rho'), sy.symbols('phi')
    x, y = sy.symbols('x'), sy.symbols('y')

    # Define coordinate transformation between polar and cartesian
    substitute_rho = sy.sqrt(x**2 + y**2)
    substitute_phi = sy.atan2(y, x)

    # Substitute rho and phi by the respective functions of x and y
    result = deepcopy(expression)
    result = result.subs(rho, substitute_rho)
    result = result.subs(phi, substitute_phi)

    return result


def derive(
    expression: sy.Expr,
    wrt: Union[str, sy.Symbol]
) -> sy.Expr:
    """
    Compute the derivative of `expression` with respect to `wrt`.

    Args:
        expression: A sympy expression.
        wrt: The variable with respect to which to take the derivative.
            Can either be a `sy.Symbol` or a string containing the name
            of the variable.

    Returns:
        The derivative of `expression` with respect to `wrt` as a
        `sy.Expr`.
    """

    # If wrt is a sympy symbol, and is part of the expression's free symbols,
    # we can directly take the derivative of the expression w.r.t. wrt
    if isinstance(wrt, sy.Symbol) and (wrt in expression.free_symbols):
        return sy.diff(expression, wrt)

    # If wrt is a string, we check if there is a symbol in the free symbols of
    # the expression whose name matches wrt, and then we take the derivative
    # w.r.t. to this symbol
    elif isinstance(wrt, str):
        for symbol in expression.free_symbols:
            if wrt == symbol.name:  # type: ignore
                return sy.diff(expression, symbol)

    # In every other case, the derivative is simply 0
    return sy.sympify(0)


def is_cartesian(
    expression: sy.Expr
) -> bool:
    """
    Check if a given `expression` is in Cartesian coordinates, that is,
    if the names of its free symbols are a subset of `{"x", "y"}`.

    Args:
        expression: A sympy expression.

    Returns:
        True if `expression` is in Cartesian coordinates; else False.
    """

    symbols = {_.name for _ in expression.free_symbols}  # type: ignore
    return symbols.issubset({'x', 'y'})


def is_polar(
    expression: sy.Expr
) -> bool:
    """
    Check if a given `expression` is in polar coordinates, that is, if
    the names of its free symbols are a subset of `{"rho", "phi"}`.

    Args:
        expression: A sympy expression.

    Returns:
        True if `expression` is in polar coordinates; else False.
    """

    symbols = {_.name for _ in expression.free_symbols}  # type: ignore
    return symbols.issubset({'rho', 'phi'})


def eval_cartesian(
    expression: sy.Expr,
    x_0: Union[float, np.ndarray],
    y_0: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    """
    Evaluate an expression that is in in Cartesian coordinates, either
    at a single position or on a grid of positions.

    Args:
        expression: A sympy expression.
        x_0: The value(s) of :math:`x` at which to evaluate the given
            `expression`. This can either be a single float, or an array
            of arbitrary size (its shape, however, must match `y_0`).
        y_0: The value(s) of :math:`y` at which to evaluate the given
            `expression`. This can either be a single float, or an array
            of arbitrary size (its shape, however, must match `x_0`).

    Returns:
        The value of `expression` at the given position(s). The type
        and shape of the output matches the one of the input: for `x_0`,
        `y_0` as floats, a float is returned; for numpy array inputs, a
        numpy array is returned.
    """

    # Make sure that expression is a function of Cartesian coordinates
    assert is_cartesian(expression), \
        '"expression" is not in Cartesian coordinates!'

    # Make sure that x_0 and y_0 have compatible shapes
    assert ((isinstance(x_0, float) and isinstance(y_0, float)) or
            (isinstance(x_0, np.ndarray) and isinstance(y_0, np.ndarray) and
             x_0.shape == y_0.shape)), \
        '"x_0" and "y_0" must be either both float, or both numpy array ' \
        'with the same shape!'

    # If the expression is not constant, we can use sympy.lambdify() to
    # generate a numpy version of the expression, which can be used to
    # evaluate the function efficiently:
    if not expression.is_constant():

        numpy_func: Callable[..., Union[float, np.ndarray]] = \
            sy.utilities.lambdify(args=sy.symbols('x, y'),
                                  expr=expression,
                                  modules='numpy')

    # Otherwise, that is, if the expression is constant, we need to define
    # the evaluation function manually because the result of sympy.lambdify()
    # does not behave as desired (it does not vectorize properly).
    else:

        # The multiplication with _ / _ makes sure that everything that is NaN
        # in the input also is NaN in the output; non-NaN values are unchanged
        def numpy_func(_: float, __: float) -> float:
            return float(expression) * _ / _ * __ / __
        numpy_func = np.vectorize(numpy_func)

    return numpy_func(x_0, y_0)


# -----------------------------------------------------------------------------
# CLASSES
# -----------------------------------------------------------------------------

class ZernikePolynomial:
    """
    Implements the Zernike polynomial :math:`Z^m_n` (in double-index
    notation), or :math:`Z_j` (in single-index notation).
    """

    def __init__(self,
                 m: Optional[int] = None,
                 n: Optional[int] = None,
                 j: Optional[int] = None):

        # Make sure that we have received *either* (m, n) *or* j
        error_msg = 'ZernikePolynomial must be instantiated either with ' \
                    'double indices (m, n) *or* a single index j!'
        if j is not None:
            assert m is None, error_msg
            assert n is None, error_msg
            self.j = j
            self.m, self.n = j_to_mn(self.j)
        else:
            assert m is not None, error_msg
            assert n is not None, error_msg
            self.m, self.n = m, n
            self.j = mn_to_j(self.m, self.n)

        # Run basic sanity checks on inputs
        assert (-self.n <= self.m <= self.n), \
            'Zernike polynomials are only defined for -n <= m <= n!'
        assert self.j >= 0, \
            'Zernike polynomials are only defined for j >= 0!'

    def __repr__(self) -> str:
        return f'Z^{self.m}_{self.n}'

    @property
    def radial_part(self) -> sy.Expr:
        r"""
        The radial polynomial :math:`R^m_n`, which is given by:

        .. math::
            R^m_n ( \rho ) = \sum_{k=0}^{\frac{n-m}{2}} (-1)^{k} \,
                {{n - k} \choose {k}} \,
                {{n - 2k} \choose {\frac{n - m}{2} - k}} \,
                \rho^{n - 2k}

        Returns:
            The radial polynomial :math:`R^m_n` as a `sympy.Expr`.
        """

        # Define a symbol for the radius (rho)
        rho = sy.Symbol('rho')

        # If n - m is odd, the radial polynomial is simply 0
        if (self.n - self.m) % 2 == 1:
            return sy.sympify(0)

        # Otherwise, things are a little more complicated
        else:
            return sum(
                sy.Pow(-1, k) * sy.binomial(int(self.n - k), int(k)) *
                sy.binomial(
                    int(self.n - 2 * k), int((self.n - self.m) / 2 - k)
                ) *
                sy.Pow(rho, self.n - 2 * k)
                for k in range(0, int((self.n - self.m) / 2) + 1)
            )

    @property
    def azimuthal_part(self) -> sy.Expr:
        r"""
        The azimuthal component of :math:`Z^m_n`, which is given by:

        .. math::
            \Phi_m ( \phi ) =
                \begin{cases}
                    \cos( m \, \phi ) & \text{for}\ m \geq 0 \\
                    \sin( m \, \phi ) & \text{for}\ m < 0
                \end{cases}

        Returns:
            The azimuthal part of :math:`Z^m_n` as a `sympy.Expr`.
        """

        # Define a symbol for the azimuthal angle (phi)
        phi = sy.Symbol('phi')

        # Return the azimuthal part, which depends only on the value of m
        if self.m > 0:
            return sy.cos(self.m * phi)
        elif self.m < 0:
            return sy.sin(-self.m * phi)
        else:
            return sy.sympify(1)

    @property
    def normalization(self) -> sy.Expr:
        r"""
        The normalization factor of :math:`Z^m_n`.

        Zernike polynomials are normalized such that:

        .. math::
            \int_0^1 d\rho \int_0^{2\pi}d\phi \
            Z^2(\rho, \phi) \, \rho = \pi

        .. note::
            Note that this choice of normalization is not universal,
            and that some authors choose different normalizations;
            for example, they normalize the above integral to 1
            instead of :math:`\pi`.

        Returns:
            The normalization factor of :math:`Z^m_n` as a `sy.Expr`.
        """

        if self.m == 0:
            return sy.sqrt(self.n + 1)
        else:
            return sy.sqrt(self.n + 1) * sy.sqrt(2)

    @property
    def polar(self) -> sy.Expr:
        r"""
        :math:`Z^m_n(\rho, \phi)`, that is, the full Zernike polynomial
        in polar coordinates :math:`\rho, \phi`.

        Returns:
            The Zernike polynomial :math:`Z^m_n(\rho, \phi)` as
            a `sy.Expr`.
        """
        return self.normalization * self.radial_part * self.azimuthal_part

    @property
    def cartesian(self) -> sy.Expr:
        r"""
        :math:`Z^m_n(x, y)`, that is, the full Zernike polynomial in
        Cartesian coordinates :math:`x, y`.

        Returns:
            The Zernike polynomial :math:`Z^m_n(x, y)` as a `sy.Expr`.
        """
        return polar_to_cartesian(self.polar)

    @property
    def fourier_transform(self) -> sy.Expr:
        r"""
        The 2D Fourier transform of :math:`Z^m_n`.

        This function essentially implements eq. (7) of [Tatulli_2013]_.

        .. note::
            Compared to [Tatulli_2013]_, we are using a slightly
            different notation. Instead of indexing the dimensions of
            the Fourier space (or :math:`k`-space) using :math:`\kappa`
            and :math:`\alpha`, we use :math:`k_1` and :math:`k_2`.

        Returns:
            The 2D Fourier transform of :math:`Z^m_n`, that is,
            :math:`\mathcal{F}\lbrace Z^m_n \rbrace (k_1, k_2)`,
            as a `sy.Expr`.
        """

        # Define symbols for k1 and k2
        k1 = sy.Symbol('k1')
        k2 = sy.Symbol('k2')

        # Define the first factor, which only depends on n
        factor_1 = (sy.Pow(-1, self.n) * sy.sqrt(self.n + 1) / (sy.pi * k1) *
                    sy.besselj(2 * sy.pi * k1, self.n + 1))

        # Define the second factor that also depends on m
        if self.m == 0:
            factor_2 = sy.Pow(-1, self.n / 2)
        elif self.m > 0:
            factor_2 = (sy.sqrt(2) * sy.Pow(-1, (self.n - self.m) / 2) *
                        sy.Pow(sy.I, self.m) * sy.cos(self.m * k2))
        else:
            factor_2 = (sy.sqrt(2) * sy.Pow(-1, (self.n + self.m) / 2) *
                        sy.Pow(sy.I, -self.m) * sy.sin(-self.m * k2))

        return sy.nsimplify(sy.simplify(factor_1 * factor_2))


class Wavefront:
    r"""
    A wavefront, expressed as a weighted sum of Zernike polynomials.
    The wavefront is returned as a `sy.Expr`. This is useful if we, for
    example, also want to compute the derivative of the wavefront.
    Both the polar and the Cartesian representation of the wavefront
    are available.

    Args:
        coefficients: The coefficients to be used as weights for the
            Zernike polynomials. There are two ways to specify this:

            1. As a sequence of floats. In this case, the :math:`j`-th
               entry of the sequence will be used as the weight for the
               :math:`j`-th Zernike polynomial. That means:

               >>> coefficients = [0, 1, 2, 0, 4]

               will produce the following wavefront:

               .. math::
                 \text{WF} = Z_1 + 2 \cdot Z_2 + 4 \cdot Z_4

            2. As a dictionary with entries of the form `(j, weight)`.
               To reproduce the previous example, we could therefore
               also write:

               >>> coefficients = {1: 1, 2: 2, 4: 4}

            Note that there is per se no limit of the number of Zernike
            polynomials that can be used for a wavefront; however,
            things of course will get slower for higher orders.
    """

    def __init__(self,
                 coefficients: Union[Sequence[float], Dict[int, float]]):

        # Store constructor arguments
        self.coefficients = coefficients

    @property
    def polar(self) -> sy.Expr:
        """
        Get the polar representation of the wavefront.

        Returns: sy.Expr
            A `sy.Expr` containing the polar representation of the
            wavefront.
        """

        if isinstance(self.coefficients, dict):
            return sum(coefficient * ZernikePolynomial(j=j).polar
                       for j, coefficient in self.coefficients.items())
        else:
            return sum(coefficient * ZernikePolynomial(j=j).polar
                       for j, coefficient in enumerate(self.coefficients))

    @property
    def cartesian(self) -> sy.Expr:
        """
        Get the Cartesian representation of the wavefront.

        Returns:
            A `sy.Expr` containing the Cartesian representation of the
            wavefront.
        """

        return polar_to_cartesian(self.polar)
