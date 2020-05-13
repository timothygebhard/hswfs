"""
Pre-compute Zernike functions to speed up computations.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from tqdm import tqdm

import black
import sympy as sy

from hswfs.zernike import j_to_mn, ZernikePolynomial, derive


# -----------------------------------------------------------------------------
# MAIN CODE
# -----------------------------------------------------------------------------

if __name__ == '__main__':

    # Define file header, with docstring, imports, function signature, ...
    lines = \
        ['"""',
         'This module provides pre-computed numpy versions of the Cartesian\n'
         'derivatives of Zernike polynomials, which can be used to '
         '(significantly)\nspeed up the fitting of wavefronts.',
         '',
         '.. warning::',
         '   **Note:** This file / module is automatically generated by the\n'
         '   ``precompute_fast_zernike.py`` script in the ``scripts``'
         ' directory!',
         '"""\n',
         '# ' + 77 * '-',
         '# IMPORTS',
         '# ' + 77 * '-' + '\n',
         'from typing import Union\n',
         'from numpy import ndarray, ones, sqrt, sin, cos, arctan2 as atan2',
         'from sympy import Symbol',
         '\n',
         '# ' + 77 * '-',
         '# FUNCTION DEFINITIONS',
         '# ' + 77 * '-',
         '',
         'def zernike_derivative_cartesian('
         '    m: int,'
         '    n: int,'
         '    x: Union[float, ndarray],'
         '    y: Union[float, ndarray],'
         '    wrt: Union[str, Symbol]) -> Union[float, ndarray]:',
         '    r"""',
         '    Evaluate the Cartesian derivative of a the Zernike polynomial ',
         '    :math:`Z^m_n` at the position(s) `x`, `y`. Fast.',
         '',
         '    Args:',
         '        m: Index :math:`m` of :math:`Z^m_n`.',
         '        n: Index :math:`n` of :math:`Z^m_n`. ',
         '            Default maximum value is :math:`n_\\text{max} = 15`.',
         '        x: The x-coordinate(s) at which to evaluate the derivative.',
         '        y: The y-coordinate(s) at which to evaluate the derivative.',
         '        wrt: A string or `sy.Symbol` that specifies with respect to',
         '            which variable the derivative is taken: `"x"` or `"y"`.',
         '',
         '    Returns:',
         '        The value(s) of the derivative of :math:`Z^m_n` with ',
         '        respect to `wrt` at the given position(s) `x`, `y`.',
         '    """\n']

    # Loop over Zernike polynomials, compute their derivatives and add them
    # to the function body
    print('Pre-computing derivatives of Zernike polynomials:', flush=True)
    for j in tqdm(range(136)):

        # Compute double indices from single index
        m, n = j_to_mn(j)

        # Loop over derivatives in x- and y-direction
        entry = f'\n    # Derivatives for j = {j}\n'
        for wrt in ('x', 'y'):

            # Compute the derivative as a sympy expression
            derivative = derive(ZernikePolynomial(m=m, n=n).cartesian, wrt=wrt)
            derivative = sy.nsimplify(sy.simplify(derivative))

            # Manually take care of the case where the derivative is constant
            # to make sure that the output always has the same shape as the
            # input (i.e., also works on numpy arrays)
            if derivative.is_constant():
                entry += \
                    (f""
                     f"    if m == {m} and n == {n} and wrt == '{wrt}':\n"
                     f"        if isinstance(x, ndarray):\n"
                     f"            return {derivative} * ones(x.shape)\n"
                     f"        return {derivative}\n")

            # Otherwise, we can directly use the sympy's string representation
            else:
                entry += \
                    (f""
                     f"    if m == {m} and n == {n} and wrt == '{wrt}':\n"
                     f"        return {derivative}\n")

        # Store lines with code for the current Zernike polynomial
        lines.append(entry)

    # Raise value error if the function is called with wrong arguments
    lines += \
        ['    # Raise value error if we have not returned yet',
         '    raise ValueError("No pre-computed derivative available for '
         'the given arguments!")',
         '']

    # Combine all lines into a single string
    code = '\n'.join(lines)

    # Run black on that string to automatically format it
    code = black.format_str(code, mode=black.FileMode(line_length=80))

    # Store the result as a Python file in the hswfs package directory
    with open('../hswfs/fast_zernike.py', 'w') as python_file:
        python_file.write(code)
