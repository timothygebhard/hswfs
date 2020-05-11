.. raw:: html

  <h2 align="center">
    hswfs: A virtual Hartmann-Shack wavefront sensor
  </h2>
  <p align="center">
    <img src="https://img.shields.io/badge/python-v3.7-blue" alt="Python 3.7">
    <img src="https://img.shields.io/badge/license-MIT-green" alt="License: MIT">
  </p>

hswfs is a Python package which provides a minimalistic simulation of a `Hartmann-Shack Wavefront Sensor <https://en.wikipedia.org/wiki/Shack%E2%80%93Hartmann_wavefront_sensor>`_.
A HSWFS is a popular and conceptually simple type of `wavefront sensor <https://en.wikipedia.org/wiki/Wavefront_sensor>`_ that can be used to measure the abberations of a `wavefront <https://en.wikipedia.org/wiki/Wavefront>`_.
It is commonly used in the `adaptive optics systems <https://en.wikipedia.org/wiki/Adaptive_optics>`_ of larger telescopes, which correct the image for the perturbations introduced by the turbulent movement of the air in the Earth's atmosphere.


**Among other things, hswfs includes the following features:**

- Generate random wavefront sensor data.
- Fit the wavefront using `Zernike polynomials <https://en.wikipedia.org/wiki/Zernike_polynomials>`_ of *arbitrary* order (hswfs uses `sympy <https://sympy.org/>`_ to compute expressions of arbitrary Zernike polynomials and their derivatives in both polar and Cartesian coordinates).
- Compute the `point spread function (PSF) <https://en.wikipedia.org/wiki/Point_spread_function>`_ from a wavefront.


⚡ Quickstart
-------------

To get started, clone this repository and install hswfs as a Python package:

.. code-block:: bash

   git clone git@github.com:timothygebhard/hswfs.git
   cd hswfs
   pip install .


🔭 Example
----------

Here are the first 15 Zernike polynomials, as computed by hswfs:

.. raw:: html

  <p align="center">
    <img src="./demo/zernike_polynomials.png" alt="Zernike Polynomials">
  </p>
