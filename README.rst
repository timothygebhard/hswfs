.. raw:: html

  <h1 align="center">
    hswfs: A virtual Hartmann-Shack wavefront sensor
  </h1>
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
------------

To get started, clone this repository and install hswfs as a Python package:

.. code-block:: bash

   git clone git@github.com:timothygebhard/hswfs.git
   cd hswfs
   pip install .

🔭 Example
----------

Here is an example result, showing a virtual wavefront sensor with the shifts measured in each subaperture, and the respective reconstructed wavefront and point spread function:

.. raw:: html

  <p align="center">
    <img src="./demo/example_result.png" alt="Example Results" width="600">
  </p>

This example was created by the ``hswfs_example.py`` script in the ``demo`` directory.


📚 Documentation
----------------

The documentation can be build by running the following command in the ``docs`` folder:

.. code-block:: bash

   make html

This will create a ``build`` folder in the ``docs`` directory which contains an HTML version of the documentation.


🐭 Tests
--------

hswfs comes with a series of unit tests, which can be run as follows:

.. code-block:: bash

   pytest tests


🎁 Contributing to hswfs
------------------------

Contributions are always very welcome!
Whether you have found a bug, or want to enhance hswfs's functionality, please feel free to open an issue here on GitHub, or send a pull request.


⚠️ License and Warranty
-----------------------

hswfs is provided under the permissive `MIT License <https://github.com/timothygebhard/hswfs/blob/master/LICENSE>`_, which gives you a lot of freedom to use it for your own work.

Please note, however, that hswfs is provided *as is*, without any guarantees regarding completeness or correctness.
If you break your telescope or AO system, that's on you! 😉
