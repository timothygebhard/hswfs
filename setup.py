"""
Set up hswfs as a Python package.
"""

# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------

from os.path import join, dirname
from setuptools import find_packages, setup


# -----------------------------------------------------------------------------
# PRELIMINARIES
# -----------------------------------------------------------------------------

# Read in current package version
with open(join(dirname(__file__), "hswfs/VERSION")) as version_file:
    version = version_file.read().strip()


# -----------------------------------------------------------------------------
# RUN setup() FUNCTION
# -----------------------------------------------------------------------------

setup(
    name='hswfs',
    version=version,
    description='A virtual Hartmann-Shack wavefront sensor',
    author='Timothy Gebhard',
    author_email='mail@timothygebhard.de',
    url='https://github.com/timothygebhard/hswfs',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'astropy>=5.1',
        'black>=22.6.0',
        'matplotlib>=3.5.3',
        'numpy>=1.23.2',
        'scipy>=1.9.0',
        'sympy>=1.10.1',
        'tqdm>=4.64.0',
    ],
    license='MIT',
    zip_safe=False,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
)
