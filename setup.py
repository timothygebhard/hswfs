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
        'astropy>=4.0.1.post1',
        'black>=19.10b0',
        'matplotlib>=3.2.1',
        'numpy>=1.18.4',
        'scipy>=1.4.1',
        'sympy>=1.5.1',
        'tqdm>=4.46.0',
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
