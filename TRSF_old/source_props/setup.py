from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension("region_expansion", ["region_expansion.pyx"],
              include_dirs=[np.get_include()])  # Add numpy include directory
]

setup(
    ext_modules=cythonize(extensions),
    include_dirs=[np.get_include()]  # Add numpy include directory
)