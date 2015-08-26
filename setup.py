#!/usr/bin/python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

_extra = ['-fopenmp', '-O3', '-ffast-math']


extensions = [
  Extension('zonemap3d',
            sources = ['./src/zonemap3d.pyx'],
            extra_compile_args = _extra
  )
]

setup(
  name = 'zonemap3d',
  version = '0.0.1',
  author = '@inconvergent',
  install_requires = ['numpy>=1.8.2', 'cython>=0.20.0'],
  license = 'MIT',
  cmdclass={'build_ext' : build_ext},
  ext_modules = cythonize(extensions)
)
