#!/usr/bin/env python

"""
setup.py file for SWIG o2scl
"""
import platform
from distutils.core import setup, Extension

idirs=[]
ldirs=[]

# Presume a linux system so add the standard HDF5 directories
if platform.system()=='Darwin':
    print('\nPresuming a OS X system and adding standard',
          'directories. If the installation fails to find either',
          'the GSL, HDF5 or O2scl libraries, then',
          'you may need to specify the directories using, e.g.',
          'python setup.py install',
          '--library-dirs=/library/dir/1:/library/dir/2',
          '--include-dirs=/include/dir/1:/include/dir/2\n')
    # This is normally automatically included, but if the user
    # specifies a --include-dirs argument then setup.py doesn't seem
    # to work without explicitly including this directory here.
    idirs.append('/usr/local/include')
    idirs.append('../../include')
    ldirs.append('/Users/awsteiner/install/o2scl-0.924/lib')
else:
    print('\nPresuming a Linux system and adding standard',
          'directories. If the installation fails to find either',
          'the GSL, HDF5 or O2scl libraries, then',
          'you may need to specify the directories using, e.g.',
          'python setup.py install',
          '--library-dirs=/library/dir/1:/library/dir/2',
          '--include-dirs=/include/dir/1:/include/dir/2\n')
    idirs.append('/usr/lib/x86_64-linux-gnu/hdf5/serial/include')
    ldirs.append('/usr/lib/x86_64-linux-gnu/hdf5/serial/lib')

o2scl_module = Extension('_o2scl',
                           sources=['python.cpp'],
                           include_dirs=idirs,
                           library_dirs=ldirs,
                           swig_opts=['-c++'],
                           libraries=['o2scl']
)

setup(name='o2scl',
      version='0.1',
      author='Andrew W. Steiner',
      description="""Simple swig o2scl from docs""",
      ext_modules=[o2scl_module],
      py_modules=['o2scl']
)
