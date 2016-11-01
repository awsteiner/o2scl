from setuptools import setup

setup(name='o2sclpy',version='0.920',
      description='Python extensions for O2scl',
      url='http://web.utk.edu/~asteine1/o2scl',
      author='Andrew W. Steiner',
      author_email='awsteiner@mykolab.com',license='GPLv3',
      packages=['o2sclpy'],install_requires=['h5py','numpy','matplotlib'],
      zip_safe=False,scripts=['bin/o2graph'],
      classifiers=[
          'Development Status :: 4 - Beta',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3 :: Only',
          'Topic :: Scientific/Engineering'
      ])
