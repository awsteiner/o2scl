Related Projects
================

:ref:`O2scl <o2scl>`

Several noteworthy related projects:

- GNU Scientific Library (GSL): The first truly free, ANSI-compliant,
  fully-tested, and well-documented scientific computing library. The
  GSL is located at https://www.gnu.org/software/gsl , and the manual
  is available at https://www.gnu.org/software/gsl/manual/html_node/ .
  Many GSL routines are included and reworked in O₂scl. GSL
  is required for installation of O₂scl. See also
  [Galassi09]_.
- Hierarchical Data Format (HDF; https://www.hdfgroup.org) A library
  developed for file I/O of scientific data. Contains C, C++, FORTRAN
  and Java APIs, several command-line tools, and a longstanding active
  community of users.
- Boost (https://www.boost.org) Free (license is compatible with GPL)
  peer-reviewed portable C++ source libraries that work well with the
  C++ Standard Library. Boost also contains uBlas, for linear algebra
  computing.
- CERNLIB: (https://cernlib.web.cern.ch/cernlib/mathlib.html) A
  collection of FORTRAN functions for scientific computing. This
  library has apparently disappeared from CERN's web pages. Several
  CERNLIB routines are rewritten completely in C++ and included in O\
  :sub:`2`\ scl.
- LAPACK and ATLAS - The gold standard for linear algebra
  Available at https://www.netlib.org/lapack and
  https://www.netlib.org/atlas . See also
  https://www.netlib.org/clapack .
- QUADPACK - FORTRAN adaptive integration library,
  https://www.netlib.org/quadpack This is the library on which the
  GSL integration routines are based.
- GNU Units - Unit conversion utility from
  https://www.gnu.org/software/units/ . This utility can be used used
  by \ref o2scl::convert_units in O₂scl to perform unit
  conversions.
- Eigen - (https://eigen.tuxfamily.org) A linear algebra library in
  C++. Eigen is licensed under MPLv2.
- Armadillo - (https://arma.sourceforge.net) A linear algebra library
  in C++. Armadillo is licensed under MPLv2.
- OOL - Open Optimization Library (https://ool.sourceforge.net)
  Constrained minimization library designed with GSL in mind. The O\
  :sub:`2`\ scl constrained minimization classes are derived from this
  library (see :ref:`Constrained Minimization`). OOL appears not to be
  under current development.
- Root - CERN's new C++ analysis package (https://root.cern.ch)
  A gargantuan library for data analysis, focused mostly on
  high-energy physics. Their histograms, graphics, file I/O and
  support for large data sets is particularly good. 
