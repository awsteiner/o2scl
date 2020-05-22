:ref:`O2scl <o2scl>`

Related Projects
================

Several noteworthy related projects:

- GSL - The GNU Scientific Library \n The first truly free,
ANSI-compliant, fully-tested, and well-documented scientific
computing library. The GSL is located at
http://www.gnu.org/software/gsl , and the manual is available at \n
http://www.gnu.org/software/gsl/manual/html_node/ . Many GSL
routines are included and reworked in \o2 (the corresponding
classes begin with the \c gsl_ prefix) and \o2 was specifically
designed to be used with GSL. GSL is a must-have for most serious
numerical work in C or C++ and is required for installation of
\o2. See also \ref Galassi09 .

- HDF - Hierarchical Data Format, (http://www.hdfgroup.org) \n A
library developed for file I/O of scientific data. Contains C,
C++, FORTRAN and Java APIs, several command-line tools, and a
longstanding active community of users.

- Boost - (http://www.boost.org) \n Free (license is compatible with
GPL) peer-reviewed portable C++ source libraries that work well
with the C++ Standard Library. Boost also contains uBlas,
for linear algebra computing. 

- CERNLIB - (http://cernlib.web.cern.ch/cernlib/mathlib.html) \n The
gold standard in FORTRAN computing libraries. Several CERNLIB
routines are rewritten completely in C++ and included in \o2 (they
begin with the \c cern_ prefix).

- LAPACK and ATLAS - The gold standard for linear algebra \n
Available at http://www.netlib.org/lapack and
http://www.netlib.org/atlas . See also
http://www.netlib.org/clapack .

- QUADPACK - FORTRAN adaptive integration library,
http://www.netlib.org/quadpack \n This is the library on which the
GSL integration routines are based (prefix \c gsl_inte_).

- GNU Units - Unit conversion utility from
http://www.gnu.org/software/units/ . This utility can be used used
by \ref o2scl::convert_units in \o2 to perform unit conversions.

- Eigen - (http://eigen.tuxfamily.org) \n A linear algebra
library in C++. Eigen is licensed under MPLv2.

- Armadillo - (http://arma.sourceforge.net) \n A linear algebra
library in C++. Armadillo is licensed under MPLv2.

- FLENS - (http://flens.sourceforge.net) \n A C++ library with
object-oriented linear algebra types and an interface to BLAS
and LAPACK. 

- MESA - Modules for Experiments in Stellar Astrophysics,
(http://mesa.sourceforge.net) \n An excellent FORTRAN library with
accurate low-density equations of state, interpolation, opacities
and other routines useful in stellar physics. Work is currently
under way to rewrite some of the MESA routines in C/C++ for \o2.
Licensed with LGPL (not GPL). 

- OOL - Open Optimization Library (http://ool.sourceforge.net) \n
Constrained minimization library designed with GSL in mind. The
\o2 constrained minimization classes are derived from this
library, though OOL appears not to be under current development.
\ref conmin_section .

- Root - CERN's new C++ analysis package (http://root.cern.ch) \n
A gargantuan library for data analysis, focused mostly on
high-energy physics. Their histograms, graphics, file I/O and
support for large data sets is particularly good. 
