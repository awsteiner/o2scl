Differentiation
===============

Differentiation is performed by descendants of \ref o2scl::deriv_base. These
allow one to calculate either first, second, and third
derivatives. A GSL-based routine is used in \ref o2scl::deriv_gsl, and
the CERNLIB routine is used in \ref o2scl::deriv_cern. For functions 
which are tabulated over equally-spaced
abscissas, the class \ref o2scl::deriv_eqi is provided which applies the
formulas from Abramowitz and Stegun at a specified order. The
class \ref o2scl::deriv_cern is slower and sometimes more accurate, but
also fails more often than \ref o2scl::deriv_gsl, which never calls the
error handler.

\b Warning: For \ref o2scl::deriv_gsl and \ref o2scl::deriv_cern,
the second and third derivatives are calculated by naive repeated
application of the code for the first derivative and can be
particularly troublesome if the function is not sufficiently
smooth. Error estimation is not provided for second and third
derivatives.

\section ex_deriv_sect Differentiation example

This example computes first and second derivatives of
\f[
y(x) = \sin (2 x) + \frac{1}{2}
\f]
with both \ref o2scl::deriv_gsl and \ref o2scl::deriv_cern .

\dontinclude ex_deriv.cpp
\skip Example:
\until End of example
