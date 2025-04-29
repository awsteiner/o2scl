Roots of Polynomials
====================

:ref:`O2scl <o2scl>`

Contents
--------

- :ref:`Introduction`
- :ref:`Polynomial solver example`

Introduction
------------

Classes are provided for solving quadratic, cubic, and quartic
equations as well as general polynomials. There is a standard
nomenclature: classes which handle polynomials with real coefficients
and real roots end with the suffix ``_real`` (:ref:`quadratic_real
<quadratic_real>`, :ref:`cubic_real <cubic_real>` and
:ref:`quartic_real <quartic_real>`), classes which handle real
coefficients and complex roots end with the suffix ``_real_coeff`` (
:ref:`quadratic_real_coeff <quadratic_real_coeff>`,
:ref:`cubic_real_coeff <cubic_real_coeff>`, :ref:`quartic_real_coeff
<quartic_real_coeff>`, and :ref:`poly_real_coeff <poly_real_coeff>`),
and classes which handle complex polynomials with complex coefficients
end with the suffix ``_complex`` (:ref:`quadratic_complex
<quadratic_complex>`, :ref:`cubic_complex <cubic_complex>`,
:ref:`quartic_complex <quartic_complex>`, and :ref:`poly_complex
<poly_complex>`). As a reminder, complex roots may not occur in
conjugate pairs if the coefficients are not real. Most of these
routines will return an error if the leading coefficient is zero.

In the public interfaces to the polynomial solvers, the
complex type ``std::complex<double>`` is used. 

.. 
  These can
  be converted to and from the GSL complex type using the 
  o2scl::complex_to_gsl() and o2scl::gsl_to_complex() functions.

For quadratics, :ref:`quadratic_real_coeff_gsl2
<quadratic_real_coeff_gsl2>` is the best if the coefficients are real,
while if the coefficients are complex, use :ref:`quadratic_complex_std
<quadratic_complex_std>`. Both of these classes work with
multiprecision floating point types.

For cubics with real coefficients, :ref:`cubic_real_coeff_gsl
<cubic_real_coeff_gsl>`, :ref:`cubic_real_coeff_gsl2
<cubic_real_coeff_gsl2>`, and :ref:`cubic_real_coeff_cern
<cubic_real_coeff_cern>` are nearly equivalent. the GSL-based classes
are a bit faster and more accurate, but the CERNLIB-based class is a
bit better at identifying the nature of the roots. If the coefficients
are complex, use :ref:`cubic_complex_std <cubic_complex_std>`, but
keep in mind this class sometimes fails to properly identify the
nature of the roots. The classes :ref:`cubic_real_coeff_gsl2
<cubic_real_coeff_gsl2>`, :ref:`cubic_real_coeff_cern
<cubic_real_coeff_cern>`, and :ref:`cubic_complex_std
<cubic_complex_std>` all work with multiprecision floating point
types.

Currently, quartics with real coefficients are best
solved with :ref:`poly_real_coeff_gsl <poly_real_coeff_gsl>`,
albeit a bit more slowly than direct solvers. The direct
solvers are faster, but fail unpredictably, and are thus
still experimental.

..
  For a quartic polynomial with real coefficients,
  :ref:`quartic_real_coeff_cern <quartic_real_coeff_cern>` is the best,
  unless the coefficients of odd powers happen to be small, in which
  case, :ref:`quartic_real_gsl2 <quartic_real_gsl2>` tends to work
  better. For quartics, generic polynomial solvers such as
  :ref:`poly_real_coeff_gsl <poly_real_coeff_gsl>` can provide more
  accurate (but slower) results. If the coefficients are complex, then
  you can use :ref:`quartic_complex_simple <quartic_complex_simple>`.

Polynomial solver example
-------------------------

This example shows how to find the roots of the 
second-, third-, fourth-, and fifth-order 
Chebyshev polynomials

.. math:: 

   \begin{eqnarray}
   &2x^2-1& \nonumber \\
   &4x^3-3 x& \nonumber \\
   &8x^4-8x^2+1& \nonumber \\
   &16x^5-20x^3+5x& \nonumber
   \end{eqnarray}
   
For the Chebyshev polynomial of order :math:`n`, the roots are given
by

.. math:: 

   \cos \left[ \frac{\pi(k-1/2)}{n}\right]

for :math:`k = 1,\ldots,n` These roots are used in
:ref:`cheb_approx_tl <cheb_approx_tl>` to approximate functions using
Chebyshev polynomials .

.. literalinclude:: ../../../examples/ex_poly.cpp
   :language: c++		    
   :start-after: sphinx-example-start
