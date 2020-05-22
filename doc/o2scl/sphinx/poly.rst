:ref:`O2scl <o2scl>`

Roots of Polynomials
====================

Classes are provided for solving quadratic, cubic, and quartic
equations as well as general polynomials. There is a standard
nomenclature: classes which handle polynomials with real
coefficients and real roots end with the suffix <tt>_real</tt>
(\ref o2scl::quadratic_real, \ref o2scl::cubic_real and \ref
o2scl::quartic_real), classes which handle real coefficients and
complex roots end with the suffix <tt>_real_coeff</tt> (\ref
o2scl::quadratic_real_coeff, \ref o2scl::cubic_real_coeff, \ref
o2scl::quartic_real_coeff, and \ref o2scl::poly_real_coeff), and
classes which handle complex polynomials with complex coefficients
end with the suffix <tt>_complex</tt> (\ref o2scl::quadratic_complex,
\ref o2scl::cubic_complex, \ref o2scl::quartic_complex, and 
\ref o2scl::poly_complex). As a reminder,
complex roots may not occur in conjugate pairs if the coefficients
are not real. Most of these routines will return an error if the
leading coefficient is zero.

In the public interfaces to the polynomial solvers, the
complex type <tt>std::complex<double></tt> is used. 
\comment
These can
be converted to and from the GSL complex type using the 
o2scl::complex_to_gsl() and o2scl::gsl_to_complex() functions.
\endcomment

For quadratics, \ref o2scl::quadratic_real_coeff_gsl is the best if the
coefficients are real, while if the coefficients are complex, use
\ref o2scl::quadratic_complex_std. For cubics with real coefficients,
\ref o2scl::cubic_real_coeff_cern is the best, while if the coefficients
are complex, use \ref o2scl::cubic_complex_std.
    
For a quartic polynomial with real coefficients, 
\ref o2scl::quartic_real_coeff_cern is the best, unless the
coefficients of odd powers happen to be small, in which case, 
\ref o2scl::quartic_real_gsl2 tends to work better. For quartics,
generic polynomial solvers such as \ref o2scl::poly_real_coeff_gsl
can provide more accurate (but slower) results. If the
coefficients are complex, then you can use \ref
o2scl::quartic_complex_simple.

\section ex_poly_sect Polynomial solver example

This example shows how to find the roots of the 
second-, third-, fourth-, and fifth-order 
Chebyshev polynomials
\f{eqnarray*}
&2x^2-1& \nonumber \\
&4x^3-3 x& \nonumber \\
&8x^4-8x^2+1& \nonumber \\
&16x^5-20x^3+5x& \nonumber
\f}
For the Chebyshev polynomial of order \f$ n \f$, the 
roots are given by
\f[
\cos \left[ \frac{\pi(k-1/2)}{n}\right]
\f]
for \f$ k = 1,\ldots,n \f$ These roots are used in
\ref o2scl::cheb_approx_tl to approximate functions using 
Chebyshev polynomials .

\dontinclude ex_poly.cpp
\skip Example:
\until End of example
