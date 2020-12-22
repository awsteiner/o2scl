Multiprecision Support
======================

:ref:`O2scl <o2scl>`

Some O\ :sub:`2`\ scl classes support floating-point types beyond
``double``. Generally, these classes support both ``long double`` and
the ``boost::multiprecision`` types. (Though some of the boost
multiprecision types require additional libraries, such as GMP.)

AWS 3/7/20: I don't use the boost-supplied ``std::abs()`` and
``std::isfinite()`` functions for multiprecision types because some
HPC compilers don't have recent versions of boost installed. I
have currently defined new functions :cpp:func:`o2scl::o2abs()`,
:cpp:func:`o2scl::o2isfinite()`, and
:cpp:func:`o2scl::o2hypot()` which work with both ``long double`` and the
``boost::multiprecision`` types.

Testing with multiprecision types is enabled by defining
``O2SCL_LD_TYPES`` at the time the configure script is run.
    
List of classes which support multiprecision types:

- Numerical differentiation: :ref:`deriv_gsl <deriv_gsl>`
- Function approximation: :ref:`cheb_approx_tl <cheb_approx_tl>`
- Root-finding: :ref:`root_brent_gsl <root_brent_gsl>` and
  :ref:`root_toms748 <root_toms748>`
- String to double conversion: :cpp:func:`o2scl::dtos()`
- Integration: :ref:`inte_gauss56_cern <inte_gauss56_cern>`,
  :ref:`inte_gauss_cern <inte_gauss_cern>`, :ref:`inte_adapt_cern
  <inte_adapt_cern>`, :ref:`inte_kronrod_boost <inte_kronrod_boost>`,
  :ref:`inte_tanh_sinh_boost <inte_tanh_sinh_boost>`,
  :ref:`inte_exp_sinh_boost <inte_exp_sinh_boost>`, and
  :ref:`inte_sinh_sinh_boost <inte_sinh_sinh_boost>`.
