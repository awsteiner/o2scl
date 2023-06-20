Multiprecision Support
======================

:ref:`O2scl <o2scl>`

Some O₂scl classes support floating-point types beyond
``double``. Generally, these classes support both ``long double`` and
the ``boost::multiprecision`` types. (Though some of the boost
multiprecision types require additional libraries, such as GMP.)

List of classes which support multiprecision types:

- Function evaluation: :ref:`funct_multip <funct_multip>`
- Numerical differentiation: :ref:`deriv_gsl <deriv_gsl>` and
  :ref:`deriv_multip_gsl <deriv_multip_gsl>`
- Function approximation: :ref:`cheb_approx_tl <cheb_approx_tl>`
- Root-finding: :ref:`root_brent_gsl <root_brent_gsl>` and
  :ref:`root_toms748 <root_toms748>`
- String to double conversion: :cpp:func:`o2scl::dtos()`
- Mathematical expression evaluation :ref:`calc_utf8 <calc_utf8>`
- Integration: :ref:`inte_gauss56_cern <inte_gauss56_cern>`,
  :ref:`inte_gauss_cern <inte_gauss_cern>`, :ref:`inte_adapt_cern
  <inte_adapt_cern>`, :ref:`inte_kronrod_boost <inte_kronrod_boost>`,
  :ref:`inte_tanh_sinh_boost <inte_tanh_sinh_boost>`,
  :ref:`inte_exp_sinh_boost <inte_exp_sinh_boost>`, and
  :ref:`inte_sinh_sinh_boost <inte_sinh_sinh_boost>`.

