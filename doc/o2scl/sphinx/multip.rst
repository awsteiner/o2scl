Multiprecision Support
======================

:ref:`O2scl <o2scl>`

Some O₂scl classes support floating-point types beyond
``double``. Generally, these classes support both ``long double`` and
the ``boost::multiprecision`` types. (Though some of the boost
multiprecision types require additional libraries, such as GMP.)
See also :ref:`Multiprecision Support for Particles and EOSs`. 

List of classes which support multiprecision types:

- Function evaluation: :ref:`funct_multip_tl <funct_multip_tl>`
- Numerical differentiation: :ref:`deriv_gsl <deriv_gsl>` and
  :ref:`deriv_multip_gsl <deriv_multip_gsl>`
- Function approximation: :ref:`cheb_approx_tl <cheb_approx_tl>`
- Root-finding: :ref:`root_brent_gsl <root_brent_gsl>` and
  :ref:`root_toms748 <root_toms748>`
- String to and from floating-point conversions: :cpp:func:`o2scl::dtos()`,
  :cpp:func:`o2scl::function_to_fp()`, and 
  :cpp:func:`o2scl::function_to_fp_nothrow()`.
- Polynomial solving:
  :ref:`quadratic_real_coeff_gsl2<quadratic_real_coeff_gsl2>`,
  :ref:`quadratic_complex_std <quadratic_complex_std>`,
  :ref:`cubic_real_coeff_cern <cubic_real_coeff_cern>`,
  :ref:`cubic_real_coeff_gsl2 <cubic_real_coeff_gsl2>` and
  :ref:`cubic_complex_std <cubic_complex_std>`.
- Mathematical expression evaluation :ref:`calc_utf8 <calc_utf8>`
- Constant library :ref:`find_constants <find_constants>`
- Integration: :ref:`inte_gauss56_cern <inte_gauss56_cern>`,
  :ref:`inte_gauss_cern <inte_gauss_cern>`, :ref:`inte_adapt_cern
  <inte_adapt_cern>`, :ref:`inte_kronrod_boost <inte_kronrod_boost>`,
  and :ref:`inte_double_exp_boost <inte_double_exp_boost>`.

Many of the vector and matrix functions in :ref:`Arrays, Vectors,
Matrices, and Tensors` section also support multiprecision.       
       
Multiprecision function typedefs
--------------------------------

.. _funct_ld:

.. doxygentypedef:: funct_ld

.. _funct_cdf25:

.. doxygentypedef:: funct_cdf25

.. _funct_mpfr25:

.. doxygentypedef:: funct_mpfr25

