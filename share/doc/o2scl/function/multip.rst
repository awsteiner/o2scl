Multiprecision Support
======================

:ref:`O2scl <o2scl>`

Some O₂scl classes support floating-point types beyond ``double``.
Generally, these classes support both ``long double`` and the
``boost::multiprecision`` types. (Though some of the boost
multiprecision types require additional libraries, such as GMP.) In
some cases, O₂scl classes use *adaptive multiprecision*, i.e., they
are able to automatically choose among several different floating
point types to perform a calculation to a specified accuracy. See also
:ref:`Multiprecision Support for Particles and EOSs`.

Unfortunately, the multiprecision support also dramatically increases
the compilation time and the associated memory requirements.
Multiprecision support can be disabled with the
``O2SCL_NO_BOOST_MULTIPRECISION`` flag.

List of classes which support adaptive multiprecision:

- Function evaluation: :ref:`funct_multip_tl <funct_multip_tl>`
- Numerical differentiation: :ref:`deriv_gsl <deriv_gsl>`
- Root-finding: :ref:`root_brent_gsl <root_brent_gsl>` and
  :ref:`root_cern <root_cern>`
- Integration: :ref:`inte_kronrod_boost <inte_kronrod_boost>`,
  :ref:`inte_double_exp_boost <inte_double_exp_boost>`.
  :ref:`inte_gauss_cern <inte_gauss_cern>`, and
  :ref:`inte_adapt_cern_tl <inte_adapt_cern_tl>`.
  
List of classes which support multiprecision types:

- Numerical differentiation: :ref:`deriv_gsl <deriv_gsl>`
- Function approximation: :ref:`cheb_approx_tl <cheb_approx_tl>`
- Root-finding: :ref:`mroot_cern <mroot_cern>`, 
  :ref:`root_brent_gsl <root_brent_gsl>`, and
  :ref:`root_toms748 <root_toms748>`
- String to and from floating-point conversions:
  :cpp:func:`o2scl::dtos()`, :cpp:func:`o2scl::function_to_fp()`, and 
  :cpp:func:`o2scl::function_to_fp_nothrow()`.
- Polynomial solving:
  :ref:`quadratic_real_coeff_gsl2<quadratic_real_coeff_gsl2>`,
  :ref:`quadratic_complex_std <quadratic_complex_std>`,
  :ref:`cubic_real_coeff_cern <cubic_real_coeff_cern>`,
  :ref:`cubic_real_coeff_gsl2 <cubic_real_coeff_gsl2>` and
  :ref:`cubic_complex_std <cubic_complex_std>`.
- Mathematical expression evaluation :ref:`calc_utf8 <calc_utf8>`
- Constant library :ref:`find_constants <find_constants>`
- Integration: :ref:`inte_gauss56_cern <inte_gauss56_cern>`.

Almost all of the vector and matrix functions in :ref:`Arrays,
Vectors, Matrices, and Tensors` section support multiprecision.
Finally, the BLAS functions documented at :ref:`BLAS functions`
support multiprecision.

Multiprecision function typedefs
--------------------------------

.. _funct_ld:

.. doxygentypedef:: funct_ld

.. _funct_cdf25:

.. doxygentypedef:: funct_cdf25

.. _funct_cdf35:

.. doxygentypedef:: funct_cdf35

.. _funct_cdf50 :

.. doxygentypedef:: funct_cdf50 

.. _funct_cdf100:

.. doxygentypedef:: funct_cdf100

.. _funct_ret_ld:

.. doxygentypedef:: funct_ret_ld

.. _funct_ret_cdf25:

.. doxygentypedef:: funct_ret_cdf25

.. _funct_ret_cdf35:

.. doxygentypedef:: funct_ret_cdf35

.. _funct_ret_cdf50 :

.. doxygentypedef:: funct_ret_cdf50 

.. _funct_ret_cdf100:

.. doxygentypedef:: funct_ret_cdf100

.. _funct_mpfr25:

.. doxygentypedef:: funct_mpfr25

.. _funct_mpfr35:

.. doxygentypedef:: funct_mpfr35

.. _funct_mpfr50:

.. doxygentypedef:: funct_mpfr50

.. _funct_mpfr100:

.. doxygentypedef:: funct_mpfr100

.. _mm_funct_ld:

.. doxygentypedef:: mm_funct_ld

.. _mm_funct_cdf25:

.. doxygentypedef:: mm_funct_cdf25

.. _mm_funct_cdf35:

.. doxygentypedef:: mm_funct_cdf35

.. _mm_funct_cdf50 :

.. doxygentypedef:: mm_funct_cdf50 

.. _mm_funct_cdf100:

.. doxygentypedef:: mm_funct_cdf100
