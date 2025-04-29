Special Functions
=================

:ref:`O2scl <o2scl>`

Many mathematical functions are now provided by the GSL and Boost
libraries, on which O₂scl depends. A few useful functions not provided
by these two are included in O₂scl.

These classes and functions depend on the Boost multiprecision
libraries, thus may not work well if O2SCL_NO_BOOST_MULTIPRECISION is
defined. 

Polylogs and Related Functions
------------------------------

Some of the Fermi-Dirac integrals can be computed by the GSL library,
and the class :ref:`fermi_dirac_integ_gsl <fermi_dirac_integ_gsl>`
provides a C++ interface for those functions. The class
:ref:`fermi_dirac_multip <fermi_dirac_multip>` uses multiprecision to
compute Fermi-Diract integrals to the full floating-point precision of
the user-specified type, up to approximately 45 digits. The class
:ref:`bose_einstein_multip <bose_einstein_multip>` works similarly to
compute Bose-Einstein integrals up to approximately 45 digits.

Approximate polylogarithms can be computed by the :ref:`polylog
<polylog>` class. This class is relatively fast, but cannot
always compute polylogs to full double precision. The
template class :ref:`polylog_multip <polylog_multip>` is a bit
slower, but can compute polylogs to high precision. The
template specialization ``polylog_multip<double,double>``,
for example, will compute polylogs to full double precision.

Scaled Bessel Functions
-----------------------

The particular combination :math:`K_{\nu}(x) \exp(x)` is useful for
particle thermodynamics The functions for :math:`\nu=1, 2` and
:math:`3` can be computed to double prectision by the GSL library, and
an interface is provided in :ref:`bessel_K_exp_integ_gsl
<bessel_K_exp_integ_gsl>`. The class :ref:`bessel_K_exp_integ_boost
<bessel_K_exp_integ_boost>` computes these functions to high precision
using Boost.

