Integration
===========

:ref:`O2scl <o2scl>`

Integration contents
--------------------

- :ref:`One-dimensional integration based on GSL, CERNLIB, and Boost`
- :ref:`GSL-based integration details`
- :ref:`GSL-based integration error messages`
- :ref:`One-dimensional integration example`
- :ref:`Gauss-Kronrod integration coefficients`
- :ref:`Non-adaptive quadrature integration coefficients`

..
  This section is commented out for now

..
  - :ref:`Multi-dimensional integration routines`
  
One-dimensional integration based on GSL, CERNLIB, and Boost
------------------------------------------------------------

Several classes integrate arbitrary one-dimensional functions:

- Integration over a finite interval:
  
  * :ref:`inte_adapt_cern_tl <inte_adapt_cern_tl>`
  * :ref:`inte_gauss_cern <inte_gauss_cern>`
  * :ref:`inte_gauss56_cern <inte_gauss56_cern>`
  * :ref:`inte_kronrod_boost <inte_kronrod_boost>`
  * :ref:`inte_double_exp_boost <inte_double_exp_boost>`
  * :ref:`inte_qag_gsl <inte_qag_gsl>`
  * :ref:`inte_qng_gsl <inte_qng_gsl>`
    
- Integration from :math:`a` to :math:`\infty`:
  
  * :ref:`inte_qagiu_gsl <inte_qagiu_gsl>`
  * :ref:`inte_double_exp_boost <inte_double_exp_boost>`
  * :ref:`inte_transform <inte_transform>`
    
- Integration from  :math:`-\infty` to :math:`b`:
  
  * :ref:`inte_qagil_gsl <inte_qagil_gsl>`
  * :ref:`inte_double_exp_boost <inte_double_exp_boost>`
  * :ref:`inte_transform <inte_transform>`
       
- Integration from  :math:`-\infty` to :math:`\infty`:
  
  * :ref:`inte_qagi_gsl <inte_qagi_gsl>`,
  * :ref:`inte_double_exp_boost <inte_double_exp_boost>`
  * :ref:`inte_transform <inte_transform>`.
       
- Integration over a finite interval for a function with
  singularities:

  * :ref:`inte_qags_gsl <inte_qags_gsl>` 
  * :ref:`inte_qaws_gsl <inte_qaws_gsl>`
  * :ref:`inte_double_exp_boost <inte_double_exp_boost>` can handle 
    singularities at either endpoint.
    
- Cauchy principal value integration over a finite interval:
  
  * :ref:`inte_cauchy_cern <inte_cauchy_cern>`
  * :ref:`inte_qawc_gsl <inte_qawc_gsl>`
    
- Integration over a function weighted by ``cos(x)`` or ``sin(x)``:
  
  * :ref:`inte_qawo_gsl_cos <inte_qawo_gsl_cos>`
  * :ref:`inte_qawo_gsl_sin <inte_qawo_gsl_sin>`
       
- Fourier integrals:

  * :ref:`inte_qawf_gsl_cos <inte_qawf_gsl_cos>`
  * :ref:`inte_qawf_gsl_sin <inte_qawf_gsl_sin>`
       
- Integration over a weight function

.. math::
   
  W(x)=(x-a)^{\alpha}(b-x)^{\beta}\log^{\mu}(x-a)\log^{\nu}(b-x)

is performed by :ref:`inte_qaws_gsl <inte_qaws_gsl>`. 

Note that some of the integrators support multiprecision,
see :ref:`Multiprecision Support`.

There are two competing factors that can slow down an adaptive
integration algorithm: (1) each evaluation of the integrand can be
numerically expensive, depending on how the function is defined,
and (2) the process of subdividing regions and recalculating
values is almost always numerically expensive in its own right.
For integrands that are very smooth (e.g., analytic functions), a
high-order Gauss-Kronrod rule (e.g., 61-point) will achieve the
desired error tolerance with relatively few subdivisions. For
integrands with discontinuities or singular derivatives, a
low-order rule (e.g., 15-point) is often more efficient. 

GSL-based integration details
-----------------------------

For the GSL-based integration routines, the variables
:cpp:var:`o2scl::inte::tol_abs` and :cpp:var:`o2scl::inte::tol_rel`
have the same role as the quantities usually denoted in the GSL
integration routines by ``epsabs`` and ``epsrel``. In particular, the
integration classes attempt to ensure that

.. math::

  |\mathrm{result}-I| \leq \mathrm{Max}(\mathrm{tol\_abs},
  \mathrm{tol\_rel}|I|)

and returns an error to attempt to ensure that

.. math::
  
  |\mathrm{result}-I| \leq \mathrm{abserr} \leq
  \mathrm{Max}(\mathrm{tol\_abs},\mathrm{tol\_rel}|I|)

where ``I`` is the integral to be evaluated. Even when the
corresponding descendant of :cpp:func:`o2scl::inte::integ()` returns
success, these inequalities may fail for sufficiently difficult
functions. All of the GSL integration routines except for
:ref:`inte_qng_gsl <inte_qng_gsl>` use a workspace given in
:ref:`inte_workspace_gsl <inte_workspace_gsl>` which holds the results
of the various subdivisions of the original interval.

The GSL routines were originally based on QUADPACK, which is
available at http://www.netlib.org/quadpack . 

For adaptive GSL integration classes, the type of Gauss-Kronrod
quadrature rule that is used to approximate the integral and estimate
the error of a subinterval is set by
:cpp:func:`o2scl::inte_kronrod_gsl::set_rule()`.
    
The number of subdivisions of the integration region is limited by the
size of the workspace, set in
:cpp:func:`o2scl::inte_kronrod_gsl::set_limit()`. The number of
subdivisions required for the most recent call to
:cpp:func:`o2scl::inte::integ()` or
:cpp:func:`o2scl::inte::integ_err()` is given in
:cpp:var:`o2scl::inte::last_iter`. This number will always be less
than or equal to the workspace size.
         
.. note:: The GSL integration routines can sometimes lose precision if
	  the integrand is everywhere much smaller than unity. Some
	  rescaling may be required in these cases.

GSL-based integration error messages
------------------------------------

The error messages given by the adaptive GSL integration routines tend
to follow a standard form and are documented here. There are several
error messages which indicate improper usage and cause the error
handler to be called regardless of the value of
:cpp:var:`o2scl::inte::err_nonconv`:

- ``Iteration limit exceeds workspace in
  class::function().`` [ :cpp:enumerator:`exc_einval` ] 
- ``Could not integrate function in class::function() (it may have
  returned a non-finite result).`` [ :cpp:enumerator:`exc_efailed` ]
  This often occurs when the user-specified function returns ``inf``
  or ``nan``.
- ``Tolerance cannot be achieved with given value of tol_abs and
  tol_rel in class::function().`` [ :cpp:enumerator:`exc_ebadtol`
  ] This occurs if the user supplies unreasonable values for
  :cpp:var:`o2scl::inte::tol_abs` and :cpp:var:`o2scl::inte::tol_rel`.
  All positive values for :cpp:var:`o2scl::inte::tol_abs` are allowed.
  If zero-tolerance for :cpp:var:`o2scl::inte::tol_abs` is desired,
  then :cpp:var:`o2scl::inte::tol_rel` must be at least :math:`50
  \cdot \epsilon_\mathrm{mach}` (:math:`\approx 1.11 \times 10^{-14}`
  ).
- ``Cannot integrate with singularity on endpoint in
  inte_qawc_gsl::qawc().`` [ :cpp:enumerator:`exc_einval` ] The
  class \ref o2scl::inte_qawc_gsl cannot handle the case when a
  singularity is one of the endpoints of the integration.

There are also convergence errors which will call the error handler
unless :cpp:var:`o2scl::inte::err_nonconv` is false. See :ref:`What is
an error?` for more discussion on convergence errors versus fatal
errors:

- ``Cannot reach tolerance because of roundoff error on first attempt
  in class::function().`` [ :cpp:enumerator:`exc_eround` ] Each
  integration attempt tests for round-off errors by comparing the
  computed integral with that of the integrand's absolute value (i.e.,
  :math:`L^1` -norm). A highly oscillatory integrand may cause this
  error.
- ``A maximum of 1 iteration was insufficient in
  class::function().`` [ :cpp:enumerator:`exc_emaxiter` ] This occurs
  if the workspace is allocated for one interval and a single
  Gauss-Kronrod integration does not yield the accuracy demanded by
  :cpp:var:`o2scl::inte::tol_abs` and :cpp:var:`o2scl::inte::tol_rel`.
- ``Bad integrand behavior in class::function().`` [
  :cpp:enumerator:`exc_esing` ] This occurs if the integrand is
  (effectively) singular in a region, causing the subdivided intervals
  to become too small for floating-point precision.
- ``Maximum number of subdivisions 'value' reached in
  class::function().`` [ :cpp:enumerator:`exc_emaxiter` ] This occurs
  if the refinement algorithm runs out of allocated workspace. The
  number of iterations required for the most recent call to
  :cpp:func:`o2scl::inte::integ()` or
  :cpp:func:`o2scl::inte::integ_err()` is given in
  :cpp:var:`o2scl::inte::last_iter`. This number will always be less
  than or equal to the workspace size.
- ``Roundoff error prevents tolerance from being achieved in
  class::function().`` [ :cpp:enumerator:`exc_eround` ] The refinement
  procedure counts round-off errors as they occur and terminates
  if too many such errors accumulate.
- ``Roundoff error detected in extrapolation table in 
  inte_singular_gsl::qags().`` [ :cpp:enumerator:`exc_eround` ] This occurs
  when error-terms from the :math:`\varepsilon` -algorithm
  are are monitored and compared with the error-terms from the
  refinement procedure. The algorithm terminates if these
  sequences differ by too many orders of magnitude. See 
  :cpp:func:`o2scl::inte_singular_gsl::qelg()`.
- ``Integral is divergent or slowly convergent in
  inte_singular_gsl::qags().`` [ :cpp:enumerator:`exc_ediverge` ] This
  occurs if the approximations produced by the refinement
  algorithm and the extrapolation algorithm differ by too many
  orders of magnitude.
- ``Exceeded limit of trigonometric table in
  inte_qawo_gsl_sin()::qawo().`` [ :cpp:enumerator:`exc_etable` ] This
  occurs if the maximum level of the table of Chebyshev
  moments is reached.

..
  This section is commented out for now
  
  Multi-dimensional integration routines
  --------------------------------------

  O₂scl reimplements the Cubature library for
  multi-dimensional integration. The h-adaptive and p-adaptive
  integration methods are implemented in :ref:`inte_hcubature
  <inte_hcubature>` and :ref:`inte_pcubature <inte_pcubature>`. See also
  the Monte Carlo integration routines in :ref:`Monte Carlo
  Integration`.

  Multi-dimensional hypercubic integration is performed by
  children of :ref:`inte_multi . Currently in O₂scl, only the 

  General multi-dimensional integration is performed by \ref
  o2scl::inte_gen_comp, the sole descendant of :ref:inte_gen.
  The user is allowed to specify a upper and lower limits which are
  functions of the variables for integrations which have not yet
  been performed, i.e. the n-dimensional integral
  \f[ 
  \int_{x_0=a_0}^{x_0=b_0} f(x_0) \int_{x_1=a_1(x_0)}^{x_1=b_1(x_0)} 
  f(x_0, x_1) ...
  \int_{x_{\mathrm{n}-1}=a_{\mathrm{n}-1}(x_0,x_1,..,x_{\mathrm{n}-2})}^
  {x_{\mathrm{n}-1}=b_{\mathrm{n}-1}(x_0,x_1,..,x_{\mathrm{n}-2})} 
  f(x_0,x_1,...,x_{\mathrm{n-1}})~d x_{\mathrm{n}-1}~...~d x_1~d x_0
  \f]
  Again, one specifies a set of inte objects to apply to
  each variable to be integrated over.

One-dimensional integration example
-----------------------------------

This example computes the integral :math:`\int_{-\infty}^{\infty}
e^{-x^2} ~dx` with :ref:`inte_qagi_gsl <inte_qagi_gsl>`, the integral
:math:`\int_0^{\infty} e^{-x^2} ~dx` with :ref:`inte_qagiu_gsl
<inte_qagiu_gsl>`, the integral :math:`\int_{-\infty}^{0} e^{-x^2}
~dx` with :ref:`inte_qagil_gsl <inte_qagil_gsl>`, and the integral
:math:`\int_0^1 \left[ \sin (2 x) + \frac{1}{2} \right]~dx` with both
:ref:`inte_qag_gsl <inte_qag_gsl>` and :ref:`inte_adapt_cern_tl
<inte_adapt_cern_tl>`, and compares the computed results with the exact
results.

.. literalinclude:: ../../../examples/ex_inte.cpp
   :language: c++		    
   :start-after: sphinx-example-start

.. 
  Multi-dimensional integration example

  This example computes the integral :math:`\int_{0}^{1} \int_{0}^{1}
  \int_{0}^{1} \sqrt{x^3+y^3+z^3+x y^2 z}~dx~dy~dz` with \ref
  o2scl::inte_multi_comp .

  \dontinclude ex_minte.cpp
  \skip Example:
  \until End of example

Gauss-Kronrod integration coefficients
--------------------------------------

:ref:`Top <Physical Constants>`

.. doxygennamespace:: o2scl_inte_gk_coeffs

Non-adaptive quadrature integration coefficients
------------------------------------------------

:ref:`Top <Physical Constants>`
		      
.. doxygennamespace:: o2scl_inte_qng_coeffs

