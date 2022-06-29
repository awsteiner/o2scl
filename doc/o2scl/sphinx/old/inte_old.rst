..

  This section is commented out for now

  - :ref:`Multi-dimensional integration routines`
  
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

