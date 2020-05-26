:ref:`O2scl <o2scl>`

Least-Squares Fitting
=====================

Linear least-squares fitting is performed by :ref:`fit_linear
<fit_linear>`. Non-linear least-squares Fitting is performed by
descendants of :ref:`fit_base <fit_base>` and fitting functions can be
specified using :ref:`fit_funct <fit_funct>`. The GSL fitting routines
(scaled and unscaled) are implemented in :ref:`fit_nonlin
<fit_nonlin>`. A generic fitting routine using a minimizer object
specified as a child of :ref:`mmin_base <mmin_base>` is implemented in
:ref:`fit_min <fit_min>`. When the :ref:`mmin_base <mmin_base>` object
is (for example) a :ref:`anneal_base <anneal_base>` object,
:ref:`fit_min <fit_min>` can avoid local minima which can occur when
fitting noisy data.

