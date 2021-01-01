Equation Solving
================

:ref:`O2scl <o2scl>`

Equation Solving Contents
-------------------------

- :ref:`One-dimensional solvers`
- :ref:`Multi-dimensional solvers`
- :ref:`Multi-dimensional solver example`

One-dimensional solvers
-----------------------

Solution of one equation in one variable is accomplished by
children of the class :ref:`root <root>`. 

For one-dimensional solving, if the root is bracketed, use
:ref:`root_bkt_cern <root_bkt_cern>` or :ref:`root_brent_gsl
<root_brent_gsl>`. The :ref:`root_bkt_cern <root_bkt_cern>` class is
typically faster (for the same accuracy) than :ref:`root_brent_gsl
<root_brent_gsl>`. If a relatively fast derivative is available, use
:ref:`root_stef <root_stef>`. If neither a bracket nor a derivative is
available, you can use :ref:`root_cern <root_cern>`.

The :ref:`root <root>` base class provides the structure for three
different solving methods:

- :cpp:func:`o2scl::root::solve()` which solves a function given an
  initial guess ``x``
- :cpp:func:`o2scl::root::solve_bkt()` which solves a function given a
  solution bracketed between ``x1`` and ``x2``. The values of the
  function at ``x1`` and ``x2`` should have different signs.
- :cpp:func:`o2scl::root::solve_de()` which solves a function given an
  initial guess ``x`` and the function's derivative.

There is an example using the one-dimensional solver at
:ref:`Second function object example`.

The :ref:`root <root>` base class also contains the relative tolerance
(:cpp:var:`o2scl::root::tol_rel`), absolute tolerance
(:cpp:var:`o2scl::root::tol_abs`), the number of iterations
(:cpp:var:`o2scl::root::ntrial`), the verbosity parameter
(:cpp:var:`o2scl::root::verbose`), and the number of iterations in the
last solve (:cpp:var:`root::last_ntrial`).

If not all of these three functions are overloaded, then the source
code in the :ref:`root <root>` base class is designed to try to
automatically provide the solution using the remaining functions. Most
of the one-dimensional solving routines, in their original form, are
written in the second or third form above. For example,
:ref:`root_brent_gsl <root_brent_gsl>` is originally a bracketing
routine of the form :cpp:func:`o2scl::root::solve_bkt()`, but calls to
either :cpp:func:`o2scl::root::solve()` or
:cpp:func:`o2scl::root::solve_de()` will attempt to automatically
bracket the function given the initial guess that is provided. Of
course, it is frequently most efficient to use the solver in the way
it was intended.

Multi-dimensional solvers
-------------------------

Solution of more than one equation is accomplished by descendants
of the class :ref:`mroot <mroot>`. The higher-level interface is
provided by the function :cpp:func:`o2scl::mroot::msolve()`.

For multi-dimensional solving, you can use either :ref:`mroot_cern
<mroot_cern>` or :ref:`mroot_hybrids <mroot_hybrids>`. While
:ref:`mroot_cern <mroot_cern>` cannot utilize user-supplied
derivatives, :ref:`mroot_hybrids <mroot_hybrids>` can use
user-supplied derivative information (as in the GSL ``hybridsj``
method) using the function
:cpp:func:`o2scl::mroot_hybrids::msolve_de()` .

A specialization of :ref:`mroot_hybrids <mroot_hybrids>` for Armadillo
is given in :ref:`mroot_hybrids_arma_qr_econ
<mroot_hybrids_arma_qr_econ>` where the QR decomposition used in the
solver is performed by the Armadillo library. A similar specialization
for Eigen is in :ref:`mroot_hybrids_eigen <mroot_hybrids_eigen>` . These
specializations will be faster than when the number of variables is
sufficiently large.

Multi-dimensional solver example
--------------------------------

This demonstrates several ways of using the multi-dimensional
solvers to solve the equations

.. math::
   
   \begin{eqnarray}
   \sin \left( x_0 - \frac{1}{4} \right) &=& 0 \nonumber \\
   \sin \left( x_1 - \frac{1}{5} \right) &=& 0
   \end{eqnarray}

.. literalinclude:: ../../../examples/ex_mroot.cpp
   :language: c++		    
   :start-after: sphinx-example-start
