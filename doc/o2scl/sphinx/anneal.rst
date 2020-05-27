:ref:`O2scl <o2scl>`

Simulated Annealing
===================

Minimization by simulated annealing is performed by descendants of
sim_anneal (see :ref:`anneal_gsl <anneal_gsl>`). 

.. 
  Because simulated
  annealing is particularly well-suited to parallelization, a
  multi-threaded minimizer analogous to \ref o2scl::anneal_gsl is
  given in \ref o2scl::anneal_mt. This header-only class uses the
  Boost libraries and requires it for use.

Simulated annealing example
---------------------------

This example minimizes the function

.. math::

   f(x,y) = J_0(x-2) J_0(y+3)

over :math:`(x,y)` where :math:`J_0(x)` is the Bessel function
given in \c gsl_sf_bessel_J0. The initial guess at :math:`(9,9)`
is far away from the global minimum. 

The plot below plots the function above, the initial guess, and
the minimum obtained by the example program.

.. image:: ../../../examples/plot/ex_anneal_plot.png
   :width: 60%	   
   :alt: alt text

.. literalinclude:: ../../../examples/ex_anneal.cpp
   :language: c++		    
   :start-after: sphinx-example-start

