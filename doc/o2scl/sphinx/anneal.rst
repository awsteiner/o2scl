Simulated Annealing
===================

:ref:`O2scl <o2scl>`

Simulated annealing contents
----------------------------

- :ref:`Simulated annealing introduction`
- :ref:`Simulated annealing example`
     
Simulated annealing introduction
--------------------------------

Minimization by simulated annealing is performed by a modified
version of the GSL algorithm in the :ref:`anneal_gsl <anneal_gsl>`
class. A version which includes OpenMP and MPI parallelization
is provided in :ref:`anneal_para <anneal_para>`.

Simulated annealing example
---------------------------

This example minimizes the function

.. math::

   f(x,y) = J_0(x-2) J_0(y+3)

over :math:`(x,y)` where :math:`J_0(x)` is the Bessel function
given in ``gsl_sf_bessel_J0``. The initial guess at :math:`(9,9)`
is far away from the global minimum. 

The plot below plots the function above, the initial guess, and
the minimum obtained by the example program.

.. image:: ../../../examples/plot/ex_anneal_plot.png
   :width: 60%	   
   :alt: A density and contour plot of the product of Bessel
         functions, the initial point, and the global minimum.

.. literalinclude:: ../../../examples/ex_anneal.cpp
   :language: c++		    
   :start-after: sphinx-example-start

