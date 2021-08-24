Function Objects
================

:ref:`O2scl <o2scl>`

Function Object Contents
------------------------

- :ref:`Lambda functions and std::mem_fn <lambda_func>`
- :ref:`First function object example`
- :ref:`General comments about function objects`
- :ref:`Second function object example`
- :ref:`Function typedefs`

Lambda functions and ``std::mem_fn``
------------------------------------

.. _lambda_func:

Functions are passed to numerical routines using template-based
function classes, sometimes called "functors". O\ :sub:`2`\ scl
classes which accept functions as parameters generally default to
types built upon ``std::function``. If the
user would like to use Boost function objects instead, these may
also be used, simply by specifying the Boost function type in the
template parameter. 

Some template aliases are defined to save typing of the function
types, e.g.

- :ref:`funct <funct>` : One function of one variable (used in 
  one-dimensional solver and minimizer classes, derivative classes,
  integration classes, etc.)
- :ref:`funct_ld <funct_ld>` : One function of one variable using
  long double 
- :ref:`multi_funct <multi_funct>` : One function of several variables (used
  in minimizer and integration classes)
- :ref:`mm_funct <mm_funct>`: :math:`n` functions of :math:`n`
  variables (used in solver classes)
- :ref:`grad_funct <grad_funct>`: gradient function for minimizer classes
- :ref:`ode_funct <ode_funct>`: :math:`n` derivatives as a function of
  :math:`n`
  function values and the value of the independent variable
- :ref:`ode_jac_funct <ode_jac_funct>`: Jacobian function for ODE classes
- :ref:`ode_it_funct <ode_it_funct>`: Function to specify ODEs for
  iterative solution
- :ref:`jac_funct <jac_funct>` : Jacobian function for solver and
  fitting classes
- :ref:`fit_funct <fit_funct>`: Fit function
- :ref:`ool_hfunct <ool_hfunct>`: Hessian matrix function for constrained
  minimization

First function object example
-----------------------------

The example below demonstrates how C++11 function objects can be used
with the :ref:`root_brent_gsl <root_brent_gsl>` solver.

.. literalinclude:: ../../../examples/ex_lambda.cpp
   :language: c++		    
   :start-after: sphinx-example-start

General comments about function objects
---------------------------------------

The C++ standard library functors employ copy construction at
various types, so one must be careful about the types involved in
creating the functor. Generally, all classes should have
constructors and structs should be avoided because they can cause
difficulties with default copy construction.

There is a small overhead associated with the indirection: a "user
class" accesses the function class which then calls function which
was specified in the constructor of the function class. In many
problems, the overhead associated with the indirection is small.
Some of this overhead can always be avoided by inheriting directly
from the function class and thus the user class will make a direct
virtual function call. To eliminate the overhead entirely, one can
specify a new type for the template parameter in the user class.

Second function object example
------------------------------

This example shows how to provide functions to O\ :sub:`2`\ scl
classes by solving the equation

.. math::

   \left\{ 1+\frac{1}{p_2} 
   \sin \left[ 50 \left( x-p_1 \right) \right] \right\}
   \tan^{-1} \left[ 4 \left( x-p_1 \right) \right] = 0

Where :math:`p_1 = 0.01` and :math:`p_2 = 1.1`. The parameter 
:math:`p_1` is stored as member data for the class, and the 
parameter :math:`p_2` is an argument to the member function.
    
The image below shows how the solver progresses to the 
solution of the example function.

.. image:: ../../../examples/plot/ex_fptr_plot.png
   :width: 60%	   
   :alt: alt text

.. literalinclude:: ../../../examples/ex_fptr.cpp
   :language: c++		    
   :start-after: sphinx-example-start

Function typedefs
-----------------

.. _funct:

.. doxygentypedef:: funct

.. _funct_ld:

.. doxygentypedef:: funct_ld

.. _multi_funct:

.. doxygentypedef:: multi_funct

.. _mm_funct:

.. doxygentypedef:: mm_funct

.. _grad_funct:

.. doxygentypedef:: grad_funct

.. _ode_funct:

.. doxygentypedef:: ode_funct

.. _ode_jac_funct:

.. doxygentypedef:: ode_jac_funct

.. _ode_it_funct:

.. doxygentypedef:: ode_it_funct

.. _jac_funct:

.. doxygentypedef:: jac_funct

.. _fit_funct:

.. doxygentypedef:: fit_funct

.. _ool_hfunct:

.. doxygentypedef:: ool_hfunct

