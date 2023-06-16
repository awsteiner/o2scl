Differentiation
===============

:ref:`O2scl <o2scl>`

Differentiation contents
------------------------

- :ref:`Differentiation introduction`
- :ref:`Differentiation example`  

Differentiation introduction
----------------------------

A GSL-based routine for derivatives is provided in :ref:`deriv_gsl
<deriv_gsl>`, and a CERNLIB-based routine is given in :ref:`deriv_cern
<deriv_cern>`. These classes allow one to calculate either first,
second, and third derivatives. For functions which are tabulated over
equally-spaced abscissas, the class :ref:`deriv_eqi <deriv_eqi>` is
provided which applies the formulas from Abramowitz and Stegun at a
specified order. The class :ref:`deriv_cern <deriv_cern>` is slower
and sometimes more accurate, but also fails more often than
:ref:`deriv_gsl <deriv_gsl>`, which never calls the error handler. The
GSL derivative class :ref:`deriv_gsl <deriv_gsl>` supports numerical
derivatives of functions which operate on multiprecision numbers (see
also :ref:`Multiprecision Support`). The class :ref:`deriv_multip_gsl
<deriv_multip_gsl>` uses multiprecision to automatically compute a
derivative to within a requested precision.

The estimatation of second and third derivatives by :ref:`deriv_gsl
<deriv_gsl>` and :ref:`deriv_cern <deriv_cern>` can be particularly
troublesome if the function is not sufficiently smooth. Error
estimation is not provided for second and third derivatives. If one
has more information about the function, then second and third
derivatives are often better computed by fitting to a model and then
taking the second or third derivative of the model instead.

Differentiation example
-----------------------

This example computes first and second derivatives of

.. math::
   
   y(x) = \sin (2 x) + \frac{1}{2}

with both :ref:`deriv_gsl <deriv_gsl>` and :ref:`deriv_cern <deriv_cern>`.

.. literalinclude:: ../../../examples/ex_deriv.cpp
   :language: c++		    
   :start-after: sphinx-example-start
