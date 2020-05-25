Differentiation
===============

:ref:`O2scl <o2scl>`

Differentiation is performed by descendants of :ref:`deriv_base
<deriv_base>`. These allow one to calculate either first, second, and
third derivatives. A GSL-based routine is used in :ref:`deriv_gsl
<deriv_gsl>`, and the CERNLIB routine is used in :ref:`deriv_cern
<deriv_cern>`. For functions which are tabulated over equally-spaced
abscissas, the class :ref:`deriv_eqi <deriv_eqi>` is provided which
applies the formulas from Abramowitz and Stegun at a specified order.
The class :ref:`deriv_cern <deriv_cern>` is slower and sometimes more
accurate, but also fails more often than :ref:`deriv_gsl <deriv_gsl>`,
which never calls the error handler.

.. warning::
   For :ref:`deriv_gsl <deriv_gsl>` and :ref:`deriv_cern <deriv_cern>`
   the second and third derivatives are calculated by naive repeated
   application of the code for the first derivative and can be
   particularly troublesome if the function is not sufficiently
   smooth. Error estimation is not provided for second and third
   derivatives.

Differentiation example
-----------------------

This example computes first and second derivatives of

.. math::
   
   y(x) = \sin (2 x) + \frac{1}{2}

with both :ref:`deriv_gsl <deriv_gsl>` and :ref:`deriv_cern <deriv_cern>`.

.. literalinclude:: ../../../examples/ex_deriv.cpp
