Constrained Minimization
========================

:ref:`O2scl <o2scl>`

Constrained minimization contents
---------------------------------

- :ref:`Constrained minimization introduction`
- :ref:`Constrained minimization example`
     
Constrained minimization introduction
-------------------------------------

.. note::

   The functionality provided by these legacy classes is probably
   better found in either https://github.com/stevengj/nlopt or
   https://petsc.org/main/docs/manual/tao/ .

O₂scl reimplements the Open Optimization Library (OOL) available at
http://ool.sourceforge.net. The associated classes allow constrained
minimization when the constraint can be expressed as a hyper-cubic
constraint on all of the independent variables. The routines have been
rewritten and reformatted for C++ in order to facilitate the use of
member functions and user-defined vector types as arguments. The base
class is :ref:`mmin_constr <mmin_constr>` and there are two different
constrained minimzation algorithms implemented in
:ref:`mmin_constr_pgrad <mmin_constr_pgrad>`, :ref:`mmin_constr_spg
<mmin_constr_spg>`. (The :ref:`mmin_constr_gencan
<mmin_constr_gencan>` minimizer is not yet finished). The O₂scl
implementation should be essentially identical to the most recently
released version of OOL.

The constrained minimization classes operate in a similar way to the
other multi-dimensional minimization classes (which are derived from
:ref:`mmin_base <mmin_base>`). The constraints are specified with the
function::

  mmin_constr::set_constraints(size_t nc, vec_t &lower, 

and the minimization can be performed by calling either
:cpp:func:`o2scl::mmin_base::mmin()` or
:cpp:func:`o2scl::mmin_base::mmin_de()` (if the gradient is provided
by the user). The method in :ref:`mmin_constr_gencan
<mmin_constr_gencan>` requires a Hessian vector product and the user
can specify this product for the minimization by using
:cpp:func:`o2scl::mmin_constr::mmin_hess()`. The Hessian product
function can be specified as an object of type :ref:`ool_hfunct
<ool_hfunct>` in a similar way to the other function objects in O₂scl.

There are five error codes defined in :ref:`mmin_constr <mmin_constr>`
which are specific to the classes derived from OOL.

The class :ref:`anneal_gsl <anneal_gsl>` can handle some kinds of
constraints by ignoring proposed steps which cause the user-specified
function to return a non-zero value.

Also, a simple way of implementing constraints is to add a function to
the original which increases the value outside of the allowed region.
This can be done with the functions :cpp:func:`o2scl::constraint()`
and :cpp:func:`o2scl::lower_bound()`. There are two analogous
functions, :cpp:func:`o2scl::cont_constraint()` and
:cpp:func:`o2scl::cont_lower_bound()`, which continuous and
differentiable versions. Where possible, it is better to use the
constrained minimization routines described above.

Constrained minimization example
--------------------------------

This example minimizes the function

.. math::

   f(x,y) = \left[x^2 \log(x)+1\right]\left[\sqrt{y}(y-1)+1\right)]

which is undefined for :math:`x<0` and :math:`y<0`. The function is
also minimized by :ref:`mmin_simp2 <mmin_simp2>`, which goes outside
the allowed region where the function is undefined.

.. literalinclude:: ../../../examples/ex_conmin.cpp
   :language: c++		    
   :start-after: sphinx-example-start
