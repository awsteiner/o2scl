Interpolation
=============

:ref:`O2scl <o2scl>`

Interpolation Contents
----------------------

- :ref:`Interpolation Introduction`
- :ref:`Lookup and binary search`

Interpolation Introduction
--------------------------

Basic interpolation of generic vector types is performed by
:ref:`interp <interp>` and :ref:`interp_vec <interp_vec>`. The vector
representing the independent variable must be monotonic, but need
not be equally-spaced. The difference between the two classes is
analogous to the difference between using ``gsl_interp_eval()`` and
``gsl_spline_eval()`` in GSL. You can create a :ref:`interp <interp>`
object and use it to interpolate among any
pair of chosen vectors. For example, cubic spline interpolation
with natural boundary conditions::

  boost::numeric::ublas::vector<double> x(20), y(20);
  // fill x and y with data
  o2scl::interp<> oi(itp_cspline);
  double y_half=oi.eval(0.5,20,x,y);

Alternatively, you can create a :ref:`interp_vec <interp_vec>`
object which can be optimized for a pair of vectors that you
specify in advance (now using linear interpolation instead)::

  boost::numeric::ublas::vector<double> x(20), y(20);
  // fill x and y with data
  o2scl::interp_vec<> oi(20,x,y,itp_linear);
  double y_half=oi.eval(0.5);

These interpolation classes require that the vector ``x`` is either
monotonically increasing or monotonically decreasing, but these
classes do not verify that this is the case. If the vectors or
not monotonic, then the interpolation functions are not defined.
These classes give identical results to the GSL interpolation
routines when the vector is monotonically increasing.

These interpolation classes will extrapolate to regions outside
those defined by the user-specified vectors and will not warn you
when they do this (this is not the same behavior as in GSL).

The different interpolation types are defined in ``src/base/interp.h``

.. doxygenenumvalue:: itp_linear

.. doxygenenumvalue:: itp_cspline

.. doxygenenumvalue:: itp_cspline_peri

.. doxygenenumvalue:: itp_akima

.. doxygenenumvalue:: itp_akima_peri

.. doxygenenumvalue:: itp_monotonic

.. doxygenenumvalue:: itp_steffen

.. doxygenenumvalue:: itp_nearest_neigh

Integrals are always computed assuming that if the limits are
ordered so that if the upper limit appears earlier in the array
``x`` in comparison to the lower limit, that the value of the integral
has the opposite sign than if the upper limit appears later in the
array ``x``.

The classes :ref:`interp <interp>` and :ref:`interp_vec <interp_vec>`
are based on the lower-level interpolation classes of type
:ref:`interp_base <interp_base>`. Also, the interpolation classes
based on :ref:`interp_base <interp_base>` and also the class
:ref:`interp_vec <interp_vec>` also have defined a function
``operator()`` which also returns the result of the interpolation.

Two specializations for C-style arrays of double-precision numbers are
provided in :ref:`interp_array <interp_array>` and
:ref:`interp_array_vec <interp_array_vec>`.

An experimental class for one-dimensional kriging is also provided in
:ref:`interp_krige <interp_krige>`.
    
Lookup and binary search
------------------------

The classes :ref:`search_vec <search_vec>` and :ref:`search_vec_ext
<search_vec_ext>` contain searching functions for generic vector types
which contain monotonic (either increasing or decreasing) data. It is
:ref:`search_vec <search_vec>` which is used internally by the
interpolation classes to perform cached binary searching. These
classes also allow one to to exahaustively search for the index of an
element in a vector without regard to any kind of ordering, e.g. 
:cpp:func:`o2scl::search_vec::ordered_lookup()`.

Interpolation example
---------------------

.. literalinclude:: ../../../examples/ex_interp.cpp
   :language: c++		    
   :start-after: sphinx-example-start

.. image ../../../examples/plot/ex_fptr_plot.png
   width: 60%	   
   alt: alt text

.. todo:: Fix the interpolation plot for this example.
   
Two and higher-dimensional interpolation
----------------------------------------

Support for multi-dimensional interpolation is documented in
:ref:`Higher-dimensional Interpolation`.

Inverse interpolation and related functions
-------------------------------------------

The equivalent to "inverse" linear interpolation, which computes all
the abcissae which have a fixed value of the ordinate, is implemented
in the template function :cpp:func:`vector_find_level()`. This
function together with \ref :cpp:func:`vector_invert_enclosed_sum()`
can be used to determine confidence limits surrounding the peak of a
1-dimensional data set using linear interpolation. To count level
crossings in a function, use :cpp:func:`vector_level_count()`. The
function \ref :cpp:func:`vector_integ_interp()` uses interpolation to
compute the integral defined by a set of vectors, and the function
:cpp:func:`o2scl::vector_region_fracint()` finds the set of regions
which gives a fraction of the integral reported by \ref
:cpp:func:`o2scl::vector_integ_interp()`.

Derivatives and integrals on a fixed grid
-----------------------------------------
    
If the indepedent variable is represented by a uniform
(equally-spaced) then the functions in ``src/deriv/vector_derint.h``
can provide faster (and occasionally more accurate) results. See:

- :cpp:func:`o2scl::vector_deriv_fivept()`
- :cpp:func:`o2scl::vector_deriv_fivept_tap()`
- :cpp:func:`o2scl::vector_deriv_interp()`
- :cpp:func:`o2scl::vector_deriv_threept()`
- :cpp:func:`o2scl::vector_deriv_threept_tap()`
- :cpp:func:`o2scl::vector_integ_durand()`
- :cpp:func:`o2scl::vector_integ_extended4()`
- :cpp:func:`o2scl::vector_integ_extended8()`
- :cpp:func:`o2scl::vector_integ_threept()`
- :cpp:func:`o2scl::vector_integ_trap()`


