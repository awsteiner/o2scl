:ref:`O2scl <o2scl>`

Interpolation
=============

Basic interpolation of generic vector types is performed by \ref
o2scl::interp and \ref o2scl::interp_vec. The vector
representing the independent variable must be monotonic, but need
not be equally-spaced. The difference between the two classes is
analogous to the difference between using \c gsl_interp_eval() and
\c gsl_spline_eval() in \c GSL. You can create a \ref
o2scl::interp object and use it to interpolate among any
pair of chosen vectors. For example, cubic spline interpolation
with natural boundary conditions
\code
boost::numeric::ublas::vector<double> x(20), y(20);
// fill x and y with data
o2scl::interp<> oi(itp_cspline);
double y_half=oi.eval(0.5,20,x,y);
\endcode
Alternatively, you can create a \ref o2scl::interp_vec
object which can be optimized for a pair of vectors that you
specify in advance (now using linear interpolation instead)
\code
boost::numeric::ublas::vector<double> x(20), y(20);
// fill x and y with data
o2scl::interp_vec<> oi(20,x,y,itp_linear);
double y_half=oi.eval(0.5);
\endcode

These interpolation classes require that the vector \c x is either
monotonically increasing or monotonically decreasing, but these
classes do not verify that this is the case. If the vectors or
not monotonic, then the interpolation functions are not defined.
These classes give identical results to the GSL interpolation
routines when the vector is monotonically increasing.

These interpolation classes will extrapolate to regions outside
those defined by the user-specified vectors and will not warn you
when they do this (this is not the same behavior as in GSL).

The different interpolation types are defined in \ref interp.h
- \ref o2scl::itp_linear - Linear interpolation
- \ref o2scl::itp_cspline - Cubic spline interpolation with natural
boundary conditions
- \ref o2scl::itp_cspline_peri - Cubic spline interpolation with periodic
boundary conditions
- \ref o2scl::itp_akima - Akima interpolation with natural
boundary conditions
- \ref o2scl::itp_akima_peri - Akima interpolation with periodic
boundary conditions
- \ref o2scl::itp_monotonic - Monotonicity-preserving interpolation
of Fritsch and Carlson (unfinished)
- \ref o2scl::itp_steffen - Monotonicity-preserving method of
\ref Steffen90
- \ref o2scl::itp_nearest_neigh - Nearest neighbor interpolation

Integrals are always computed assuming that if the limits are
ordered so that if the upper limit appears earlier in the array \c
x in comparison to the lower limit, that the value of the integral
has the opposite sign than if the upper limit appears later in the
array \c x .

The classes \ref o2scl::interp and \ref
o2scl::interp_vec are based on the lower-level interpolation
classes of type \ref o2scl::interp_base. Also, the interpolation
classes based on \ref o2scl::interp_base and also the class \ref
o2scl::interp_vec also have defined a function
<tt>operator()</tt> which also returns the result of the
interpolation.

Two specializations for C-style arrays of double-precision numbers
are provided in \ref o2scl::interp_array and \ref
o2scl::interp_array_vec. 
    
An experimental class for one-dimensional kriging is also 
provided in \ref o2scl::interp_krige .
    
Lookup and binary search
------------------------

The classes \ref o2scl::search_vec and \ref o2scl::search_vec_ext
contain searching functions for generic vector types which contain
monotonic (either increasing or decreasing) data. It is \ref
o2scl::search_vec which is used internally by the interpolation
classes to perform cached binary searching. These classes also
allow one to to exahaustively search for the index of an element
in a vector without regard to any kind of ordering, e.g. \ref
o2scl::search_vec::ordered_lookup() .

Two and higher-dimensional interpolation
----------------------------------------

Support for multi-dimensional interpolation is documented in
\ref tintp_section.

Inverse interpolation and related functions
-------------------------------------------

The equivalent to "inverse" linear interpolation, which computes
all the abcissae which have a fixed value of the ordinate, is
implemented in the template function \ref
o2scl::vector_find_level() which is documented in \ref interp.h .
This function together with \ref
o2scl::vector_invert_enclosed_sum() can be used to determine
confidence limits surrounding the peak of a 1-dimensional data set
using linear interpolation. To count level crossings in a
function, use \ref o2scl::vector_level_count(). The function \ref
o2scl::vector_integ_interp() uses interpolation to compute
the integral defined by a set of vectors, and the function \ref
o2scl::vector_region_fracint() finds the set of regions which gives
a fraction of the integral reported by \ref
o2scl::vector_integ_interp().

Derivatives and integrals on a fixed grid
-----------------------------------------
    
If the indepedent variable is represented by a uniform
(equally-spaced) then the functions in \ref vector_derint.h
can provide faster (and occasionally more accurate) results.

