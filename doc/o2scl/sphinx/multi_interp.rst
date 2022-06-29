Higher-dimensional Interpolation
================================

:ref:`O2scl <o2scl>`

Higher-dimensional Interpolation Contents
-----------------------------------------
     
- :ref:`Two-dimensional interpolation`
- :ref:`Multi-dimensional interpolation`
- :ref:`Interpolation on a rectangular grid`
- :ref:`Contour lines`  

Two-dimensional interpolation
-----------------------------

There are two types of two-dimensional interpolation classes, the
first is based on a function defined on a two-dimensional grid
(though the spacings between grid points need not be equal). The
class :ref:`interp2_direct <interp2_direct>` implements bilinear or bicubic
interpolation, and is based on D. Zaslavsky's routines at
https://github.com/diazona/interp2d (licensed under GPLv3).
A slightly slower (but a bit more flexible) alternative is 
successive use of :ref:`interp_base <interp_base>` objects, implemented
in :ref:`interp2_seq <interp2_seq>` . 

If data is arranged without a grid, then :ref:`interp2_neigh
performs nearest-neighbor interpolation. At present, the only way
to compute :ref:`contour <contour>` lines on data which is not defined on a grid
is to use this class or one of the multi-dimensional interpolation
classes described below the data on a grid and then use :ref:`contour
<contour>` afterwards.

.. 
  7/10/19: I removed the reference to interp2_planar because
  it's unstable and I don't recommend using it. 

Multi-dimensional interpolation
-------------------------------

Multi-dimensional interpolation for table defined on a grid is
possible with :ref:`tensor_grid <tensor_grid>`. See the documentation
for :cpp:func:`o2scl::tensor_grid::interpolate()`,
:cpp:func:`o2scl::tensor_grid::interp_linear()` and
:cpp:func:`o2scl::tensor_grid::rearrange_and_copy()`. Also, if you
want to interpolate ``rank-1`` indices to get a vector result, you can
use :cpp:func:`o2scl::tensor_grid::interp_linear_vec()` .

If the data is not on a grid, then inverse distance weighted
interpolation is performed by :ref:`interpm_idw <interpm_idw>`.

An experimental class for multidimensional-dimensional kriging is also 
provided in :ref:`interpm_krige <interpm_krige>` .
    
Interpolation on a rectangular grid
-----------------------------------

.. literalinclude:: ../../../examples/ex_interp2.cpp
   :language: c++		    
   :start-after: sphinx-example-start

This example creates a sample 3 by 3 grid of data with the 
function :math:`\left[ \sin \left( x/10 + 3 y/10 \right) \right]^2`
and performs some interpolations and compares them with the 
exact result.

.. literalinclude:: ../../../examples/ex_interp2.scr

..
  AWS: 6/6/19: I'm commenting this out because interp2_planar is
  unstable and probably not recommended.

  \section ex_interp2_planar_sect Interpolation of randomly spaced points
    
  For example, with 10 random points in the x-y plane with \f$
  -1<x<1 \f$ and \f$ -1<y<1 \f$, the figure contains several
  polygonal regions, each of which represents the set of all points
  in the domain which will be mapped to the same plane in order to
  to approximate the original function.

  \image html ex_planar_plot.png "Planes from interp2_planar class"
  \image latex ex_planar_plot.pdf "Planes from interp2_planar class" width=9cm

Contour lines
-------------

This example generates contour lines of the function

.. math::

   z = f(x,y) = 15 \exp \left[ - \frac{1}{20^2}\left( x-20 \right)^2 
   - \frac{1}{5^2}\left(y-5\right)^2\right] + 
   40 \exp \left[ - \frac{1}{500}\left( x-70 \right)^2 
   - \frac{1}{2^2}\left(y-2\right)^2\right]

.. literalinclude:: ../../../examples/ex_contour.cpp
   :language: c++		    
   :start-after: sphinx-example-start

The figure below shows contour lines in the region :math:`x\in(0,121),
y\in(0,9)`. The data grid is represented by plus signs, and the
associated generated contours. The figure clearly shows the peaks at
:math:`(20,5)` and :math:`(70,2)`.

.. image:: ../../../examples/plot/ex_contour_plot1.png
   :width: 60%	   
   :alt: alt text

The :ref:`contour <contour>` class can also use interpolation to
attempt to refine the data grid. The new contours after a refinement
of a factor of 5 is given in the figure below.
    
.. image:: ../../../examples/plot/ex_contour_plot2.png
   :width: 60%	   
   :alt: alt text

