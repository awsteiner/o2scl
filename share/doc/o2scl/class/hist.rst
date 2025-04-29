Histograms
==========

:ref:`O2scl <o2scl>`

One- and two-dimensional histograms are implemented in 
:ref:`hist <hist>` and :ref:`hist_2d <hist_2d>`.
    
O₂scl histograms are made up of bins, where the size of the
histogram is equal to the number of bins. Each bin has a upper
boundary and a lower boundary, also called edges. O₂scl
histograms require that the upper edge of each bin is equal to the
lower edge of the following bin. Each bin stores one value or
"weight". Thus if a histogram has size ``n``, it has ``n`` bins but
``n+1`` unique edges.

By convention, most functions which take a bin index as an
argument have an extra ``"_i"`` suffix to distinguish them from
functions which take a floating-point value to be binned as their
argument.

Empty histograms have zero size. Also, the histogram weights are
automatically initialized to zero when memory is allocated for them.
By default, histograms do not allow one to add weights corresponding
to values of the independent variable which are less than the smallest
bin or greater than the largest bin.

Each bin also has a representative coordinate (or set of
coordinates), typically corresponding to some value between the
lower and upper edges which represents the independent variable(s)
for each bin. These bin representative values are automatically
created and can be used by the function evaluation and
interpolation functions for some particular interpolation type. By
default, these representative values are taken to be the midpoint
of each bin, but this option is configurable and the
representative values may be set by the user for each individual
bin.

Histograms can be read and written to HDF files using the
``hdf_output`` and ``hdf_input`` functions in \ref o2scl_hdf.

.. warning::

   The bin edges may be either increasing or decreasing, but should be
   monotonic. Also, no two adjacent bin edges should be equal. The
   histogram classes do not exhaustively check that this is the case.

.. todo::

   Create a 1- and 2-D histogram example.
