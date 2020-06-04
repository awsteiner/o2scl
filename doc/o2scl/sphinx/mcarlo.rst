Monte Carlo Integration
=======================

:ref:`O2scl <o2scl>`

.. contents:: 

Monte Carlo integration is performed by descendants of
:ref:`mcarlo <mcarlo>` (:ref:`mcarlo_plain <mcarlo_plain>`,
:ref:`mcarlo_miser <mcarlo_miser>`, and :ref:`mcarlo_vegas
<mcarlo_vegas>`. These routines are generally superior to the direct
methods for integrals over regions with large numbers of spatial
dimensions.

Monte Carlo integration example
-------------------------------

This example computes the integral

.. math::

   \int_0^{\pi} \int_0^{\pi} \int_0^{\pi} \frac{1}{\pi^3}
   \left(1-\cos x \cos y \cos z\right)^{-1}

and compares the result with the exact value, 1.3932039296.

.. literalinclude:: ../../../examples/ex_mcarlo.cpp

Analysis of results from numerical simulations
----------------------------------------------
    
The base :ref:`expval_base <expval_base>` and its children form a set
of classes useful for recording the outputs of several iterations of a
numerical simulation, and summarizing the average, standard deviation,
and error in the average over all iterations. After each iteration,
some measurement is performed, and that measurement is added to the
class with an ``add()`` functions (e.g.
:cpp:func:`expval_scalar::add()` ).
    
Autocorrelations are common in numerical simulations, and these
classes allow the user to assess the impact of correlations by
breaking the measurements up into blocks.
    
Blocks are filled one by one, moving to the next block when the
requested of measurements (``n_per_block``) in the previous
block has been provided. If more measurements are given than the
number of blocks times ``n_per_block``, then the number of
filled blocks is halved, each block is filled with twice as much
data, and ``n_per_block`` is doubled. The blocks are always
combined in such a way as to preserve the original ordering, i.e.
the first block is combined with the second, the third with the
fourth, and so on. Subsequent measurements are added to newly
empty blocks.
    
The constructor for these classes need as input the number of
blocks to record (at least one) and the number of measurements per
block (also at least one). If either of these is specified to be
zero then the error handler is called.
    
In order to properly handle autocorrelations, The number of
measurements over which autocorrelation occurs should be less
than the number of measurements per block. 
    
Current averages are always reported as long as at least one
measurement has been added, though not all data is always used if
there are incomplete blocks or the blocks have not yet been filled
evenly. Two functions ``current_avg()`` and ``current_avg_stats()``
are provided in children. The latter provides the number of blocks and
number of measurements per block used in the currently reported
statistics. If only one (full or partially full) block is available,
the standard deviation and uncertainty on the average are reported as
zero. If no data is available, then calls to ``current_avg()`` will
call the error handler. See, e.g.
:cpp:func:`o2scl::expval_scalar::current_avg()` and
:cpp:func:`o2scl::expval_scalar::current_avg_stats()`.

Many of the :ref:`expval_base <expval_base>` children also support HDF
I/O.

Note that the traditional Monte Carlo integration routines
(:ref:`mcarlo_plain <mcarlo_plain>`, :ref:`mcarlo_miser
<mcarlo_miser>`, and :ref:`mcarlo_vegas <mcarlo_vegas>`) have no
significant autocorrelations by construction, so long as the random
number generator does not have any.
    
For Markov chain Monte Carlo, see the discussion in :ref:`Markov chain
Monte Carlo`.
