.. _o2scl:

Object-oriented Scientific Computing Library: v0.927a1
======================================================

O\ :sub:`2`\ scl is a C++ library for object-oriented scientific
computing. This is a beta version. The library should install and test
successfully, and most of the classes are ready for production use.
Some of the interfaces may change slightly in future versions. There
are a few classes which are more experimental, and this is stated at
the top of the documentation for these classes.

There are sub-libraries for handling thermodynamics of particles and
nuclei, `O2scl_part <../part/html/index.html>`_ and equations of state
of dense matter, `O2scl_eos <../eos/html/index.html>`_ .

Gallery
=======

.. |fptr| image:: ../../../examples/plot/ex_fptr_plot.png
		  
.. |anneal| image:: ../../../examples/plot/ex_anneal_plot.png
		    
.. |ode| image:: ../../../examples/plot/ex_ode_bessel.png
		 
.. |nucprof| image:: ../../../examples/plot/ex_nuc_prof.png
		     
.. |nucmass| image:: ../../../examples/plot/ex_nucmass_dz96.png

.. csv-table::

   Root finding,Minimization,ODEs,Nuclear Structure,Nuclear masses
   |fptr|,|anneal|,|ode|,|nucprof|,|nucmass|

User's Guide
============
   
.. toctree::
   :maxdepth: 1

   download
   install
   usage
   vecmat
   linalg
   interp
   const
   funct
   table
   string
   diff
   inte
   poly
   solve
   min
   mcarlo
   mcmc
   conmin
   anneal
   hist
   fit
   ode
   rng
   multi_interp
   algebraic
   cheb
   unitconv
   multip
   hdf
   acol
   python
   para
   other
   lset
   yanic
   design
   dev_guide
   license
   related
   ack
   ref
   class_list
   function_list
   todos

* :ref:`genindex`
* :ref:`search`

Test sphinx references
----------------------

Class reference: :ref:`table <table>` or
:cpp:class:`o2scl::table` or :cpp:class:`table <o2scl::table>` .

Function reference: :cpp:func:`o2scl::tensor_grid::set_grid`
or :cpp:func:`tensor_grid::set_grid() <o2scl::tensor_grid::set_grid>` 

