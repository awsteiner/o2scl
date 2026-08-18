File I/O with HDF5
==================

:ref:`O2scl <o2scl>`

The class :ref:`hdf_file <hdf_file>` facilitates I/O of data to hdf
files. This class implements a simple way to read and write basic data
types and O₂scl data types. All files created by the
:cpp:class:`o2scl_hdf::hdf_file` class are normal HDF5 files, and can
be manipulated in the usual way, for example with ``h5dump``
command-line tool. Users can easily mix code which performs I/O with
:cpp:class:`o2scl_hdf::hdf_file` and other O₂scl functions with their
own HDF code. The caveat is that O₂scl has limited support for parsing
generic HDF5 files, so that HDF files which contain data not output
through O₂scl cannot always be read by O₂scl.

HDF5 files which are written by O₂scl can be viewed and manipulated on
the command line using the ``acol`` utility. See :ref:`The acol
Command Line Utility` for more information.

Objects are stored by referring to their dataset name. I/O for basic
objects is provided directly in the :ref:`hdf_file <hdf_file>` class.
Some classes also provide their own I/O functions named
``hdf_output``, ``hdf_input``, and ``hdf_input_n`` based on the
:ref:`hdf_file <hdf_file>` class. Some of the current classes which
provide I/O are :ref:`hist <hist>`, :ref:`hist_2d <hist_2d>`,
:ref:`kde_python <kde_python>`, :ref:`prob_dens_mdim_amr
<prob_dens_mdim_amr>`, :ref:`prob_dens_mdim_gaussian
<prob_dens_mdim_gaussian>`, :ref:`prob_dens_mdim_gmm
<prob_dens_mdim_gmm>`, :ref:`table <table>`, :ref:`table_units
<table_units>`, :ref:`table3d <table3d>`, :ref:`tensor_grid
<tensor_grid>`, and :ref:`uniform_grid <uniform_grid>`.
   
O₂scl formats complicated data types for HDF I/O by combining basic
data into groups. For that reason, one cannot use O₂scl to read or
write HDF files where groups have the same name as a dataset in the
current HDF id. All O₂scl groups in HDF files come with a fixed-length
string named ``o2scl_type``, which refers to the type of object which
has been written to the HDF file as a group. Once an object of a
particular type is written to an HDF5 file, O₂scl does not support
overwriting an object with the same name but a different type to
the same HDF5 file. O₂scl does support overwriting an object with
different data if it has the same type.

.. note:: Vector I/O from HDF5 files can be performed directly only if
	  the vector object provides a pointer to a contiguous chunk
	  of memory. This works for ``std::vector`` objects, because
	  the C++ standard guarantees this. It is not necessarily
	  possible for uBlas vector objects (nor desirable for vectors
	  built upon slices of matrices or tensors), and thus HDF5 I/O
	  with uBlas vectors or matrices requires an extra copy.

.. note:: There are some current limitations regarding the matching of
	  error handling policies between O₂scl and the HDF
	  library. HDF functions do not always call the O₂scl error
          handler and thus do not always throw O₂scl exceptions.
    
.. todo:: (Future) Create an HDF file I/O example

