File I/O with HDF5
==================

:ref:`O2scl <o2scl>`

The class :ref:`hdf_file <hdf_file>` facilitates I/O of data to hdf
files. This class implements a simple way to I/O basic data types and
O\ :sub:`2`\ scl data types. All files created by the 
:cpp:class:`o2scl_hdf::hdf_file` class are normal HDF5 files, and can be
manipulated in the usual way, for example with ``h5dump``
command-line tool. Users can easily mix code which performs I/O with
:cpp:class:`o2scl_hdf::hdf_file` and other O\ :sub:`2`\ scl functions with
their own HDF code. The sole caveat is that O\ :sub:`2`\ scl cannot
parse generic HDF5 files, so that HDF files which contain data not
output through O\ :sub:`2`\ scl cannot always be read by O\ :sub:`2`\
scl.

Objects are stored by refering to their dataset name. I/O for basic
objects is provided directly in the :cpp:class`o2scl_hdf::hdf_file`
class. Some classes also provide their own I/O functions named
``hdf_output`` and ``hdf_input`` based on the
:cpp:class`o2scl_hdf::hdf_file` class. Some of the current classes
which provide I/O are :ref:`table <table>`, :ref:`table_units
<table_units>`, :ref:`table3d <table3d>`, :ref:`tensor_grid
<tensor_grid>`, :ref:`hist <hist>`, and :ref:`hist_2d <hist_2d>`.
    
O\ :sub:`2`\ scl formats complicated data types for HDF I/O by
combining basic data into groups. For that reason, one cannot use O\
:sub:`2`\ scl to read or write HDF files where groups have the same
name as a dataset in the current HDF id. All O\ :sub:`2`\ scl groups
in HDF files come with a fixed-length string named
``o2scl_type``, which refers to the type of object which has been
written to the HDF file as a group.

.. note:: Vector I/O from HDF5 files can be done directly only if the
	  vector object provides a pointer to a contiguous chunk of
	  memory. This works for ``std::vector`` objects, because the
	  C++ standard guarantees this. It is not necessarily possible
	  for uBlas vector objects (nor desirable for vectors built
	  upon slices of matrices or tensors), and thus HDF5 I/O with
	  uBlas vectors or matrices requires an extra copy. See also
	  the discussion \ref vec_io_cont_subsect in the User's guide.

.. note:: There are some current limitations regarding the matching of
	  error handling policies between O\ :sub:`2`\ scl and the HDF
	  library. HDF functions do not always call the O\ :sub:`2`\
	  scl error handler and thus do not always throw O\ :sub:`2`\
	  scl exceptions.
    
\future Create an HDF file I/O example

