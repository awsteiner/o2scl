Arrays, Vectors, Matrices and Tensors
=====================================

:ref:`O2scl <o2scl>`

.. contents:: 

Because many such vector and matrix objects are defined elsewhere, O\
:sub:`2`\ scl no longer includes native vector and matrix classes.
Internally, O\ :sub:`2`\ scl most often uses Boost uBLAS vector and
matrix objects: ``boost::numeric::ublas::vector<>``,
``boost::numeric::ublas::matrix<>``, and other related class
templates. Many O\ :sub:`2`\ scl routines are templates which are
compatible with a wide range of vector and matrix types. See the
:ref:`Multi-dimensional solver example` which shows how an O\
:sub:`2`\ scl class can be used with Boost, Eigen, or Armadillo
objects.

The O\ :sub:`2`\ scl library uses a standard nomenclature to distinguish a
couple different concepts. The word "array" is always used to
refer to C-style arrays, i.e. ``double[]``. If there are two
dimensions in the array, it is a "two-dimensional array", i.e.
``double[][]`` . The word "vector" is reserved generic
objects with array-like semantics.

In general, there are many vector types (STL, Boost, etc.) and
they can be characterized by whether or not they satisfy certain
"concepts" like ``DefaultConstructible``. O\ :sub:`2`\ scl classes which
operate on vector types are designed to be as flexible as
possible, so that they can be used with almost any vector type.
Eventually, all O\ :sub:`2`\ scl classes with template vector and matrix types
will specify exactly which concepts are required to be satisified,
but this is still in progress.

The word "matrix" is reserved for the a generic object which has
matrix-like semantics and can be accessed using ``operator(,)``. C++
matrix types typically prefer ``operator(,)`` over ``operator[][]``.
This is because ``operator[][]`` implies the creation of a temporary
row object, and it is thus difficult to implement ``operator[]``
without incurring an overhead. Nevertheless, some O\ :sub:`2`\ scl
classes have separate versions which operate on matrix types which are
only accessible with ``operator[][]`` (like two-dimensional arrays).
See the :ref:`Linear Algebra` section of the User's guide for examples
of this.

With ``std::function<>`` and the new lambda function support in C++11,
it is important to notice that ``std::function<double
&(size_t,size_t)>`` is also a matrix type (the ampersand is important
unless the matrix is read-only). This means that some matrices (e.g.
slices of tensors) can be trivially constructed from ``std::bind`` and
``std::mem_fn``. An example of this in O\ :sub:`2`\ scl_eos is how
``o2scl::eos_sn_base::slice`` generates a matrix from a 3-D tensor.

A matrix type is distinct from a "vector of vectors" or a "list of
vectors", such as that implied by
``std::vector<std::vector<double> >``. In some cases, There
are places where a vector of vectors is preferable to a matrix,
and O\ :sub:`2`\ scl expects that elements in a vector of vectors can be
accessed by ``operator[][]``. A :ref:`table <table>` object can
be thought of as a vector of vectors in this sense. The function
:cpp:func:`o2scl::tensor_grid::set_grid` also accepts a vector of
vectors, and for this function, none of the vectors needs to have
the same size. A vector of vectors can also be used to specify a
scattered list of points in a multi-dimensional space. Thus, a
vector of vectors is what is used for the argument of 
:ref:`interpm_idw <interpm_idw>`.

The word "tensor" is used for a generic object which has rank ``n``
and then has ``n`` associated indices. A vector is just a \tensor of
rank 1 and a matrix is just a \tensor of rank 2. Tensors are
implemented in O\ :sub:`2`\ scl by :ref:`tensor <tensor>`. A
multivariate function specified on a grid can be implemented in O\
:sub:`2`\ scl with :ref:`tensor_grid <tensor_grid>`. See more
discussion in the tensor section below.

Rows and columns vs. x and y
----------------------------

The most common convention is that the first index
of a matrix is the row index, i.e. to print a matrix
to the screen one uses something like::

  for(size_t row=0;row<n_rows;row++) {
    for(size_t col=0;col<n_cols;col++) {
      cout << M(row,col) << " ";
    }
    cout << endl;
  }

This is the form used in :cpp:func:`o2scl::matrix_out()` and
:cpp:func:`o2scl::array_2d_out()`. To reverse the rows and columns use
:cpp:func:`o2scl::matrix_trans_out()` and
:cpp:func:`o2scl::array_2d_trans_out()`.

A related issue is how matrices are stored. In C, two-dimensional
arrays are stored in row-major order, and the distance from the first
element to the element at ``(row,col)`` is given by
``row*n_cols+col``. In row-major order storage, the matrix elements
are stored in the same order in which they are output by the functions
:cpp:func:`o2scl::matrix_out()` and :cpp:func:`o2scl::array_2d_out()`.
The alternative is column-major order where the distance from the
first element to the element at ``(row,col)`` is given by
``col*n_rows+row``. The :ref:`tensor <tensor>` class uses a simple
generalization of row-major order. O\ :sub:`2`\ scl classes and
functions which use ``operator(,)`` operate independently of how the
data is represented in memory.

Sometimes its useful to think about the rows and columns in a
matrix as referring to a elements of a grid, and the matrix
indices refer to points in a grid in :math:`(x,y)`. It might seem
intuitive to think of a matrix as ``A[ix][iy]`` where ``ix``
and ``iy`` are the :math:`x` and :math:`y` indices because the
ordering of the indices is alphabetical. However, it is useful to
note that because functions like :cpp:func:`o2scl::matrix_out()` print
the first "row" first rather than the first column, a matrix
constructed as ``A[ix][iy]`` will be printed out with x on
the "vertical axis" and y on the "horizontal axis", which is
backwards from the usual convention when plotting data.

O\ :sub:`2`\ scl classes which interpret matrix data on a grid
(:ref:`table3d <table3d>`, :ref:`contour <contour>`, :ref:`interp2_seq
<interp2_seq>` and :ref:`interp2_direct <interp2_direct>`) use ``x``
to denote the row index and ``y`` to denote the column index by
convention.

Generic vector and matrix functions
-----------------------------------
    
GSL convenience wrappers: :ref:`gsl_vector_wrap <gsl_vector_wrap>` and 
:ref:`gsl_matrix_wrap <gsl_matrix_wrap>`.
    
Vector equality testing:
- :cpp:func:`bool vectors_equal(size_t, const vec_t &, const vec2_t &)`
- :cpp:func:`bool vectors_equal(const vec_t &, const vec2_t &)`

There are a couple functions which operate on generic vectors of any
type in ``src/base/vector.h``. This header contains functions for
sorting, summing, searching, swapping, reversing, monotonicity
testing, rotating, copying, constructing ranges, and computations of
minima and maxima. This header also contains similar operations for
matrices. For more statistically-oriented operations, see also
``src/other/vec_stats.h``. For generic functions which compute
derivatives and integrals of data specified in vectors, see
``src/deriv/vector_derint.h``. There are a few generic vector functions
related to interpolation in ``src/base/interp.h``.
    
Vector and matrix output
------------------------

For writing generic vectors to a stream, you can use
:cpp:func:`vector_out()` which is defined in ``src/base/vector.h``.
Pretty matrix output is performed by global template functions
:cpp:func:`matrix_out()` which is defined in ``src/base/columnify.h``
since it internally uses a :ref:`columnify <columnify>` object to
format the output.

Tensors
-------

Some preliminary support is provided for tensors of arbitrary rank and
size in the class :ref:`tensor <tensor>`. Classes :ref:`tensor1
<tensor1>`, :ref:`tensor2 <tensor2>`, :ref:`tensor3 <tensor3>`, and
:ref:`tensor4 <tensor4>` are rank-specific versions for 1-, 2-, 3- and
4-rank tensors. For n-dimsional data defined on a grid,
:ref:`tensor_grid <tensor_grid>` provides a space to define a
hyper-cubic grid in addition to the the tensor data. This class
:ref:`tensor_grid <tensor_grid>` also provides simple n-dimensional
interpolation of the data defined on the specified grid. There are
functions in ``src/hdf/hdf_io.h`` which provide HDF5 I/O for tensor
objects.

I/O and contiguous storage
--------------------------

O\ :sub:`2`\ scl uses HDF5 for file I/O, and in order to perform I/O
of vector-like data, HDF5 works with bare pointers. In order to
efficiently read and write vectors and other objects to HDF5 files, it
is thus important to ensure that these objects are stored contiguously
in memory. The standard template library objects, e.g. ``std::vector``
have this property as part of the recent C++ standard. The ublas
objects, so far as I know, do not necessarily have this property. For
this reason, ``o2scl::hdf_file::getd_vec`` and
``o2scl::hdf_file::setd_vec`` are efficient when working with
``std::vector`` objects, but otherwise require an extra copy upon
reading from and writing to an HDF5 file. The same holds for matrix
and tensor I/O. It is the efficiency of this I/O which motivated the
default choice of ``std::vector`` objects as the default vector type
in :ref:`table <table>` and :ref:`tensor <tensor>`. Also because of
this issue, O\ :sub:`2`\ scl does not currently provide HDF I/O
functions for :ref:`tensor <tensor>` classes unless it is built upon
``std::vector``.

Test
----

.. doxygenfile:: src/base/vector.h
