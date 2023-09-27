Arrays, Vectors, Matrices, and Tensors
======================================

:ref:`O2scl <o2scl>`

Arrays, vectors, matrices, and tensors contents
-----------------------------------------------

- :ref:`Vector and matrix introduction`
- :ref:`Rows and columns vs. x and y`
- :ref:`Assignment and copying`
- :ref:`Differentiation and integration`
- :ref:`Interpolating vectors`
- :ref:`Searching`
- :ref:`Vector properties`
- :ref:`Maximum and minimum functions`
- :ref:`Vector rearrangement functions`
- :ref:`Statistics functions`
- :ref:`Statistics with weighted vectors`
- :ref:`Matrix assignment and copying`
- :ref:`Matrix properties`
- :ref:`Matrix maximum and minimum functions`
- :ref:`Matrix searching`
- :ref:`Matrix rearrangement functions`
- :ref:`Tensors`
- :ref:`I/O and contiguous storage`

Vector and matrix introduction
------------------------------
     
Many useful vector and matrix objects are defined elsewhere, thus
O₂scl does not include native vector and matrix classes. Internally,
O₂scl uses ``std::vector``, Boost uBLAS vector and matrix objects:
``boost::numeric::ublas::vector<>``,
``boost::numeric::ublas::matrix<>``, and other related class
templates. Many O₂scl routines are templates which are compatible with
a wide range of vector and matrix types. See the
:ref:`Multi-dimensional solver example` which shows how an O₂scl class
can be used with Boost, Eigen, or Armadillo objects.

The O₂scl library uses a standard nomenclature to distinguish a couple
different concepts. The word "array" is typically used to refer to
C-style arrays, i.e. ``double[]``. If there are two dimensions in the
array, it is a "two-dimensional array", i.e. ``double[][]`` . The word
"vector" is reserved generic objects with array-like semantics.

In general, there are many vector types (STL, Boost, etc.) and they
can be characterized by whether or not they satisfy certain "concepts"
like ``DefaultConstructible``. O₂scl classes which operate on vector
types are designed to be as flexible as possible, so that they can be
used with almost any vector type. Eventually, all O₂scl classes with
template vector and matrix types should specify exactly which concepts
are required to be satisified, but this is still in progress.

The word "matrix" is reserved for the a generic object which has
matrix-like semantics and can be accessed using ``operator(,)``. C++
matrix types typically prefer ``operator(,)`` over ``operator[][]``.
This is because ``operator[][]`` implies the creation of a temporary
row object, and it is thus difficult to implement ``operator[]``
without incurring an overhead. Nevertheless, some O₂scl classes have
separate versions which operate on matrix types which are only
accessible with ``operator[][]`` (like two-dimensional arrays). See
:ref:`Linear Algebra` for examples of this distinction.

With ``std::function<>`` and the new lambda function support in C++11,
it is important to notice that ``std::function<double
&(size_t,size_t)>`` is also a matrix type (the ampersand is important
unless the matrix is read-only). This means that some matrices (e.g.
slices of tensors) can be trivially constructed from ``std::bind`` and
``std::mem_fn``. An example of this in O₂scl_eos is how the
:cpp:class:`o2scl::eos_sn_base::slice` class generates a matrix from a
3-D tensor.

A matrix type is distinct from a "vector of vectors" or a "list of
vectors", such as that implied by ``std::vector<std::vector<double>
>`` because in the latter, not all of the vectors in the list need to
have the same size. In some cases, There are places where a list of
vectors is preferable to a matrix, and O₂scl expects that
elements in a list of vectors can be accessed by ``operator[][]``.
The function :cpp:func:`o2scl::tensor_grid::set_grid` 
accepts a list of vectors, and for this function, none of the vectors
needs to have the same size. 

The word "tensor" is used for a generic object which has rank ``n``
and then has ``n`` associated indices. A vector is just a \tensor of
rank 1 and a matrix is just a \tensor of rank 2. Tensors are
implemented in O₂scl with class :ref:`tensor <tensor>`. A multivariate
function specified on a grid can be implemented in O₂scl with
:ref:`tensor_grid <tensor_grid>`. See more discussion in the tensor
section below.

Rows and columns vs. x and y
----------------------------

The most common convention is that the first index of a matrix is the
row index, i.e. to print a matrix to the screen one uses something
like::

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
generalization of row-major order. O₂scl classes and
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

O₂scl classes which interpret matrix data on a grid
(:ref:`table3d <table3d>`, :ref:`contour <contour>`, :ref:`interp2_seq
<interp2_seq>` and :ref:`interp2_direct <interp2_direct>`) use ``x``
to denote the row index and ``y`` to denote the column index by
convention.

Assignment and copying
----------------------

O₂scl has several functions which perform various operations
on generic vector and matrix types. These functions are listed
in the next few sections. 

- :cpp:func:`o2scl::vector_copy()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_copy_jackknife()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_set_all()` [``src/base/vector.h``]

Differentiation and integration
-------------------------------

- :cpp:func:`o2scl::vector_deriv2_interp()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_deriv2_xy_interp()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_deriv_fivept()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_deriv_fivept_tap()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_deriv_interp()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_deriv_threept()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_deriv_threept_tap()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_deriv_xy_interp()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_integ_durand()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_integ_extended4()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_integ_extended8()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_integ_interp()` [``src/base/interp``]
- :cpp:func:`o2scl::vector_integ_threept()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_integ_trap()` [``src/base/vector_derint.h``]
- :cpp:func:`o2scl::vector_integ_ul_interp()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_integ_ul_xy_interp()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_integ_xy_interp()` [``src/base/interp.h``]

Interpolating vectors
---------------------

- :cpp:func:`o2scl::vector_level_count()` [``src/base/interp.h``] 
- :cpp:func:`o2scl::vector_invert_enclosed_sum()` [``src/base/interp.h``] 
- :cpp:func:`o2scl::vector_find_level()` [``src/base/interp.h``] 
- :cpp:func:`o2scl::vector_bound_int()` [``src/base/interp.h``] 
- :cpp:func:`o2scl::vector_bound_fracint()` [``src/base/interp.h``] 
- :cpp:func:`o2scl::vector_refine()` [``src/base/interp.h``] 
- :cpp:func:`o2scl::vector_region_int()` [``src/base/interp.h``] 
- :cpp:func:`o2scl::vector_region_fracint()` [``src/base/interp.h``] 

Searching
---------

- :cpp:func:`o2scl::vector_lookup()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_search()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_bsearch()` [``src/base/search_vec.h``]
- :cpp:func:`o2scl::vector_bsearch_dec()` [``src/base/search_vec.h``]
- :cpp:func:`o2scl::vector_bsearch_inc()` [``src/base/search_vec.h``]

Vector properties
-----------------

- :cpp:func:`o2scl::vector_diffs()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vectors_equal()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vectors_equal_tol()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_is_finite()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_is_monotonic()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_is_strictly_montonic()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_largest()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_norm()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_norm_double()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_smallest()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_smallest_index()` [``src/base/vector.h``]

Maximum and minimum functions
-----------------------------

- :cpp:func:`o2scl::vector_max()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_max_index()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_max_quad()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_max_quad_loc()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_max_value()` [``src/base/vector.h``]

- :cpp:func:`o2scl::vector_min()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_min_index()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_min_quad()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_min_quad_loc()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_min_value()` [``src/base/vector.h``]

- :cpp:func:`o2scl::vector_minmax()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_minmax_index()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_minmax_value()` [``src/base/vector.h``]

Vector rearrangement functions
------------------------------

- :cpp:func:`o2scl::vector_grid()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_linear_or_log_chi2()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_linear_or_log()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_range()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_range_copy()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_rebin_xy()` [``src/base/interp.h``]
- :cpp:func:`o2scl::vector_reverse()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_reverse_double()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_rotate()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_sort()` [``src/base/vector.h``] This
  function is typically only useful for types which cannot be
  sorted with ``std::sort()``.
- :cpp:func:`o2scl::vector_sort_double()` [``src/base/vector.h``] This
  function is typically only useful for types which cannot be
  sorted with ``std::sort()``.
- :cpp:func:`o2scl::vector_spec()` [``src/hdf/hdf_io.h``]
- :cpp:func:`o2scl::vector_sum()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_sum_double()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_swap()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_swap_double()` [``src/base/vector.h``]
- :cpp:func:`o2scl::vector_to_bins()` [``src/base/vector.h``]

Statistics functions
--------------------

- :cpp:func:`o2scl::vector_absdev()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_acor()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_autocorr_tau()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_autocorr_tau_vector()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_autocorr_tau_vector_mult()`
  [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_bin_size_freedman()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_bin_size_scott()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_correlation()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_covariance()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_kurtosis()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_lag1_autocorr()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_lagk_autocorr()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_lagk_autocorr_mult()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_pvariance()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_quantile_sorted()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_roll_avg()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_skew()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_sorted_quantile()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_stddev()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_stddev_fmean()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_variance()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::vector_variance_fmean()` [``src/other/vec_stats.h``]

Statistics with weighted vectors
--------------------------------

- :cpp:func:`o2scl::wvector_absdev()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_covariance()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_factor()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_kurtosis()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_mean()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_skew()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_stddev()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_stddev_fmean()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_variance()` [``src/other/vec_stats.h``]
- :cpp:func:`o2scl::wvector_variance_fmean()` [``src/other/vec_stats.h``]

Matrix assignment and copying
-----------------------------

- :cpp:func:`o2scl::matrix_copy()`
- :cpp:func:`o2scl::matrix_set_all()`
- :cpp:func:`o2scl::matrix_set_identity()`

Matrix properties
-----------------

- :cpp:func:`o2scl::matrix_is_finite()`
- :cpp:func:`o2scl::matrix_is_lower()`
- :cpp:func:`o2scl::matrix_is_upper()`
- :cpp:func:`o2scl::matrix_make_lower()`
- :cpp:func:`o2scl::matrix_make_upper()`
- :cpp:func:`o2scl::matrix_sum()`
  
Matrix maximum and minimum functions
------------------------------------

- :cpp:func:`o2scl::matrix_max()`
- :cpp:func:`o2scl::matrix_max_index()`
- :cpp:func:`o2scl::matrix_max_value()`
- :cpp:func:`o2scl::matrix_max_value_double()`
- :cpp:func:`o2scl::matrix_min()`
- :cpp:func:`o2scl::matrix_min_index()`
- :cpp:func:`o2scl::matrix_min_value()`
- :cpp:func:`o2scl::matrix_min_value_double()`
- :cpp:func:`o2scl::matrix_minmax()`
- :cpp:func:`o2scl::matrix_minmax_index()`

Matrix searching
----------------

- :cpp:func:`o2scl::matrix_lookup()`
  
Matrix rearrangement functions
------------------------------

- :cpp:func:`o2scl::matrix_column()` [``src/base/vector.h``]
- :cpp:func:`o2scl::matrix_row()`
- :cpp:func:`o2scl::matrix_swap()`
- :cpp:func:`o2scl::matrix_swap_cols()`
- :cpp:func:`o2scl::matrix_swap_cols_double()`
- :cpp:func:`o2scl::matrix_swap_double()`
- :cpp:func:`o2scl::matrix_swap_rows()`
- :cpp:func:`o2scl::matrix_swap_rows_double()`
- :cpp:func:`o2scl::matrix_transpose()`

- :cpp:func:`o2scl::matrix_out()`
- :cpp:func:`o2scl::matrix_trans_out()`

Vector and matrix output
------------------------

For writing generic vectors to a stream, you can use the
:cpp:func:`o2scl::vector_out()` functions, which are defined in
``src/base/vector.h``. Pretty matrix output is performed by the
:cpp:func:`o2scl::matrix_out()` functions, which are defined in
``src/base/columnify.h``. The matrix output function uses a
:ref:`columnify <columnify>` object to format the output.

Tensors
-------

Tensors of arbitrary rank and size can be stored in the class
:ref:`tensor <tensor>`. Classes :ref:`tensor1 <tensor1>`,
:ref:`tensor2 <tensor2>`, :ref:`tensor3 <tensor3>`, and :ref:`tensor4
<tensor4>` are rank-specific versions for 1-, 2-, 3- and 4-rank
tensors. For n-dimsional data defined on a grid, :ref:`tensor_grid
<tensor_grid>` provides a space to define a hyper-cubic grid in
addition to the the tensor data. This class :ref:`tensor_grid
<tensor_grid>` also provides simple n-dimensional interpolation of the
data defined on the specified grid. See :ref:`File I/O with HDF5` for
functions in which provide HDF5 I/O for tensor objects.

Tensor example
--------------

.. literalinclude:: ../../../examples/ex_tensor.cpp
   :language: c++		    
   :start-after: sphinx-example-start

I/O and contiguous storage
--------------------------

O₂scl uses HDF5 for file I/O, and in order to perform I/O of
vector-like data, HDF5 works with bare pointers. In order to
efficiently read and write vectors and other objects to HDF5 files, it
is thus important to ensure that these objects are stored contiguously
in memory. The standard template library objects, e.g. ``std::vector``
have this property as part of the recent C++ standard. The ublas
objects, so far as I know, do not necessarily have this property. For
this reason, :cpp:func:`o2scl_hdf::hdf_file::getd_vec()` and
:cpp:func:`o2scl_hdf::hdf_file::setd_vec()` are efficient when working
with ``std::vector`` objects. For other vector types, one must use
:cpp:func:`o2scl_hdf::hdf_file::getd_vec_copy()` or
:cpp:func:`o2scl_hdf::hdf_file::setd_vec_copy()` which require an
extra copy upon reading from and writing to an HDF5 file. The same
holds for matrix and tensor I/O. It is the efficiency of this I/O
which motivated the default choice of ``std::vector`` objects as the
default vector type in :ref:`table <table>` and :ref:`tensor
<tensor>`. Also because of this issue, O₂scl does not currently
provide HDF I/O functions for :ref:`tensor <tensor>` classes unless
they are built upon ``std::vector``.

