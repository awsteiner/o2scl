Linear Algebra
==============

:ref:`O2scl <o2scl>`

Contents
--------

- :ref:`Introduction`
- :ref:`Specializations for Armadillo and Eigen`
- :ref:`Linear algebra enums`

Introduction
------------
  
There is a small set of linear algebra routines. These are not
intended to be a replacement for higher performance linear algebra
libraries. In the case that O\ :sub:`2`\ scl was compiled with support
for either the Armadillo or Eigen libraries, some O\ :sub:`2`\ scl
template functions are overloaded with the respective Armadillo or
Eigen versions.

The fallback O\ :sub:`2`\ scl linear algebra routines offer a more
generic and flexible interface: they work for almost all vector and
matrix types. For matrix types using ``operator(,)``, the BLAS and
linear algebra routines routines are inside the ``o2scl_cblas`` and
``o2scl_linalg`` namespaces. For matrix types using ``operator[][]``,
the BLAS and linear algebra routines are inside the
``o2scl_cblas_bracket`` and ``o2scl_linalg_bracket`` namespaces.

The linear algebra classes and functions include:

- Householder transformations, ``src/linalg/householder.h``

  * :cpp:func:`o2scl_linalg::householder_transform()`
  * :cpp:func:`o2scl_linalg::householder_transform_subcol()`
  * :cpp:func:`o2scl_linalg::householder_transform_subrow()`
  * :cpp:func:`o2scl_linalg::householder_hm()`
  * :cpp:func:`o2scl_linalg::householder_hm_subcol()`
  * :cpp:func:`o2scl_linalg::householder_hm_subrow()`
  * :cpp:func:`o2scl_linalg::householder_hv()`
  * :cpp:func:`o2scl_linalg::householder_hv_subcol()`
  * :cpp:func:`o2scl_linalg::householder_hm1()`
  * :cpp:func:`o2scl_linalg::householder_hm1_sub()`
  * :cpp:func:`o2scl_linalg::householder_mh()`
  * :cpp:func:`o2scl_linalg::householder_mh_subrow()`
  
- Householder solver, ``src/linalg/hh.h``

  * :cpp:func:`o2scl_linalg::HH_svx()`
  * :cpp:func:`o2scl_linalg::HH_solve()`
  
- LU decomposition and solver, ``src/linalg/lu.h``

  * :cpp:func:`o2scl_linalg::diagonal_has_zero()`
  * :cpp:func:`o2scl_linalg::LU_decomp()`
  * :cpp:func:`o2scl_linalg::LU_svx()`
  * :cpp:func:`o2scl_linalg::LU_decomp_alt()`
  * :cpp:func:`o2scl_linalg::LU_solve()`
  * :cpp:func:`o2scl_linalg::LU_refine()`
  * :cpp:func:`o2scl_linalg::LU_invert()`
  * :cpp:func:`o2scl_linalg::LU_det()`
  * :cpp:func:`o2scl_linalg::LU_lndet()`
  * :cpp:func:`o2scl_linalg::LU_sgndet()`

- Cholesky decomposition, ``src/linalg/cholesky.h``

  * :cpp:func:`o2scl_linalg::cholesky_decomp()`
  * :cpp:func:`o2scl_linalg::cholesky_det()`
  * :cpp:func:`o2scl_linalg::cholesky_solve()`
  * :cpp:func:`o2scl_linalg::cholesky_invert()`
  * :cpp:func:`o2scl_linalg::cholesky_decomp_unit()`
  
- QR decomposition, ``src/linalg/qr.h``

  * :cpp:func:`o2scl_linalg::QR_decomp()`
  * :cpp:func:`o2scl_linalg::QR_QTvec()`
  * :cpp:func:`o2scl_linalg::QR_unpack()`
  * :cpp:func:`o2scl_linalg::QR_svx()`
  * :cpp:func:`o2scl_linalg::QR_solve()`
  * :cpp:func:`o2scl_linalg::QR_update()`
  * :cpp:func:`o2scl_linalg::QR_decomp_unpack()`

- QR solver, ``src/linalg/qrpt.h``

  * :cpp:func:`o2scl_linalg::QRPT_decomp()`
  
- Solve tridiagonal systems, ``src/linalg/tridiag.h``

  * :ref:`ubvector_2_mem <ubvector_2_mem>`
  * :ref:`ubvector_4_mem <ubvector_4_mem>`
  * :ref:`ubvector_5_mem <ubvector_5_mem>`
  * :cpp:func:`o2scl_linalg::solve_tridiag_sym()`
  * :cpp:func:`o2scl_linalg::solve_tridiag_nonsym()`
  * :cpp:func:`o2scl_linalg::solve_cyc_tridiag_sym()`
  * :cpp:func:`o2scl_linalg::solve_cyc_tridiag_nonsym()`

- Givens rotations, ``src/linalg/givens.h``

  * :cpp:func:`o2scl_linalg::apply_givens_qr()`
  * :cpp:func:`o2scl_linalg::apply_givens_lq()`
  * :cpp:func:`o2scl_linalg::apply_givens_vec()`

- Bidiagonalizations, ``src/linalg/bidiag.h``

  * :cpp:func:`o2scl_linalg::bidiag_decomp()`
  * :cpp:func:`o2scl_linalg::bidiag_unpack()`
  * :cpp:func:`o2scl_linalg::bidiag_unpack2()`
  * :cpp:func:`o2scl_linalg::bidiag_unpack_B()`

- Singular value decompositions, ``src/linalg/svd.h``
  
- Singular value decompositions step, ``src/linalg/svdstep.h``
  
- Lanczos diagonalization in :ref:`lanczos <lanczos>`
  which also can compute the eigenvalues of a tridiagonal matrix.

There is also a set of linear solvers for generic matrix and
vector types which descend from :ref:`linear_solver <linear_solver>` .
These classes provide GSL-like solvers, but are generalized so
that they are compatible with vector and matrix types which allow
access through ``operator[]``.
    
Specializations for Armadillo and Eigen
---------------------------------------

Armadillo and Eigen linear solvers are wrapped to have a consistent
interface with the fallback O\ :sub:`2`\ scl linear solvers. See
:cpp:class:`o2scl_linalg::linear_solver_arma`,
:cpp:class:`o2scl_linalg::linear_solver_eigen_houseQR`,
:cpp:class:`o2scl_linalg::linear_solver_eigen_colQR`,
:cpp:class:`o2scl_linalg::linear_solver_eigen_fullQR`,
:cpp:class:`o2scl_linalg::linear_solver_eigen_partLU`,
:cpp:class:`o2scl_linalg::linear_solver_eigen_fullLU`,
:cpp:class:`o2scl_linalg::linear_solver_eigen_LLT`, and
:cpp:class:`o2scl_linalg::linear_solver_eigen_LDLT`.

Specializations for :cpp:func:`o2scl_linalg::QR_decomp_unpack()` are
 ... (see qr.h)

Linear algebra enums
--------------------

.. doxygenenum:: o2cblas_order

.. doxygenenum:: o2cblas_transpose

.. doxygenenum:: o2cblas_uplo

.. doxygenenum:: o2cblas_diag

.. doxygenenum:: o2cblas_side		 		 

