Linear Algebra
==============

:ref:`O2scl <o2scl>`

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

- Householder transformations (\ref householder_base.h)
- Householder solver (\ref hh_base.h)
- LU decomposition and solver (\ref lu.h and \ref lu_base.h)
- Cholesky decomposition (\ref cholesky_base.h)
- QR decomposition and solver (\ref qr_base.h and \ref qrpt_base.h)
- Solve tridiagonal systems (\ref tridiag_base.h)
- Givens rotations (\ref givens_base.h)
- Bidiagonalizations (\ref bidiag_base.h)
- Singular value decompositions (\ref svd_base.h and \ref svdstep_base.h)
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
\ref o2scl_linalg::linear_solver_arma, 
\ref o2scl_linalg::linear_solver_eigen_houseQR, 
\ref o2scl_linalg::linear_solver_eigen_colQR, 
\ref o2scl_linalg::linear_solver_eigen_fullQR, 
\ref o2scl_linalg::linear_solver_eigen_partLU, 
\ref o2scl_linalg::linear_solver_eigen_fullLU, 
\ref o2scl_linalg::linear_solver_eigen_LLT, and
\ref o2scl_linalg::linear_solver_eigen_LDLT.

Specializations for \ref o2scl_linalg::QR_decomp_unpack() are
given in the documentation for \ref qr.h . 

Linear algebra enums
--------------------

.. doxygenenum:: o2cblas_order

.. doxygenenum:: o2cblas_transpose

.. doxygenenum:: o2cblas_uplo

.. doxygenenum:: o2cblas_diag

.. doxygenenum:: o2cblas_side		 		 

