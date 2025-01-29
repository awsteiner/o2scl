Linear Algebra
==============

:ref:`O2scl <o2scl>`

Linear algebra contents
-----------------------

- :ref:`Linear algebra introduction`
- :ref:`BLAS functions`
- :ref:`Linear algebra enums`

Linear algebra introduction
---------------------------
  
O₂scl contains a small set of linear algebra routines but is also is
designed to be used with `Armadillo <http://arma.sourceforge.net>`_
and/or `Eigen <https://eigen.tuxfamily.org>`_, both of which are high
performance C++ linear algebra libraries. In the case that O₂scl was
compiled with support for either the Armadillo or Eigen libraries,
some O₂scl template functions are overloaded with the respective
Armadillo or Eigen versions.

The fallback O₂scl linear algebra routines offer a more generic and
flexible interface: they work for almost all vector and matrix types.
For matrix types which use ``operator(,)``, the linear algebra
routines routines are inside the ``o2scl_linalg`` namespaces. For
matrix types using ``operator[][]``, the linear algebra routines are
inside the ``o2scl_linalg_bracket`` namespace.

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
  * :cpp:func:`o2scl_linalg::QR_decomp_unpack()` (there is also an
    Eigen specialization with the same class name)

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

- Singular value decompositions (experimental), ``src/linalg/svd.h``

  * :cpp:func:`o2scl_linalg::balance_columns()`
  * :cpp:func:`o2scl_linalg::SV_decomp()`
  * :cpp:func:`o2scl_linalg::SV_decomp_mod()`
  * :cpp:func:`o2scl_linalg::SV_solve()`
  * :cpp:func:`o2scl_linalg::SV_decomp_jacobi()`

- Singular value decompositions step (experimental), ``src/linalg/svdstep.h``

  * :cpp:func:`o2scl_linalg::chop_small_elements()`
  * :cpp:func:`o2scl_linalg::trailing_eigenvalue()`
  * :cpp:func:`o2scl_linalg::create_schur()`
  * :cpp:func:`o2scl_linalg::svd2()`
  * :cpp:func:`o2scl_linalg::svd2_sub()`
  * :cpp:func:`o2scl_linalg::chase_out_intermediate_zero()`
  * :cpp:func:`o2scl_linalg::chase_out_trailing_zero()`
  * :cpp:func:`o2scl_linalg::chase_out_trailing_zero_sub()`
  * :cpp:func:`o2scl_linalg::qrstep()`
  * :cpp:func:`o2scl_linalg::qrstep_sub()`

- Lanczos diagonalization, ``src/linalg/lanczos.h``

  * :cpp:class:`o2scl_linalg::lanczos` (largest eigenvalues of a
    symmetric matrix)

- Linear solvers (``o2scl_linalg`` only), ``src/linalg/linear_solver.h``

  * :cpp:class:`o2scl_linalg::linear_solver_LU` (by LU decomposition)
  * :cpp:class:`o2scl_linalg::linear_solver_QR` (by QR decomposition)
  * :cpp:class:`o2scl_linalg::linear_solver_HH` (using Householder
    transformations)

- Armadillo linear solver wrapper (``o2scl_linalg`` only),
  ``src/linalg/linear_solver.h``
  
  * :cpp:class:`o2scl_linalg::linear_solver_arma`

- Eigen linear solver wrappers (``o2scl_linalg`` only),
  ``src/linalg/linear_solver.h``
    
  * :cpp:class:`o2scl_linalg::linear_solver_eigen_houseQR`
  * :cpp:class:`o2scl_linalg::linear_solver_eigen_colQR`
  * :cpp:class:`o2scl_linalg::linear_solver_eigen_fullQR`
  * :cpp:class:`o2scl_linalg::linear_solver_eigen_partLU`
  * :cpp:class:`o2scl_linalg::linear_solver_eigen_fullLU`
  * :cpp:class:`o2scl_linalg::linear_solver_eigen_LLT`
  * :cpp:class:`o2scl_linalg::linear_solver_eigen_LDLT`

- Matrix inversion interface (``o2scl_linalg`` only),
  ``src/linalg/invert.h``

  * :cpp:class:`o2scl_linalg::matrix_invert_det_LU` (inversion using
    LU decomposition)
  * :cpp:class:`o2scl_linalg::matrix_invert_det_cholesky` (inversion
    of positive symmetric matrices using the Cholesky decomposition)
  * :cpp:class:`o2scl_linalg::matrix_invert_det_arma` (Armadillo
    wrapper)
  * :cpp:class:`o2scl_linalg::matrix_invert_det_sympd_arma`
    (Armadillo wrapper for positive symmetric matrices)
  * :cpp:class:`o2scl_linalg::matrix_invert_det_eigen` (Eigen wrapper)
  * :cpp:class:`o2scl_linalg::matrix_invert_det_eigen_decomp`
    (Eigen wrapper allowing specification of the decomposition)

BLAS functions
--------------

The fallback O₂scl BLAS routines work for almost all vector and matrix
types. They also support any floating point type, so long as the
vector and matrix types are built on the same floating point type.
Specializations for double-precision numbers are preceded by the
letter ``d``.

Similar to the linear algebra functions described above, BLAS
functions for matrix types using ``operator(,)`` are inside the
``o2scl_cblas``. For matrix types using ``operator[][]``, the BLAS
functions are in the ``o2scl_cblas_bracket`` namespace.

- Level 1
  
  * :cpp:func:`o2scl_cblas::asum()`
  * :cpp:func:`o2scl_cblas::dasum()`
  * :cpp:func:`o2scl_cblas::axpy()`
  * :cpp:func:`o2scl_cblas::daxpy()`
  * :cpp:func:`o2scl_cblas::dot()`
  * :cpp:func:`o2scl_cblas::ddot()`
  * :cpp:func:`o2scl_cblas::nrm2()`
  * :cpp:func:`o2scl_cblas::dnrm2()`
  * :cpp:func:`o2scl_cblas::scal()`
  * :cpp:func:`o2scl_cblas::dscal()`

- Level 2

  * :cpp:func:`o2scl_cblas::gemv()`
  * :cpp:func:`o2scl_cblas::dgemv()`
  * :cpp:func:`o2scl_cblas::trmv()`
  * :cpp:func:`o2scl_cblas::dtrmv()`
  * :cpp:func:`o2scl_cblas::trsv()`
  * :cpp:func:`o2scl_cblas::dtrsv()`

- Level 3

  * :cpp:func:`o2scl_cblas::gemm()`
  * :cpp:func:`o2scl_cblas::dgemm()`
  * :cpp:func:`o2scl_cblas::trsm()`
  * :cpp:func:`o2scl_cblas::dtrsm()`

- Level 1 helper functions -- subvectors

  * :cpp:func:`o2scl_cblas::axpy_subvec()`
  * :cpp:func:`o2scl_cblas::daxpy_subvec()`
  * :cpp:func:`o2scl_cblas::dot_subvec()`
  * :cpp:func:`o2scl_cblas::ddot_subvec()`
  * :cpp:func:`o2scl_cblas::nrm2_subvec()`
  * :cpp:func:`o2scl_cblas::dnrm2_subvec()`
  * :cpp:func:`o2scl_cblas::scal_subvec()`
  * :cpp:func:`o2scl_cblas::dscal_subvec()`

- Level 1 helper functions -- matrix subcolumns

  * :cpp:func:`o2scl_cblas::asum_subcol()`
  * :cpp:func:`o2scl_cblas::dasum_subcol()`
  * :cpp:func:`o2scl_cblas::axpy_subcol()`
  * :cpp:func:`o2scl_cblas::daxpy_subcol()`
  * :cpp:func:`o2scl_cblas::dot_subcol()`
  * :cpp:func:`o2scl_cblas::ddot_subcol()`
  * :cpp:func:`o2scl_cblas::nrm2_subcol()`
  * :cpp:func:`o2scl_cblas::dnrm2_subcol()`
  * :cpp:func:`o2scl_cblas::scal_subcol()`
  * :cpp:func:`o2scl_cblas::dscal_subcol()`

- Level 1 helper functions -- matrix subrows

  * :cpp:func:`o2scl_cblas::axpy_subrow()`
  * :cpp:func:`o2scl_cblas::daxpy_subrow()`
  * :cpp:func:`o2scl_cblas::dot_subrow()`
  * :cpp:func:`o2scl_cblas::ddot_subrow()`
  * :cpp:func:`o2scl_cblas::nrm2_subrow()`
  * :cpp:func:`o2scl_cblas::dnrm2_subrow()`
  * :cpp:func:`o2scl_cblas::scal_subrow()`
  * :cpp:func:`o2scl_cblas::dscal_subrow()`
  
Linear algebra enums
--------------------

.. doxygenenum:: o2cblas_order

.. doxygenenum:: o2cblas_transpose

.. doxygenenum:: o2cblas_uplo

.. doxygenenum:: o2cblas_diag

.. doxygenenum:: o2cblas_side		 		 

