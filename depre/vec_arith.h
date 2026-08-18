/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#ifndef O2SCL_VEC_ARITH_H
#define O2SCL_VEC_ARITH_H

/** \file vec_arith.h
    \brief Vector and matrix arithmetic

    By default, the following operators are defined:
    \code
    ovector operator+(ovector_const_view &x, ovector_const_view &y);
    ovector operator+(ovector_const_view &x, uvector_const_view &y);
    ovector operator+(uvector_const_view &x, ovector_const_view &y);
    uvector operator+(uvector_const_view &x, uvector_const_view &y);

    ovector_cx operator+(ovector_cx_view &x, ovector_cx_view &y);
    ovector_cx operator+(ovector_cx_view &x, uvector_cx_view &y);
    ovector_cx operator+(uvector_cx_view &x, ovector_cx_view &y);
    uvector_cx operator+(uvector_cx_view &x, uvector_cx_view &y);
    
    ovector operator-(ovector_const_view &x, ovector_const_view &y);
    ovector operator-(ovector_const_view &x, uvector_const_view &y);
    ovector operator-(uvector_const_view &x, ovector_const_view &y);
    uvector operator-(uvector_const_view &x, uvector_const_view &y);

    ovector_cx operator-(ovector_cx_view &x, ovector_cx_view &y);
    ovector_cx operator-(ovector_cx_view &x, uvector_cx_view &y);
    ovector_cx operator-(uvector_cx_view &x, ovector_cx_view &y);
    uvector_cx operator-(uvector_cx_view &x, uvector_cx_view &y);
    
    ovector operator*(omatrix_const_view &x, ovector_const_view &y);
    ovector operator*(omatrix_const_view &x, uvector_const_view &y);
    ovector operator*(umatrix_const_view &x, ovector_const_view &y);
    uvector operator*(umatrix_const_view &x, uvector_const_view &y);

    ovector_cx operator*(omatrix_cx_view &x, ovector_cx_view &y);
    ovector_cx operator*(omatrix_cx_view &x, uvector_cx_view &y);
    ovector_cx operator*(umatrix_cx_view &x, ovector_cx_view &y);
    uvector_cx operator*(umatrix_cx_view &x, uvector_cx_view &y);

    ovector operator*(ovector_const_view &x, omatrix_const_view &y);
    ovector operator*(ovector_const_view &x, umatrix_const_view &y);
    ovector operator*(uvector_const_view &x, omatrix_const_view &y);
    uvector operator*(ovector_const_view &x, umatrix_const_view &y);

    ovector trans_mult(ovector_const_view &x, omatrix_const_view &y);
    ovector trans_mult(ovector_const_view &x, umatrix_const_view &y);
    ovector trans_mult(uvector_const_view &x, omatrix_const_view &y);
    uvector trans_mult(ovector_const_view &x, umatrix_const_view &y);

    double dot(ovector_const_view &x, ovector_const_view &y);
    double dot(ovector_const_view &x, uvector_const_view &y);
    double dot(uvector_const_view &x, ovector_const_view &y);
    double dot(uvector_const_view &x, uvector_const_view &y);

    double dot(ovector_cx_view &x, ovector_cx_view &y);
    double dot(ovector_cx_view &x, uvector_cx_view &y);
    double dot(uvector_cx_view &x, ovector_cx_view &y);
    double dot(uvector_cx_view &x, uvector_cx_view &y);

    ovector operator*(double x, ovector_const_view &y);
    uvector operator*(double x, uvector_const_view &y);
    ovector operator*(ovector_const_view &x, double y);
    uvector operator*(uvector_const_view &x, double y);
    
    ovector pair_prod(ovector_const_view &x, ovector_const_view &y);
    ovector pair_prod(ovector_const_view &x, uvector_const_view &y);
    ovector pair_prod(uvector_const_view &x, ovector_const_view &y);
    uvector pair_prod(uvector_const_view &x, uvector_const_view &y);

    ovector_cx pair_prod(ovector_cx_view &x, ovector_const_view &y);
    ovector_cx pair_prod(ovector_cx_view &x, uvector_const_view &y);
    ovector_cx pair_prod(uvector_cx_view &x, ovector_const_view &y);
    uvector_cx pair_prod(uvector_cx_view &x, uvector_const_view &y);
    ovector_cx pair_prod(ovector_const_view &x, ovector_cx_view &y);
    ovector_cx pair_prod(ovector_const_view &x, uvector_cx_view &y);
    ovector_cx pair_prod(uvector_const_view &x, ovector_cx_view &y);
    uvector_cx pair_prod(uvector_const_view &x, uvector_cx_view &y);
    ovector_cx pair_prod(ovector_cx_view &x, ovector_cx_view &y);
    ovector_cx pair_prod(ovector_cx_view &x, uvector_cx_view &y);
    ovector_cx pair_prod(uvector_cx_view &x, ovector_cx_view &y);
    uvector_cx pair_prod(uvector_cx_view &x, uvector_cx_view &y);

    bool operator==(ovector_const_view &x, ovector_const_view &y);
    bool operator==(ovector_const_view &x, uvector_const_view &y);
    bool operator==(uvector_const_view &x, ovector_const_view &y);
    bool operator==(uvector_const_view &x, uvector_const_view &y);

    bool operator!=(ovector_const_view &x, ovector_const_view &y);
    bool operator!=(ovector_const_view &x, uvector_const_view &y);
    bool operator!=(uvector_const_view &x, ovector_const_view &y);
    bool operator!=(uvector_const_view &x, uvector_const_view &y);
    \endcode      
    
    \note This used to be in a separate namespace, called 
    <tt>o2scl_arith</tt>, but this causes problems with Koenig
    lookup in template classes for operator*() when defined
    for vector addition (for example).

    \future Define operators for complex vector * real matrix
    \future Define == and != for complex vectors
    \future These should be replaced by the BLAS routines where possible?
*/

#include <iostream>
#include <complex>
#include <o2scl/ovector_tlate.h>
#include <o2scl/omatrix_tlate.h>
#include <o2scl/uvector_tlate.h>
#include <o2scl/umatrix_tlate.h>
#include <o2scl/ovector_cx_tlate.h>
#include <o2scl/omatrix_cx_tlate.h>
#include <o2scl/uvector_cx_tlate.h>
#include <o2scl/umatrix_cx_tlate.h>

#ifndef DOXYGENP
namespace o2scl {
#endif
  
  /** \brief The header macro for vector-vector addition

      Given types \c vec1, \c vec2, and \c vec_3, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      vec1 operator+(const vec2 &x, const vec3 &y);
      \endcode

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_VEC_VEC_ADD.
      
      \note This used to be in a separate namespace, called 
      <tt>o2scl_arith</tt>, but this causes problems with Koenig
      lookup in template classes for operator*() when defined
      for vector addition (for example).
  */
#define O2SCL_OP_VEC_VEC_ADD(vec1,vec2,vec3) vec1 operator+	\
    (const vec2 &x, const vec3 &y);

#ifndef DOXYGENP
  
  O2SCL_OP_VEC_VEC_ADD(o2scl::ovector,o2scl::ovector_const_view,
		       o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_ADD(o2scl::ovector,o2scl::ovector_const_view,
			 o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_ADD(o2scl::ovector,o2scl::uvector_const_view,
			 o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_ADD(o2scl::uvector,o2scl::uvector_const_view,
			 o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_ADD(o2scl::ovector_cx,o2scl::ovector_cx_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_ADD(o2scl::ovector_cx,o2scl::ovector_cx_view,
			 o2scl::uvector_cx_view)
    O2SCL_OP_VEC_VEC_ADD(o2scl::ovector_cx,o2scl::uvector_cx_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_ADD(o2scl::uvector_cx,o2scl::uvector_cx_view,
			 o2scl::uvector_cx_view)

#endif
    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif
  
  /** \brief The header macro for vector-vector subtraction

      Given types \c vec1, \c vec2, and \c vec_3, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      vec1 operator-(const vec2 &x, const vec3 &y);
      \endcode

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_VEC_VEC_SUB.
  */
#define O2SCL_OP_VEC_VEC_SUB(vec1,vec2,vec3) vec1 operator-	\
    (const vec2 &x, const vec3 &y);
  
#ifndef DOXYGENP

  O2SCL_OP_VEC_VEC_SUB(o2scl::ovector,o2scl::ovector_const_view,
		       o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_SUB(o2scl::ovector,o2scl::ovector_const_view,
			 o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_SUB(o2scl::ovector,o2scl::uvector_const_view,
			 o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_SUB(o2scl::uvector,o2scl::uvector_const_view,
			 o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_SUB(o2scl::ovector_cx,o2scl::ovector_cx_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_SUB(o2scl::ovector_cx,o2scl::ovector_cx_view,
			 o2scl::uvector_cx_view)
    O2SCL_OP_VEC_VEC_SUB(o2scl::ovector_cx,o2scl::uvector_cx_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_SUB(o2scl::uvector_cx,o2scl::uvector_cx_view,
			 o2scl::uvector_cx_view)

#endif
          
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  /** \brief The header macro for matrix-vector (right) multiplication

      Given types \c vec1, \c vec2, and \c mat, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      vec1 operator*(const mat &m, const vec3 &x);
      \endcode
      
      See also the blas version of this function \ref dgemv().

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_MAT_VEC_MULT.

  */
#define O2SCL_OP_MAT_VEC_MULT(vec1,vec2,mat) vec1 operator*	\
    (const mat &m, const vec2 &x);
  
#ifndef DOXYGENP

  O2SCL_OP_MAT_VEC_MULT(o2scl::ovector,o2scl::ovector_const_view,
			o2scl::omatrix_const_view)
    O2SCL_OP_MAT_VEC_MULT(o2scl::ovector,o2scl::ovector_const_view,
			  o2scl::umatrix_const_view)
    O2SCL_OP_MAT_VEC_MULT(o2scl::ovector,o2scl::uvector_const_view,
			  o2scl::omatrix_const_view)
    O2SCL_OP_MAT_VEC_MULT(o2scl::uvector,o2scl::uvector_const_view,
			  o2scl::umatrix_const_view)

#endif
    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  /** \brief The header macro for complex matrix-vector (right) 
      multiplication

      Given types \c vec1, \c vec2, and \c mat, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      vec1 operator*(const mat &m, const vec3 &x);
      \endcode

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_CMAT_CVEC_MULT.

  */
#define O2SCL_OP_CMAT_CVEC_MULT(vec1,vec2,mat) vec1 operator*	\
    (const mat &m, const vec2 &x);
  
#ifndef DOXYGENP

  O2SCL_OP_CMAT_CVEC_MULT(o2scl::ovector_cx,o2scl::ovector_cx_view,
			  o2scl::omatrix_cx_view)
    O2SCL_OP_CMAT_CVEC_MULT(o2scl::ovector_cx,o2scl::ovector_cx_view,
			    o2scl::umatrix_cx_view)
    O2SCL_OP_CMAT_CVEC_MULT(o2scl::ovector_cx,o2scl::uvector_cx_view,
			    o2scl::omatrix_cx_view)
    O2SCL_OP_CMAT_CVEC_MULT(o2scl::uvector_cx,o2scl::uvector_cx_view,
			    o2scl::umatrix_cx_view)

#endif

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  /** \brief The header macro for vector-matrix (left) 
      multiplication

      Given types \c vec1, \c vec2, and \c mat, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      vec1 operator*(const vec3 &x, const mat &m);
      \endcode
      
      The corresponding definition is given in 
      \ref O2SCL_OPSRC_VEC_MAT_MULT.
  */
#define O2SCL_OP_VEC_MAT_MULT(vec1,vec2,mat) vec1 operator*	\
    (const vec2 &x, const mat &m);
  
#ifndef DOXYGENP

  O2SCL_OP_VEC_MAT_MULT(o2scl::ovector,o2scl::ovector_const_view,
			o2scl::omatrix_const_view)
    O2SCL_OP_VEC_MAT_MULT(o2scl::ovector,o2scl::ovector_const_view,
			  o2scl::umatrix_const_view)
    O2SCL_OP_VEC_MAT_MULT(o2scl::ovector,o2scl::uvector_const_view,
			  o2scl::omatrix_const_view)
    O2SCL_OP_VEC_MAT_MULT(o2scl::uvector,o2scl::uvector_const_view,
			  o2scl::umatrix_const_view)

#endif

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  /// The header macro for the \c trans_mult form of vector * matrix
#define O2SCL_OP_TRANS_MULT(vec1,vec2,mat) vec1 trans_mult	\
    (const vec2 &x, const mat &m);
  
#ifndef DOXYGENP

  O2SCL_OP_TRANS_MULT(o2scl::ovector,o2scl::ovector_const_view,
		      o2scl::omatrix_const_view)
    O2SCL_OP_TRANS_MULT(o2scl::ovector,o2scl::ovector_const_view,
			o2scl::umatrix_const_view)
    O2SCL_OP_TRANS_MULT(o2scl::ovector,o2scl::uvector_const_view,
			o2scl::omatrix_const_view)
    O2SCL_OP_TRANS_MULT(o2scl::uvector,o2scl::uvector_const_view,
			o2scl::umatrix_const_view)

#endif
	    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  /** \brief The header macro for vector scalar (dot) product

      Given types \c vec1, \c vec2, and \c dtype, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      dtype operator*(const vec1 &x, const vec2 &y);
      \endcode

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_DOT_PROD.
  */
#define O2SCL_OP_DOT_PROD(dtype,vec1,vec2) dtype dot	\
    (const vec1 &x, const vec2 &y);
	    
#ifndef DOXYGENP

  O2SCL_OP_DOT_PROD(double,o2scl::ovector_const_view,
		    o2scl::ovector_const_view)
    O2SCL_OP_DOT_PROD(double,o2scl::ovector_const_view,
		      o2scl::uvector_const_view)
    O2SCL_OP_DOT_PROD(double,o2scl::uvector_const_view,
		      o2scl::ovector_const_view)
    O2SCL_OP_DOT_PROD(double,o2scl::uvector_const_view,
		      o2scl::uvector_const_view)

#endif

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  /** \brief The header macro for complex vector scalar (dot) 
      product
      
      Given types \c vec1, \c vec2, and \c dtype, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      dtype operator*(const vec1 &x, const vec2 &y);
      \endcode

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_CX_DOT_PROD.
  */
#define O2SCL_OP_CX_DOT_PROD(dtype,vec1,vec2) dtype dot	\
    (const vec1 &x, const vec2 &y);

#ifndef DOXYGENP

  O2SCL_OP_CX_DOT_PROD(gsl_complex,o2scl::ovector_cx_view,
		       o2scl::ovector_cx_view)
    O2SCL_OP_CX_DOT_PROD(gsl_complex,o2scl::ovector_cx_view,
			 o2scl::uvector_cx_view)
    O2SCL_OP_CX_DOT_PROD(gsl_complex,o2scl::uvector_cx_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_CX_DOT_PROD(gsl_complex,o2scl::uvector_cx_view,
			 o2scl::uvector_cx_view)

#endif

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  /** \brief The header macro for scalar-vector multiplication
      
      Given types \c vecv, \c vec, and \c dtype, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      vec operator*(const dtype &x, const vecv &y);
      \endcode

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_SCA_VEC_MULT.
  */
#define O2SCL_OP_SCA_VEC_MULT(dtype,vecv,vec) vec operator*	\
    (const dtype &x, const vecv &y);
  
#ifndef DOXYGENP

  O2SCL_OP_SCA_VEC_MULT(double,o2scl::ovector_const_view,o2scl::ovector)
    O2SCL_OP_SCA_VEC_MULT(double,o2scl::uvector_const_view,o2scl::uvector)

#endif
    
#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif

  /** \brief The header macro for vector-scalar multiplication
      
      Given types \c vecv, \c vec, and \c dtype, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      vec operator*(const vecv &x, const dtype &y);
      \endcode

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_VEC_SCA_MULT.
  */
#define O2SCL_OP_VEC_SCA_MULT(dtype,vecv,vec) vec operator*	\
    (const vecv &x, const dtype &y);
  
#ifndef DOXYGENP

  O2SCL_OP_VEC_SCA_MULT(double,o2scl::ovector_const_view,o2scl::ovector)
    O2SCL_OP_VEC_SCA_MULT(double,o2scl::uvector_const_view,o2scl::uvector)

#endif

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif
	      
  /** \brief The header macro for pairwise vector * vector (where either 
      vector can be real or complex)
      
      Given types \c vec1, \c vec2, and \c vec3, this macro
      provides the function declaration for adding two vectors
      using the form
      \code
      vec1 pair_prod(const vec2 &x, const vec3 &y);
      \endcode

      The corresponding definition is given in 
      \ref O2SCL_OPSRC_VEC_VEC_PRO.
  */
#define O2SCL_OP_VEC_VEC_PRO(vec1,vec2,vec3) vec1 pair_prod	\
    (const vec2 &x, const vec3 &y);
  
#ifndef DOXYGENP

  O2SCL_OP_VEC_VEC_PRO(o2scl::ovector,o2scl::ovector_const_view,
		       o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector,o2scl::ovector_const_view,
			 o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector,o2scl::uvector_const_view,
			 o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::uvector,o2scl::uvector_const_view,
			 o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_cx_view,
			 o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_cx_view,
			 o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::uvector_cx_view,
			 o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::uvector_cx,o2scl::uvector_cx_view,
			 o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_const_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_const_view,
			 o2scl::uvector_cx_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::uvector_const_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::uvector_cx,o2scl::uvector_const_view,
			 o2scl::uvector_cx_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_cx_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::ovector_cx_view,
			 o2scl::uvector_cx_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::ovector_cx,o2scl::uvector_cx_view,
			 o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_PRO(o2scl::uvector_cx,o2scl::uvector_cx_view,
			 o2scl::uvector_cx_view)

#endif
    
    /** \brief The source code macro for vector-vector addition

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_VEC_VEC_ADD
    */
#define O2SCL_OPSRC_VEC_VEC_ADD(vec1,vec2,vec3) vec1 o2scl::operator+ \
      (const vec2 &x, const vec3 &y) {					\
      size_t m=x.size();						\
      if (y.size()<m) m=y.size();					\
      vec1 r(m);							\
      for(size_t i=0;i<m;i++) {						\
	r[i]=x[i]+y[i];							\
      }									\
      return r;								\
    }
    
    /** \brief The source code macro for vector-vector subtraction

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_VEC_VEC_SUB
    */
#define O2SCL_OPSRC_VEC_VEC_SUB(vec1,vec2,vec3) vec1 o2scl::operator- \
      (const vec2 &x, const vec3 &y) {					\
      size_t m=x.size();						\
      if (y.size()<m) m=y.size();					\
      vec1 r(m);							\
      for(size_t i=0;i<m;i++) {						\
	r[i]=x[i]-y[i];							\
      }									\
      return r;								\
    }

    /** \brief The source code macro for matrix * vector

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_MAT_VEC_MULT
    */
#define O2SCL_OPSRC_MAT_VEC_MULT(vec1,vec2,mat) vec1 o2scl::operator* \
      (const mat &m, const vec2 &x) {					\
      size_t nr=m.rows();						\
      size_t nc=m.cols();						\
      vec1 res(nr);							\
      for(size_t i=0;i<nr;i++) {					\
	double r=0.0;							\
	for(size_t j=0;j<nc;j++) {					\
	  r+=m[i][j]*x[j];						\
	}								\
	res[i]=r;							\
      }									\
      return res;							\
    } 
  
    /** \brief The source code macro for complex matrix * complex vector

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_CMAT_CVEC_MULT
    */
#define O2SCL_OPSRC_CMAT_CVEC_MULT(vec1,vec2,mat) vec1 o2scl::operator* \
      (const mat &m, const vec2 &x) {					\
      size_t nr=m.rows();						\
      size_t nc=m.cols();						\
      vec1 res(nr);							\
      for(size_t i=0;i<nr;i++) {					\
	double re=0.0;							\
	double im=0.0;							\
	for(size_t j=0;j<nc;j++) {					\
	  gsl_complex g=m[i][j]*x[j];					\
	  re+=g.dat[0];							\
	  im+=g.dat[1];							\
	}								\
	res[i].dat[0]=re;						\
	res[i].dat[1]=im;						\
      }									\
      return res;							\
    } 
  
    /** \brief The source code macro for the operator form of vector * matrix

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_VEC_MAT_MULT
    */
#define O2SCL_OPSRC_VEC_MAT_MULT(vec1,vec2,mat) vec1 o2scl::operator* \
      (const vec2 &x, const mat &m) {					\
      size_t nr=m.rows();						\
      size_t nc=m.cols();						\
      vec1 res(nr);							\
      for(size_t j=0;j<nc;j++) {					\
	double r=0.0;							\
	for(size_t i=0;i<nr;i++) {					\
	  r+=x[i]*m[i][j];						\
	}								\
	res[j]=r;							\
      }									\
      return res;							\
    } 
  
    /** \brief The source code macro for the \c trans_mult form of 
	vector * matrix

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_TRANS_MULT
    */
#define O2SCL_OPSRC_TRANS_MULT(vec1,vec2,mat) vec1 o2scl::trans_mult \
      (const vec2 &x, const mat &m) {					\
      size_t nr=m.rows();						\
      size_t nc=m.cols();						\
      vec1 res(nr);							\
      for(size_t j=0;j<nc;j++) {					\
	double r=0.0;							\
	for(size_t i=0;i<nr;i++) {					\
	  r+=x[i]*m[i][j];						\
	}								\
	res[j]=r;							\
      }									\
      return res;							\
    } 
  
    /** \brief The source code macro for a vector dot product

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_DOT_PROD
    */
#define O2SCL_OPSRC_DOT_PROD(dtype,vec1,vec2) dtype o2scl::dot	\
      (const vec1 &x, const vec2 &y) {					\
      size_t m=x.size();						\
      if (y.size()<m) m=y.size();					\
      dtype r=0;							\
      for(size_t i=0;i<m;i++) {						\
	r+=x[i]*y[i];							\
      }									\
      return r;								\
    }
	    
    /** \brief The source code macro for a complex vector dot product

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_CX_DOT_PROD
    */
#define O2SCL_OPSRC_CX_DOT_PROD(dtype,vec1,vec2) dtype o2scl::dot	\
      (const vec1 &x, const vec2 &y) {					\
      size_t m=x.size();						\
      if (y.size()<m) m=y.size();					\
      dtype r={{0.0,0.0}};						\
      for(size_t i=0;i<m;i++) {						\
	r+=x[i]*y[i];							\
      }									\
      return r;								\
    }

    /** \brief The source code macro for vector=scalar*vector

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_SCA_VEC_MULT
    */
#define O2SCL_OPSRC_SCA_VEC_MULT(dtype,vecv,vec) vec o2scl::operator* \
      (const dtype &x, const vecv &y) {					\
      size_t m=y.size();						\
      vec r(m);								\
      for(size_t i=0;i<m;i++) {						\
	r[i]=x*y[i];							\
      }									\
      return r;								\
    }
  
    /** \brief The source code macro for vector=vector*scalar

	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_VEC_SCA_MULT
    */
#define O2SCL_OPSRC_VEC_SCA_MULT(dtype,vecv,vec) vec o2scl::operator* \
      (const vecv &x, const dtype &y) {					\
      size_t m=x.size();						\
      vec r(m);								\
      for(size_t i=0;i<m;i++) {						\
	r[i]=x[i]*y;							\
      }									\
      return r;								\
    }
    
    /** \brief The source code macro for pairwise vector * vector (where
	either vector can be real or complex)
	
	This define macro generates the function definition. See
	the function declaration \ref O2SCL_OP_VEC_VEC_PRO
    */
#define O2SCL_OPSRC_VEC_VEC_PRO(vec1,vec2,vec3) vec1 o2scl::pair_prod \
      (const vec2 &x, const vec3 &y) {					\
      size_t m=x.size();						\
      if (y.size()<m) m=y.size();					\
      vec1 r(m);							\
      for(size_t i=0;i<m;i++) {						\
	r[i]=x[i]*y[i];							\
      }									\
      return r;								\
    }
    
    /** \brief The header macro for vector==vector 

	Given types \c vec1 and \c vec2, this macro provides the
	function declaration for vector equality comparisons using
	\code
	bool operator==(const vec1 &x, const vec2 &y);
	\endcode

	\note Two vectors with different sizes are defined to 
	be not equal, no matter what their contents.

	The corresponding definition is given in 
	\ref O2SCL_OPSRC_VEC_VEC_EQUAL.
    */
#define O2SCL_OP_VEC_VEC_EQUAL(vec1,vec2) bool operator==	\
      (const vec1 &x, const vec2 &y);
    
#ifndef DOXYGENP
    
    O2SCL_OP_VEC_VEC_EQUAL(o2scl::ovector_const_view,
			   o2scl::ovector_const_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::ovector_const_view,
			     o2scl::uvector_const_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::uvector_const_view,
			     o2scl::ovector_const_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::uvector_const_view,
			     o2scl::uvector_const_view)
      /*
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::ovector_cx_view,
			     o2scl::ovector_const_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::ovector_cx_view,
			     o2scl::uvector_const_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::uvector_cx_view,
			     o2scl::ovector_const_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::uvector_cx_view,
			     o2scl::uvector_const_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::ovector_const_view,
			     o2scl::ovector_cx_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::ovector_const_view,
			     o2scl::uvector_cx_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::uvector_const_view,
			     o2scl::ovector_cx_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::uvector_const_view,
			     o2scl::uvector_cx_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::ovector_cx_view,
			     o2scl::ovector_cx_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::ovector_cx_view,
			     o2scl::uvector_cx_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::uvector_cx_view,
			     o2scl::ovector_cx_view)
      O2SCL_OP_VEC_VEC_EQUAL(o2scl::uvector_cx_view,
			     o2scl::uvector_cx_view)
      */

#endif

#ifdef O2SCL_NEVER_DEFINED
      }
{
#endif
	      
  /** \brief The source code macro vector==vector
	
      \note Two vectors with different sizes are defined to 
      be not equal, no matter what their contents.

      This define macro generates the function definition. See
      the function declaration \ref O2SCL_OP_VEC_VEC_EQUAL
  */
#define O2SCL_OPSRC_VEC_VEC_EQUAL(vec1,vec2) bool o2scl::operator== \
    (const vec1 &x, const vec2 &y) {					\
    size_t m=x.size();							\
    size_t n=y.size();							\
    if (m!=n) return false;						\
    for(size_t i=0;i<m;i++) {						\
      if (x[i]!=y[i]) return false;					\
    }									\
    return true;							\
  }
    
  /** \brief The header macro for vector!=vector 

      Given types \c vec1 and \c vec2, this macro provides the
      function declaration for vector inequality comparisons using
      \code
      bool operator==(const vec1 &x, const vec2 &y);
      \endcode

      \note Two vectors with different sizes are defined to 
      be not equal, no matter what their contents.
      
      The corresponding definition is given in 
      \ref O2SCL_OPSRC_VEC_VEC_NEQUAL.
  */
#define O2SCL_OP_VEC_VEC_NEQUAL(vec1,vec2) bool operator!=	\
    (const vec1 &x, const vec2 &y);
    
#ifndef DOXYGENP
    
  O2SCL_OP_VEC_VEC_NEQUAL(o2scl::ovector_const_view,
			  o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::ovector_const_view,
			    o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::uvector_const_view,
			    o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::uvector_const_view,
			    o2scl::uvector_const_view)
    /*
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::ovector_cx_view,
			    o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::ovector_cx_view,
			    o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::uvector_cx_view,
			    o2scl::ovector_const_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::uvector_cx_view,
			    o2scl::uvector_const_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::ovector_const_view,
			    o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::ovector_const_view,
			    o2scl::uvector_cx_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::uvector_const_view,
			    o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::uvector_const_view,
			    o2scl::uvector_cx_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::ovector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::ovector_cx_view,
			    o2scl::uvector_cx_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::uvector_cx_view,
			    o2scl::ovector_cx_view)
    O2SCL_OP_VEC_VEC_NEQUAL(o2scl::uvector_cx_view,
			    o2scl::uvector_cx_view)
    */

#endif

#ifdef O2SCL_NEVER_DEFINED
    }
{
#endif
	      
  /** \brief The source code macro vector!=vector
	
      \note Two vectors with different sizes are defined to 
      be not equal, no matter what their contents.

      This define macro generates the function definition. See
      the function declaration \ref O2SCL_OP_VEC_VEC_NEQUAL
  */
#define O2SCL_OPSRC_VEC_VEC_NEQUAL(vec1,vec2) bool o2scl::operator!= \
    (const vec1 &x, const vec2 &y) {					\
    size_t m=x.size();							\
    size_t n=y.size();							\
    if (m!=n) return true;						\
    for(size_t i=0;i<m;i++) {						\
      if (x[i]!=y[i]) return true;					\
    }									\
    return false;							\
  }

#ifndef DOXYGENP
  // end of o2scl namespace    
}
#endif

#endif
