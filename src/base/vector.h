/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#ifndef O2SCL_VECTOR_H
#define O2SCL_VECTOR_H

/** \file vector.h
    \brief Assorted generic vector functions

    This file contains a set of template functions which can be
    applied to almost any vector or matrix type which allow element
    access through <tt>operator[]</tt> (for vectors) or
    <tt>operator(,)</tt> for matrices. Detailed requirements on the
    template parameters are given in the functions below.

    For a general discussion of vectors and matrices in \o2, see the
    \ref vecmat_section of the User's Guide.

    For statistics operations not included here, see \ref vec_stats.h
    in the directory \c src/other . Also related are the matrix output
    functions, \ref o2scl::matrix_out(), which is defined in \ref
    columnify.h because they utilize the class \ref o2scl::columnify to
    format the output.

    For functions which search for a value in an ordered (either
    increasing or decreasing) vector, see the class \ref
    o2scl::search_vec .
*/
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_matrix.h>

#include <o2scl/misc.h>
#include <o2scl/uniform_grid.h>

namespace o2scl {

  /** \brief A simple convenience wrapper for GSL vector objects

      \warning This uses typecasts on externally allocated GSL 
      pointers and is not safe or fully const-correct. 
  */
  class gsl_vector_wrap {
    const double *d;
    size_t sz;
  public:
    gsl_vector_wrap(gsl_vector *m) {
      d=(const double *)m->data;
      sz=m->size;
    }
    double operator[](size_t i) const {
      return d[i];
    }
    size_t size() {
      return sz;
    }
  };

  /** \brief A simple convenience wrapper for GSL matrix objects

      \warning This uses typecasts on externally allocated GSL 
      pointers and is not safe or fully const-correct. 
  */
  class gsl_matrix_wrap {
    const double *d;
    size_t sz1;
    size_t sz2;
    size_t tda;
  public:
    gsl_matrix_wrap(gsl_matrix *m) {
      d=(const double *)m->data;
      sz1=m->size1;
      sz2=m->size2;
      tda=m->tda;
    }
    double operator()(size_t i, size_t j) const {
      return *(d+i*tda+j);
    }
    size_t size1() const {
      return sz1;
    }
    size_t size2() const {
      return sz2;
    }
  };

  /// \name Copying vectors and matrices in src/base/vector.h
  //@{
  /** \brief Simple vector copy

      Copy \c src to \c dest, resizing \c dest if it is too small
      to hold <tt>src.size()</tt> elements.

      This function will work for any classes \c vec_t and
      \c vec2_t which have suitably defined <tt>operator[]</tt>,
      <tt>size()</tt>, and <tt>resize()</tt> methods.
  */
  template<class vec_t, class vec2_t> 
    void vector_copy(const vec_t &src, vec2_t &dest) {
    size_t N=src.size();
    if (dest.size()<N) dest.resize(N);
    size_t i, m=N%4;
    for(i=0;i<m;i++) {
      dest[i]=src[i];
    }
    for(i=m;i+3<N;i+=4) {
      dest[i]=src[i];
      dest[i+1]=src[i+1];
      dest[i+2]=src[i+2];
      dest[i+3]=src[i+3];
    }
    return;
  }
  
  /** \brief Simple vector copy of the first N elements

      Copy the first \c N elements of \c src to \c dest.
      It is assumed that the memory allocation for \c dest
      has already been performed.

      This function will work for any class <tt>vec2_t</tt> which has
      an operator[] which returns a reference to the corresponding
      element and class <tt>vec_t</tt> with an operator[] which
      returns either a reference or the value of the corresponding
      element.
  */
  template<class vec_t, class vec2_t> 
    void vector_copy(size_t N, const vec_t &src, vec2_t &dest) {
    size_t i, m=N%4;
    for(i=0;i<m;i++) {
      dest[i]=src[i];
    }
    for(i=m;i+3<N;i+=4) {
      dest[i]=src[i];
      dest[i+1]=src[i+1];
      dest[i+2]=src[i+2];
      dest[i+3]=src[i+3];
    }
    return;
  }

  /** \brief Simple matrix copy
      
      Copy \c src to \c dest, resizing \c dest if it is too small.
      
      This function will work for any classes \c mat_t and
      \c mat2_t which have suitably defined <tt>operator()</tt>,
      <tt>size()</tt>, and <tt>resize()</tt> methods.
  */
  template<class mat_t, class mat2_t> 
    void matrix_copy(mat_t &src, mat2_t &dest) {
    size_t m=src.size1();
    size_t n=src.size2();
    if (dest.size1()<m || dest.size2()<n) dest.resize(m,n);
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	dest(i,j)=src(i,j);
      }
    }
  }

  /** \brief Simple matrix copy of the first \f$ (M,N) \f$ 
      matrix elements

      Copy the first <tt>(M,N)</tt> elements of \c src to \c dest. It
      is assumed that the memory allocation for \c dest has already
      been performed.

      This function will work for any class <tt>vec2_t</tt> which has
      an operator[][] which returns a reference to the corresponding
      element and class <tt>vec_t</tt> with an operator[][] which
      returns either a reference or the value of the corresponding
      element.
  */
  template<class mat_t, class mat2_t> 
    void matrix_copy(size_t M, size_t N, mat_t &src, mat2_t &dest) {
    for(size_t i=0;i<M;i++) {
      for(size_t j=0;j<N;j++) {
	dest(i,j)=src(i,j);
      }
    }
  }
  //@}

  /// \name Tranpositions in src/base/vector.h
  //@{
  /** \brief Simple transpose
      
      Copy the transpose of \c src to \c dest, resizing \c dest if it
      is too small.
      
      This function will work for any classes \c mat_t and
      \c mat2_t which have suitably defined <tt>operator()</tt>,
      <tt>size()</tt>, and <tt>resize()</tt> methods.
  */
  template<class mat_t, class mat2_t> 
    void matrix_transpose(mat_t &src, mat2_t &dest) {
    size_t m=src.size1();
    size_t n=src.size2();
    if (dest.size1()<n || dest.size2()<m) dest.resize(n,m);
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	dest(i,j)=src(j,i);
      }
    }
  }

  /** \brief Simple transpose of the first \f$ (m,n) \f$
      matrix elements

      Copy the transpose of the first \c m rows and the first \c cols
      of the matrix \c src into the matrix \c dest
      
      This function will work for any classes \c mat_t and \c mat2_t
      which has a suitably defined <tt>operator()</tt> method.
  */
  template<class mat_t, class mat2_t> 
    void matrix_transpose(size_t m, size_t n, mat_t &src, mat2_t &dest) {
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	dest(i,j)=src(j,i);
      }
    }
  }

  /** \brief Simple in-place transpose 
      
      Transpose the matrix \c src . If the matrix is not square,
      only the upper-left square part of the matrix will be
      transposed.
      
      This function will work for any classes \c mat_t and
      \c mat2_t which have suitably defined <tt>operator()</tt>,
      <tt>size()</tt>, and <tt>resize()</tt> methods.
  */
  template<class mat_t, class data_t> 
    void matrix_transpose(mat_t &src) {
    size_t m=src.size1();
    size_t n=src.size2();
    // Choose the smaller of n and m
    if (m<n) n=m;
    data_t tmp;
    for(size_t i=0;i<n;i++) {
      for(size_t j=0;j<n;j++) {
	tmp=src(i,j);
	src(i,j)=src(j,i);
	src(j,i)=tmp;
      }
    }
  }

  /** \brief Simple in-place transpose of the first \f$ (m,n) \f$ 
      matrix elements

      Copy the transpose of the first \c m rows and the first \c cols
      of the matrix \c src into the matrix \c dest
      
      This function will work for any classes \c mat_t and \c mat2_t
      which has a suitably defined <tt>operator()</tt> method.
  */
  template<class mat_t, class data_t> 
    void matrix_transpose(size_t m, size_t n, mat_t &src) {
    data_t tmp;
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	tmp=src(i,j);
	src(i,j)=src(j,i);
	src(j,i)=tmp;
      }
    }
  }
  //@}

  /// \name Upper and lower triangular functions in src/base/vector.h
  //@{
  /** \brief Simple test that a matrix is lower triangular
   */
  template<class mat_t> bool matrix_is_lower(mat_t &src) {
    size_t m=src.size1();
    size_t n=src.size2();
    bool ret=true;
    for(size_t i=0;i<m;i++) {
      for(size_t j=i+1;j<n;j++) {
	if (src(i,j)!=0.0) return false;
      }
    }
    return ret;
  }

  /** \brief Simple test that a matrix is upper triangular
   */
  template<class mat_t> bool matrix_is_upper(mat_t &src) {
    size_t m=src.size1();
    size_t n=src.size2();
    bool ret=true;
    for(size_t j=0;j<n;j++) {
      for(size_t i=j+1;i<m;i++) {
	if (src(i,j)!=0.0) return false;
      }
    }
    return ret;
  }
  
  /** \brief Make a matrix lower triangular by setting the upper
      triangular entries to zero
  */
  template<class mat_t> void matrix_make_lower(mat_t &src) {
    size_t m=src.size1();
    size_t n=src.size2();
    for(size_t i=0;i<m;i++) {
      for(size_t j=i+1;j<n;j++) {
	src(i,j)=0.0;
      }
    }
    return;
  }
  
  /** \brief Make a matrix upper triangular by setting the lower
      triangular entries to zero
  */
  template<class mat_t> void matrix_make_upper(mat_t &src) {
    size_t m=src.size1();
    size_t n=src.size2();
    for(size_t j=0;j<n;j++) {
      for(size_t i=j+1;i<m;i++) {
	src(i,j)=0.0;
      }
    }
    return;
  }

  /** \brief Simple test that a matrix is lower triangular
      for the first \c m rows and \c n columns
  */
  template<class mat_t> bool matrix_is_lower(size_t m, size_t n, 
					     mat_t &src) {
    bool ret=true;
    for(size_t i=0;i<m;i++) {
      for(size_t j=i+1;j<n;j++) {
	if (src(i,j)!=0.0) return false;
      }
    }
    return ret;
  }

  /** \brief Simple test that a matrix is upper triangular
      for the first \c m rows and \c n columns
  */
  template<class mat_t> bool matrix_is_upper(size_t m, size_t n, 
					     mat_t &src) {
    bool ret=true;
    for(size_t j=0;j<n;j++) {
      for(size_t i=j+1;i<m;i++) {
	if (src(i,j)!=0.0) return false;
      }
    }
    return ret;
  }
  
  /** \brief Make the first \c m rows and \c n columns of a matrix
      lower triangular by setting the upper triangular entries to zero
  */
  template<class mat_t> void matrix_make_lower(size_t m, size_t n, 
					       mat_t &src) {
    for(size_t i=0;i<m;i++) {
      for(size_t j=i+1;j<n;j++) {
	src(i,j)=0.0;
      }
    }
    return;
  }
  
  /** \brief Make the first \c m rows and \c n columns of a matrix
      upper triangular by setting the lower triangular entries to zero
  */
  template<class mat_t> void matrix_make_upper(size_t m, size_t n, 
					       mat_t &src) {
    for(size_t j=0;j<n;j++) {
      for(size_t i=j+1;i<m;i++) {
	src(i,j)=0.0;
      }
    }
    return;
  }
  //@}

  /// \name Swapping parts of vectors and matrices in src/base/vector.h
  //@{
  /** \brief Swap the first \c N elements of two vectors

      This function swaps the elements of \c v1 and \c v2, one element
      at a time. 
  */
  template<class vec_t, class vec2_t, class data_t> 
    void vector_swap(size_t N, vec_t &v1, vec2_t &v2) {
    data_t temp;
    size_t i, m=N%4;
    for(i=0;i<m;i++) {
      temp=v1[i];
      v1[i]=v2[i];
      v2[i]=temp;
    }
    for(i=m;i+3<N;i+=4) {
      temp=v1[i];
      v1[i]=v2[i];
      v2[i]=temp;
      temp=v1[i+1];
      v1[i+1]=v2[i+1];
      v2[i+1]=temp;
      temp=v1[i+2];
      v1[i+2]=v2[i+2];
      v2[i+2]=temp;
      temp=v1[i+3];
      v1[i+3]=v2[i+3];
      v2[i+3]=temp;
    }
    return;
  }

  /** \brief Swap all elements in two vectors

      This function swaps the elements of \c v1 and \c v2, one element
      at a time.

      \note It is almost always better to use <tt>std::swap</tt>
      than this function, which is provided only in cases where
      one knows one is going to be forced to use a vector type
      without a properly defined <tt>std::swap</tt> method.
  */
  template<class vec_t, class vec2_t, class data_t> 
    void vector_swap(vec_t &v1, vec2_t &v2) {
    size_t N=v1.size();
    if (v2.size()<N) N=v2.size();
    data_t temp;
    size_t i, m=N%4;
    for(i=0;i<m;i++) {
      temp=v1[i];
      v1[i]=v2[i];
      v2[i]=temp;
    }
    for(i=m;i+3<N;i+=4) {
      temp=v1[i];
      v1[i]=v2[i];
      v2[i]=temp;
      temp=v1[i+1];
      v1[i+1]=v2[i+1];
      v2[i+1]=temp;
      temp=v1[i+2];
      v1[i+2]=v2[i+2];
      v2[i+2]=temp;
      temp=v1[i+3];
      v1[i+3]=v2[i+3];
      v2[i+3]=temp;
    }
    return;
  }

  /** \brief Swap of of the first N elements of two
      double-precision vectors

      This function swaps the elements of \c v1 and \c v2, one element
      at a time.
  */
  template<class vec_t, class vec2_t>
    void vector_swap_double(size_t N, vec_t &v1, vec2_t &v2) {
    return vector_swap<vec_t,vec2_t,double>(N,v1,v2);
  }

  /** \brief Swap of all the elements in two
      double-precision vectors

      This function swaps the elements of \c v1 and \c v2, one element
      at a time.

      \note It is almost always better to use <tt>std::swap</tt>
      than this function, which is provided only in cases where
      one knows one is going to be forced to use a vector type
      without a properly defined <tt>std::swap</tt> method.
  */
  template<class vec_t, class vec2_t>
    void vector_swap_double(vec_t &v1, vec2_t &v2) {
    return vector_swap<vec_t,vec2_t,double>(v1,v2);
  }

  /** \brief Swap two elements in a vector
      
      This function swaps the element \c i and element \c j of vector
      \c v1. 
  */
  template<class vec_t, class data_t> 
    void vector_swap(vec_t &v, size_t i, size_t j) {
    data_t temp=v[i];
    v[i]=v[j];
    v[j]=temp;
    return;
  }
  
  /** \brief Swap two elements in a double-precision vector
      
      This function swaps the element \c i and element \c j of vector
      \c v1. 
      
      This function is used in \ref o2scl_linalg::QRPT_decomp() .
  */
  template<class vec_t>
    void vector_swap_double(vec_t &v, size_t i, size_t j) {
    return vector_swap<vec_t,double>(v,i,j);
  }
  
  /** \brief Swap of the first \f$ (M,N) \f$ elements in two
      matrices
      
      This function swaps the elements of \c v1 and \c v2, one element
      at a time.
  */
  template<class mat_t, class mat2_t, class data_t> 
    void matrix_swap(size_t M, size_t N, mat_t &v1, mat2_t &v2) {
    data_t temp;
    for(size_t i=0;i<M;i++) {
      for(size_t j=0;j<N;j++) {
	temp=v1[i][j];
	v1[i][j]=v2[i][j];
	v2[i][j]=temp;
      }
    }
    return;
  }

  /** \brief Swap of the first \f$ (M,N) \f$ elements in two
      double-precision matrices
      
      This function swaps the elements of \c m1 and \c m2, one element
      at a time.
  */
  template<class mat_t, class mat2_t, class data_t> 
    void matrix_swap_double(size_t M, size_t N, mat_t &m1, mat2_t &m2) {
    return matrix_swap<mat_t,mat2_t,double>(M,N,m1,m2);
  }

  /** \brief Swap two elements in a matrix

      This function swaps the element <tt>(i1,j1)</tt> and 
      element <tt>(i2,j2)</tt> of matrix \c m1. 
  */
  template<class mat_t, class data_t> 
    void matrix_swap(mat_t &m, size_t i1, size_t j1, size_t i2, size_t j2) {
    data_t temp=m(i1,j1);
    m(i1,j1)=m(i2,j2);
    m(i2,j2)=temp;
    return;
  }
  
  /** \brief Swap two elements in a double-precision matrix
      
      This function swaps the element \c i and element \c j of matrix
      \c v1. 
  */
  template<class mat_t>
    void matrix_swap_double(mat_t &m, size_t i1, size_t j1, 
			    size_t i2, size_t j2) {
    return matrix_swap<mat_t,double>(m,i1,j1,i2,j2);
  }
  
  /** \brief Swap the first \c M rows of two columns in a matrix

      This function swaps the element <tt>(i1,j1)</tt> and 
      element <tt>(i2,j2)</tt> of matrix \c m1. 
  */
  template<class mat_t, class data_t> 
    void matrix_swap_cols(size_t M, mat_t &m, size_t j1, size_t j2) {
    data_t temp;
    for(size_t i=0;i<M;i++) {
      temp=m(i,j1);
      m(i,j1)=m(i,j2);
      m(i,j2)=temp;
    }
    return;
  }
  
  /** \brief Swap the first \c M rows of two columns in a 
      double-precision matrix
      
      This function swaps the element \c i and element \c j of matrix
      \c v1. 
  */
  template<class mat_t>
    void matrix_swap_cols_double(size_t M, mat_t &m, size_t j1, size_t j2) {
    return matrix_swap_cols<mat_t,double>(M,m,j1,j2);
  }
  
  /** \brief Swap the first \c N columns of two rows in a 
      matrix

      This function swaps the element <tt>(i1,j1)</tt> and 
      element <tt>(i2,j2)</tt> of matrix \c m1. 
  */
  template<class mat_t, class data_t> 
    void matrix_swap_rows(size_t N, mat_t &m, size_t i1, size_t i2) {
    data_t temp;
    for(size_t j=0;j<N;j++) {
      temp=m(i1,j);
      m(i1,j)=m(i2,j);
      m(i2,j)=temp;
    }
    return;
  }
  
  /** \brief Swap the first \c N columns of two rows in a 
      double-precision matrix
      
      This function swaps the element \c i and element \c j of matrix
      \c v1. 
  */
  template<class mat_t>
    void matrix_swap_rows_double(size_t N, mat_t &m, size_t i1, size_t i2) {
    return matrix_swap_rows<mat_t,double>(N,m,i1,i2);
  }
  //@}

  /// \name Sorting vectors in src/base/vector.h
  //@{
  /** \brief Provide a downheap() function for vector_sort()
   */
  template<class vec_t, class data_t>
    void sort_downheap(vec_t &data, size_t n, size_t k) {
    
    data_t v=data[k];
    
    while (k<=n/2) {
      size_t j=2*k;
      
      if (j<n && data[j] < data[j+1]) j++;
      if (!(v < data[j])) break;
      data[k]=data[j];
      k=j;
    }
    
    data[k]=v;

    return;
  }

  /** \brief Sort a vector (in increasing order)

      This is a generic sorting template function using a heapsort
      algorithm. It will work for any types \c data_t and \c vec_t for
      which
      - \c data_t has a non-const version of <tt>operator=</tt>
      - \c data_t has a less than operator to compare elements
      - <tt>vec_t::operator[]</tt> returns a non-const reference
      to an object of type \c data_t 
      
      In particular, it will work with the STL template class
      <tt>std::vector</tt>, and arrays and pointers of numeric,
      character, and string objects.

      For example,
      \code
      std::string list[3]={"dog","cat","fox"};
      vector_sort<std::string[3],std::string>(3,list);
      \endcode

      \note With this function template alone, the user cannot avoid
      explicitly specifying the template types for this function
      because there is no parameter of type \c data_t, and function
      templates cannot handle default template types. For this
      reason, the function template \ref o2scl::vector_sort_double() was
      also created which provides the convenience of not requiring
      the user to specify the vector template type.

      \note This sorting routine is not stable, i.e. equal elements
      have arbtrary final ordering

      \note If \c n is zero, this function will do nothing and will
      not call the error handler.

      This works similarly to the GSL function <tt>gsl_sort_vector()</tt>.
  */
  template<class vec_t, class data_t>
    void vector_sort(size_t n, vec_t &data) {
    
    size_t N;
    size_t k;
    
    if (n==0) return;
    
    N=n-1;
    k=N/2;
    k++;
    do {
      k--;
      sort_downheap<vec_t,data_t>(data,N,k);
    } while (k > 0);
    
    while (N > 0) {
      data_t tmp=data[0];
      data[0]=data[N];
      data[N]=tmp;
      N--;
      sort_downheap<vec_t,data_t>(data,N,0);
    }
    
    return;
  }

  /** \brief Provide a downheap() function for vector_sort_index()
   */
  template<class vec_t, class vec_size_t> 
    void sort_index_downheap(size_t N, const vec_t &data, vec_size_t &order,
			     size_t k) {

    const size_t pki = order[k];
    
    while (k <= N / 2) {
      size_t j = 2 * k;
      
      if (j < N && data[order[j]] < data[order[j + 1]]) {
	j++;
      }
      
      // [GSL] Avoid infinite loop if nan
      if (!(data[pki] < data[order[j]])) {
	break;
      }
      
      order[k] = order[j];
      
      k = j;
    }
    
    order[k] = pki;
    
    return;
  }

  /** \brief Create a permutation which sorts 
      the first \c n elements of a vector (in increasing order)

      This function takes a vector \c data and arranges a list of
      indices in \c order, which give a sorted version of the vector.
      The value <tt>order[i]</tt> gives the index of entry in in \c
      data which corresponds to the <tt>i</tt>th value in the sorted
      vector. The vector \c data is unchanged by this function, and
      the initial values in \c order are ignored. Before calling this
      function, \c order must already be allocated as a vector of size
      \c n.

      For example, after calling this function, a sorted version the
      vector can be output with
      \code
      size_t n=5;
      double data[5]={3.1,4.1,5.9,2.6,3.5};
      permutation order(n);
      vector_sort_index(n,data,order);
      for(size_t i=0;i<n;i++) {
      cout << data[order[i]] << endl;
      }
      \endcode

      To create a permutation which stores as its <tt>i</tt>th element,
      the index of <tt>data[i]</tt> in the sorted vector, you can
      invert the permutation created by this function.

      This is a generic sorting template function. It will work for
      any types \c vec_t and \c vec_size_t for which
      - \c vec_t has an <tt>operator[]</tt>, and
      - \c vec_size_t has an <tt>operator[]</tt> which returns 
      a \c size_t .
      One possible type for \c vec_size_t is \ref o2scl::permutation.

      This works similarly to the GSL function <tt>gsl_sort_index()</tt>.

  */
  template<class vec_t, class vec_size_t> 
    void vector_sort_index(size_t n, const vec_t &data, vec_size_t &order) {
    size_t N;
    size_t i, k;
    
    if (n == 0) return;

    // [GSL] Set permutation to identity

    for (i = 0 ; i < n ; i++) {
      order[i] = i;
    }
    
    /* [GSL] We have n_data elements, last element is at 'n_data-1',
       first at '0' Set N to the last element number.
    */
    N = n - 1;
    
    k = N / 2;
    // [GSL] Compensate the first use of 'k--'
    k++;                          
    do {
      k--;
      sort_index_downheap<vec_t,vec_size_t>(N,data,order,k);
    } while (k > 0);
    
    while (N > 0) {
      
      // [GSL] First swap the elements 
      size_t tmp = order[0];
      order[0] = order[N];
      order[N] = tmp;
      
      // [GSL] Then process the heap
      N--;
      
      sort_index_downheap<vec_t,vec_size_t>(N,data,order,0);
    }

    return;
  }

  /** \brief Create a permutation which sorts a vector 
      (in increasing order)

      This function takes a vector \c data and arranges a list of
      indices in \c order, which give a sorted version of the vector.
      The value <tt>order[i]</tt> gives the index of entry in in \c
      data which corresponds to the <tt>i</tt>th value in the sorted
      vector. The vector \c data is unchanged by this function, and
      the initial values in \c order are ignored. Before calling this
      function, \c order must already be allocated as a vector of size
      \c n.

      For example, after calling this function, a sorted version the
      vector can be output with
      \code
      size_t n=5;
      double data[5]={3.1,4.1,5.9,2.6,3.5};
      permutation order(n);
      vector_sort_index(n,data,order);
      for(size_t i=0;i<n;i++) {
      cout << data[order[i]] << endl;
      }
      \endcode

      To create a permutation which stores as its <tt>i</tt>th element,
      the index of <tt>data[i]</tt> in the sorted vector, you can
      invert the permutation created by this function.

      This is a generic sorting template function. It will work for
      any types \c vec_t and \c vec_size_t for which
      - \c vec_t has an <tt>operator[]</tt>, and
      - \c vec_size_t has an <tt>operator[]</tt> which returns 
      a \c size_t .
      One possible type for \c vec_size_t is \ref o2scl::permutation.

      This works similarly to the GSL function <tt>gsl_sort_index()</tt>.

  */
  template<class vec_t, class vec_size_t> 
    void vector_sort_index(const vec_t &data, vec_size_t &order) {
    vector_sort_index(data.size(),data,order);
    return;
  }
  
  /** \brief Sort a vector of doubles (in increasing order)

      This function is just a wrapper for
      \code
      vector_sort<vec_t,double>(n,data);
      \endcode
      See the documentation of \ref o2scl::vector_sort() for more 
      details.
  */
  template<class vec_t>
    void vector_sort_double(size_t n, vec_t &data) {
    return vector_sort<vec_t,double>(n,data);
  }
  //@}
  
  /// \name Smallest or largest subset functions in src/base/vector.h
  //@{
  /** \brief Find the k smallest entries of the first \c n elements
      of a vector

      Given a vector \c data of size \c n, this function sets the
      first \c k entries of the vector \c smallest to the k smallest
      entries from vector \c data in ascending order. The vector \c
      smallest must be allocated beforehand to hold at least \c k
      elements.

      This works similarly to the GSL function <tt>gsl_sort_smallest()</tt>.

      \note This \f$ {\cal O}(k N) \f$ algorithm is useful only when 
      \f$ k << N \f$.

      If \c k is zero, then this function does nothing and
      returns \ref o2scl::success .
  */
  template<class vec_t, class data_t>
    void vector_smallest(size_t n, vec_t &data, size_t k, vec_t &smallest) {
    if (k>n) {
      O2SCL_ERR2("Subset length greater than size in ",
		 "function vector_smallest().",exc_einval);
    }
    if (k==0 || n==0) {
      O2SCL_ERR2("Vector size zero or k zero in ",
		 "function vector_smallest().",exc_einval);
    }

    // Take the first element
    size_t j=1;
    data_t xbound=data[0];
    smallest[0]=xbound;

    // Examine the remaining elements
    for(size_t i=1;i<n;i++) {
      data_t xi=data[i];
      if (j<k) {
	j++;
      } else if (xi>=xbound) {
	continue;
      }
      size_t i1;
      for(i1=j-1;i1>0;i1--) {
	if (xi>smallest[i1-1]) break;
	smallest[i1]=smallest[i1-1];
      }
      smallest[i1]=xi;
      xbound=smallest[j-1];
    }
    return;
  }

  /** \brief Find the k smallest entries of a vector
      of a vector

      Given a vector \c data, this function sets the first \c k
      entries of the vector \c smallest to the k smallest entries from
      vector \c data in ascending order. The vector \c smallest
      is resized if necessary to hold at least \c k elements.

      This works similarly to the GSL function <tt>gsl_sort_smallest()</tt>.

      \note This \f$ {\cal O}(k N) \f$ algorithm is useful only when 
      \f$ k << N \f$.

      If \c k is zero, then this function does nothing and
      returns \ref o2scl::success .
  */
  template<class vec_t, class data_t>
    void vector_smallest(vec_t &data, size_t k, vec_t &smallest) {
    size_t n=data.size();
    if (smallest.size()<k) smallest.resize(k);
    return vector_smallest(n,data,k,smallest);
  }

  /** \brief Find the indexes of the k smallest entries among the
      first \c n entries of a vector

      Given a vector \c data, this function sets the first \c k
      entries of the vector \c smallest equal to the indexes of the k
      smallest entries from vector \c data in ascending order. The
      vector \c smallest is resized if necessary to hold at least \c k
      elements.

      \note This \f$ {\cal O}(k N) \f$ algorithm is useful only when 
      \f$ k << N \f$.

      If \c k is zero or \c n is zero or \f$ k > n\f$, then this
      function calls the error handler.
  */
  template<class vec_t, class data_t, class vec_size_t>
    void vector_smallest_index(size_t n, const vec_t &data, size_t k,
			       vec_size_t &index) {
    if (k>n) {
      O2SCL_ERR2("Subset length greater than size in ",
		 "function vector_smallest_index().",exc_einval);
    }
    if (k==0 || n==0) {
      O2SCL_ERR2("Vector size zero or k zero in ",
		 "function vector_smallest_index().",exc_einval);
    }

    index.resize(k);

    // [GSL] Take the first element
    size_t j=1;
    data_t xbound=data[0];
    index[0]=0;

    // [GSL] Examine the remaining elements
    for(size_t i=1;i<n;i++) {
      data_t xi=data[i];
      if (j<k) {
	j++;
      } else if (xi>=xbound) {
	continue;
      }
      size_t i1;
      for(i1=j-1;i1>0;i1--) {
	if (xi>data[index[i1-1]]) break;
	index[i1]=index[i1-1];
      }
      index[i1]=i;
      xbound=data[index[j-1]];
    }
    return;
  }

  /** \brief Find the indexes of the k smallest entries of a vector
   */
  template<class vec_t, class data_t, class vec_size_t>
    void vector_smallest_index(const vec_t &data, size_t k,
			       vec_size_t &index) {
    size_t n=data.size();
    if (index.size()<k) index.resize(k);
    return o2scl::vector_smallest_index<vec_t,data_t,vec_size_t>
      (n,data,k,index);
  }

  /** \brief Find the k largest entries of the first \c n elements
      of a vector

      Given a vector \c data of size \c n this sets the first \c k
      entries of the vector \c largest to the k largest entries from
      vector \c data in descending order. The vector \c largest must
      be allocated beforehand to hold at least \c k elements.

      This works similarly to the GSL function <tt>gsl_sort_largest()</tt>.

      \note This \f$ {\cal O}(k N) \f$ algorithm is useful only when 
      \f$ k << N \f$.

      If \c k is zero, then this function does nothing and
      returns \ref o2scl::success .
  */
  template<class vec_t, class data_t>
    void vector_largest(size_t n, vec_t &data, size_t k, vec_t &largest) {
    if (k>n) {
      O2SCL_ERR2("Subset length greater than size in ",
		 "function vector_largest().",exc_einval);
    }
    if (k==0 || n==0) {
      O2SCL_ERR2("Vector size zero or k zero in ",
		 "function vector_largest().",exc_einval);
    }

    // Take the first element
    size_t j=1;
    data_t xbound=data[0];
    largest[0]=xbound;

    // Examine the remaining elements
    for(size_t i=1;i<n;i++) {
      data_t xi=data[i];
      if (j<k) {
	j++;
      } else if (xi<=xbound) {
	continue;
      }
      size_t i1;
      for(i1=j-1;i1>0;i1--) {
	if (xi<largest[i1-1]) break;
	largest[i1]=largest[i1-1];
      }
      largest[i1]=xi;
      xbound=largest[j-1];
    }
    return;
  }
  
  /** \brief Find the k largest entries of a vector
      of a vector
      
      Given a vector \c data, this function sets the first \c k
      entries of the vector \c largest to the k largest entries from
      vector \c data in ascending order. The vector \c largest
      is resized if necessary to hold at least \c k elements.
      
      This works similarly to the GSL function <tt>gsl_sort_largest()</tt>.
      
      \note This \f$ {\cal O}(k N) \f$ algorithm is useful only when 
      \f$ k << N \f$.
      
      If \c k is zero, then this function does nothing and
      returns \ref o2scl::success .
  */
  template<class vec_t, class data_t>
    void vector_largest(vec_t &data, size_t k, vec_t &largest) {
    size_t n=data.size();
    if (largest.size()<k) largest.resize(k);
    return vector_largest(n,data,k,largest);
  }
  //@}
  
  /// \name Vector minimum and maximum functions in src/base/vector.h
  //@{
  /** \brief Compute the maximum of the first \c n elements of a vector
   */
  template<class vec_t, class data_t>
    data_t vector_max_value(size_t n, const vec_t &data) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_max_value().",exc_efailed);
    }
    data_t max=data[0];
    for(size_t i=1;i<n;i++) {
      if (data[i]>max) {
	max=data[i];
      }
    }
    return max;
  }

  /** \brief Compute the maximum value of a vector
   */
  template<class vec_t, class data_t>
    data_t vector_max_value(const vec_t &data) {

    size_t n=data.size();
    if (n==0) {
      O2SCL_ERR("Sent empty vector to vector_max_value().",exc_efailed);
    }
    data_t max=data[0];
    for(size_t i=1;i<n;i++) {
      if (data[i]>max) {
	max=data[i];
      }
    }
    return max;
  }

  /** \brief Compute the index which holds the 
      maximum of the first \c n elements of a vector
  */
  template<class vec_t, class data_t>
    size_t vector_max_index(size_t n, const vec_t &data) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_max_index().",exc_efailed);
    }
    data_t max=data[0];
    size_t ix=0;
    for(size_t i=1;i<n;i++) {
      if (data[i]>max) {
	max=data[i];
	ix=i;
      }
    }
    return ix;
  }

  /** \brief Compute the maximum of the first \c n elements of a vector
   */
  template<class vec_t, class data_t>
    void vector_max(size_t n, const vec_t &data, size_t &index, 
		    data_t &val) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_max().",exc_efailed);
    }
    val=data[0];
    index=0;
    for(size_t i=1;i<n;i++) {
      if (data[i]>val) {
	val=data[i];
	index=i;
      }
    }
    return;
  }

  /** \brief Compute the minimum of the first \c n elements of a vector
   */
  template<class vec_t, class data_t>
    data_t vector_min_value(size_t n, const vec_t &data) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_min_value().",exc_efailed);
    }
    data_t min=data[0];
    for(size_t i=1;i<n;i++) {
      if (data[i]<min) {
	min=data[i];
      }
    }
    return min;
  }

  /** \brief Compute the minimum value in a vector
   */
  template<class vec_t, class data_t>
    data_t vector_min_value(const vec_t &data) {

    size_t n=data.size();
    if (n==0) {
      O2SCL_ERR("Sent empty vector to vector_min_value().",exc_efailed);
    }
    data_t min=data[0];
    for(size_t i=1;i<n;i++) {
      if (data[i]<min) {
	min=data[i];
      }
    }
    return min;
  }

  /** \brief Compute the index which holds the 
      minimum of the first \c n elements of a vector
  */
  template<class vec_t, class data_t>
    size_t vector_min_index(size_t n, const vec_t &data) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_min_index().",exc_efailed);
    }
    data_t min=data[0];
    size_t ix=0;
    for(size_t i=1;i<n;i++) {
      if (data[i]<min) {
	min=data[i];
	ix=i;
      }
    }
    return ix;
  }

  /** \brief Compute the minimum of the first \c n elements of a vector
   */
  template<class vec_t, class data_t>
    void vector_min(size_t n, const vec_t &data, size_t &index, 
		    data_t &val) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_min().",exc_efailed);
    }
    val=data[0];
    index=0;
    for(size_t i=1;i<n;i++) {
      if (data[i]<val) {
	val=data[i];
	index=i;
      }
    }
    return;
  }

  /** \brief Compute the minimum and maximum of the first 
      \c n elements of a vector
  */
  template<class vec_t, class data_t>
    void vector_minmax_value(size_t n, vec_t &data, 
			     data_t &min, data_t &max) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_min().",exc_efailed);
    }
    min=data[0];
    max=min;
    for(size_t i=1;i<n;i++) {
      if (data[i]<min) {
	min=data[i];
      }
      if (data[i]>max) {
	max=data[i];
      }
    }
    return;
  }

  /** \brief Compute the minimum and maximum of the first 
      \c n elements of a vector
  */
  template<class vec_t, class data_t>
    void vector_minmax_index(size_t n, vec_t &data, 
			     size_t &ix_min, size_t &ix_max) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_min().",exc_efailed);
    }
    data_t min=data[0];
    data_t max=min;
    ix_min=0;
    ix_max=0;
    for(size_t i=1;i<n;i++) {
      if (data[i]<min) {
	min=data[i];
	ix_min=i;
      }
      if (data[i]>max) {
	max=data[i];
	ix_max=i;
      }
    }
    return;
  }

  /** \brief Compute the minimum and maximum of the first 
      \c n elements of a vector
  */
  template<class vec_t, class data_t>
    void vector_minmax(size_t n, vec_t &data, 
		       size_t &ix_min, data_t &min, 
		       size_t &ix_max, data_t &max) {
    
    if (n==0) {
      O2SCL_ERR("Sent size=0 to vector_min().",exc_efailed);
    }
    min=data[0];
    max=min;
    ix_min=0;
    ix_max=0;
    for(size_t i=1;i<n;i++) {
      if (data[i]<min) {
	min=data[i];
	ix_min=i;
      }
      if (data[i]>max) {
	max=data[i];
	ix_max=i;
      }
    }
    return;
  }
  //@}
  
  /// \name Extrema of vectors through quadratic fit in src/base/vector.h
  //@{
  /** \brief Maximum of vector by quadratic fit
   */
  template<class vec_t, class data_t>
    data_t vector_max_quad(size_t n, const vec_t &data) {
    size_t ix=vector_max_index<vec_t,data_t>(n,data);
    if (ix==0) {
      return quadratic_extremum_y<data_t>(0,1,2,data[0],data[1],data[2]);
    } else if (ix==n-1) {
      return quadratic_extremum_y<data_t>
	(n-3,n-2,n-1,data[n-3],data[n-2],data[n-1]);
    } 
    return quadratic_extremum_y<data_t>
      (ix-1,ix,ix+1,data[ix-1],data[ix],data[ix+1]);
  }

  /** \brief Maximum of vector by quadratic fit
   */
  template<class vec_t, class data_t>
    data_t vector_max_quad(size_t n, const vec_t &x, const vec_t &y) {
    size_t ix=vector_max_index<vec_t,data_t>(n,y);
    if (ix==0 || ix==n-1) return y[ix];
    return quadratic_extremum_y<data_t>(x[ix-1],x[ix],x[ix+1],
					y[ix-1],y[ix],y[ix+1]);
  }

  /** \brief Location of vector maximum by quadratic fit
   */
  template<class vec_t, class data_t>
    data_t vector_max_quad_loc(size_t n, const vec_t &x, const vec_t &y) {
    size_t ix=vector_max_index<vec_t,data_t>(n,y);
    if (ix==0 || ix==n-1) return y[ix];
    return quadratic_extremum_x<data_t>(x[ix-1],x[ix],x[ix+1],
					y[ix-1],y[ix],y[ix+1]);
  }

  /** \brief Minimum of vector by quadratic fit
   */
  template<class vec_t, class data_t>
    data_t vector_min_quad(size_t n, const vec_t &data) {
    size_t ix=vector_min_index<vec_t,data_t>(n,data);
    if (ix==0) {
      return quadratic_extremum_y<data_t>(0,1,2,data[0],data[1],data[2]);
    } else if (ix==n-1) {
      return quadratic_extremum_y<data_t>
	(n-3,n-2,n-1,data[n-3],data[n-2],data[n-1]);
    } 
    return quadratic_extremum_y<data_t>
      (ix-1,ix,ix+1,data[ix-1],data[ix],data[ix+1]);
  }

  /** \brief Minimum of vector by quadratic fit
   */
  template<class vec_t, class data_t>
    data_t vector_min_quad(size_t n, const vec_t &x, const vec_t &y) {
    size_t ix=vector_min_index<vec_t,data_t>(n,y);
    if (ix==0 || ix==n-1) return y[ix];
    return quadratic_extremum_y<data_t>(x[ix-1],x[ix],x[ix+1],
					y[ix-1],y[ix],y[ix+1]);
  }

  /** \brief Location of vector minimum by quadratic fit
   */
  template<class vec_t, class data_t>
    data_t vector_min_quad_loc(size_t n, const vec_t &x, const vec_t &y) {
    size_t ix=vector_min_index<vec_t,data_t>(n,y);
    if (ix==0 || ix==n-1) return y[ix];
    return quadratic_extremum_x<data_t>(x[ix-1],x[ix],x[ix+1],
					y[ix-1],y[ix],y[ix+1]);
  }
  //@}

  /// \name Matrix minimum and maximum functions in src/base/vector.h
  //@{
  /** \brief Compute the maximum of the lower-left part of a matrix
   */
  template<class mat_t, class data_t>
    data_t matrix_max_value(size_t m, const size_t n, const mat_t &data) {
    
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::itos(m)+","+o2scl::itos(n)+") in "+
	"matrix_max_value().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    data_t max=data(0,0);
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)>max) {
	  max=data(i,j);
	}
      }
    }
    return max;
  }

  /** \brief Compute the maximum of a matrix
   */
  template<class mat_t, class data_t> data_t
    matrix_max_value(const mat_t &data) {
    size_t m=data.size1();
    size_t n=data.size2();
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_max_value().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    data_t max=data(0,0);
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)>max) {
	  max=data(i,j);
	}
      }
    }
    return max;
  }

  /** \brief Compute the maximum of a matrix
   */
  template<class mat_t> double
    matrix_max_value_double(const mat_t &data) {
    size_t m=data.size1();
    size_t n=data.size2();
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_max_value_double().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    double max=data(0,0);
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)>max) {
	  max=data(i,j);
	}
      }
    }
    return max;
  }

  /** \brief Compute the maximum of a matrix and return 
      the indices of the maximum element
  */
  template<class mat_t, class data_t>
    void matrix_max_index(size_t m, size_t n, const mat_t &data,
			  size_t &i_max, size_t &j_max, data_t &max) {
    
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_max_index().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    max=data(0,0);
    i_max=0;
    j_max=0;
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)>max) {
	  max=data(i,j);
	  i_max=i;
	  j_max=j;
	}
      }
    }
    return;
  }

  /** \brief Compute the maximum of a matrix and return 
      the indices of the maximum element
  */
  template<class mat_t, class data_t>
    void matrix_max_index(const mat_t &data,
			  size_t &i_max, size_t &j_max, data_t &max) {

    size_t m=data.size1();
    size_t n=data.size2();
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_max_index().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    max=data(0,0);
    i_max=0;
    j_max=0;
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)>max) {
	  max=data(i,j);
	  i_max=i;
	  j_max=j;
	}
      }
    }
    return;
  }

  /** \brief Compute the minimum of a matrix
   */
  template<class mat_t, class data_t>
    data_t matrix_min_value(size_t m, size_t n, const mat_t &data) {
    
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_min_value().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    data_t min=data(0,0);
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)<min) {
	  min=data(i,j);
	}
      }
    }
    return min;
  }

  /** \brief Compute the minimum of a matrix
   */
  template<class mat_t, class data_t>
    data_t matrix_min_value(const mat_t &data) {
    
    size_t m=data.size1();
    size_t n=data.size2();
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_min_value().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    data_t min=data(0,0);
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)<min) {
	  min=data(i,j);
	}
      }
    }
    return min;
  }

  /** \brief Compute the minimum of a matrix
   */
  template<class mat_t>
    double matrix_min_value_double(const mat_t &data) {
    
    size_t m=data.size1();
    size_t n=data.size2();
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_min_value().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    double min=data(0,0);
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)<min) {
	  min=data(i,j);
	}
      }
    }
    return min;
  }

  /** \brief Compute the minimum of a matrix and return 
      the indices of the minimum element
  */
  template<class mat_t, class data_t>
    void matrix_min_index(size_t n, size_t m, const mat_t &data,
			  size_t &i_min, size_t &j_min, data_t &min) {
    
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_min_index().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    min=data(0,0);
    i_min=0;
    j_min=0;
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)<min) {
	  min=data(i,j);
	  i_min=i;
	  j_min=j;
	}
      }
    }
    return;
  }

  /** \brief Compute the minimum of a matrix and return 
      the indices of the minimum element
  */
  template<class mat_t, class data_t>
    void matrix_min_index(const mat_t &data,
			  size_t &i_min, size_t &j_min, data_t &min) {
    
    size_t m=data.size1();
    size_t n=data.size2();
    if (m==0 || n==0) {
      std::string str=((std::string)"Matrix with zero size (")+
	o2scl::szttos(m)+","+o2scl::szttos(n)+") in "+
	"matrix_min_index().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }
    min=data(0,0);
    i_min=0;
    j_min=0;
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	if (data(i,j)<min) {
	  min=data(i,j);
	  i_min=i;
	  j_min=j;
	}
      }
    }
    return;
  }

  /** \brief Compute the minimum and maximum of a matrix
   */
  template<class mat_t, class data_t>
    void matrix_minmax(size_t n, size_t m, const mat_t &data,
		       data_t &min, data_t &max) {
    
    if (n==0 || m==0) {
      O2SCL_ERR("Sent size=0 to matrix_min().",exc_efailed);
    }
    min=data(0,0);
    max=data(0,0);
    for(size_t i=0;i<n;i++) {
      for(size_t j=0;j<m;j++) {
	if (data(i,j)<min) {
	  min=data(i,j);
	} else if (data(i,j)>max) {
	  max=data(i,j);
	}
      }
    }
    return;
  }

  /** \brief Compute the minimum and maximum of a matrix
   */
  template<class mat_t, class data_t>
    void matrix_minmax(const mat_t &data,
		       data_t &min, data_t &max) {
    return matrix_minmax<mat_t,data_t>
      (data.size1(),data.size2(),data,min,max);
  }
    
  /** \brief Compute the minimum and maximum of a matrix and
      return their locations
  */
  template<class mat_t, class data_t>
    void matrix_minmax_index(size_t n, size_t m, const mat_t &data,
			     size_t &i_min, size_t &j_min, data_t &min, 
			     size_t &i_max, size_t &j_max, data_t &max) {
    
    if (n==0 || m==0) {
      O2SCL_ERR2("Sent size=0 to function ",
		 "matrix_minmax_index().",exc_efailed);
    }
    min=data(0,0);
    i_min=0;
    j_min=0;
    max=data(0,0);
    i_max=0;
    j_max=0;
    for(size_t i=0;i<n;i++) {
      for(size_t j=0;j<m;j++) {
	if (data(i,j)<min) {
	  min=data(i,j);
	  i_min=i;
	  j_min=j;
	} else if (data(i,j)>max) {
	  max=data(i,j);
	  i_max=i;
	  j_max=j;
	}
      }
    }
    return;
  }

  /** \brief Compute the sum of matrix elements
   */
  template<class mat_t, class data_t>
    data_t matrix_sum(size_t m, size_t n, const mat_t &data) {
    
    data_t sum=0.0;
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	sum+=data(i,j);
      }
    }
    return sum;
  }

  /** \brief Compute the sum of matrix elements
   */
  template<class mat_t, class data_t>
    data_t matrix_sum(const mat_t &data) {
    return matrix_sum<mat_t,data_t>(data.size1(),data.size2(),data);
  }
  //@}

  /// \name Searching vectors and matrices in src/base/vector.h
  //@{
  /** \brief Lookup the value \c x0 in the first \c n elements of 
      vector \c x

      The function finds the element among the first \c n elements of
      \c x which is closest to the value \c x0. It ignores all
      elements in \c x which are not finite. If the vector is empty,
      or if all of the first \c n elements in \c x are not finite,
      then the error handler will be called.
      
      This function works for all classes \c vec_t where an operator[]
      is defined which returns a double (either as a value or a
      reference).
  */
  template<class vec_t>
    size_t vector_lookup(size_t n, const vec_t &x, double x0) {
    if (n==0) {
      O2SCL_ERR("Empty vector in function vector_lookup().",
		exc_einval);
      return 0;
    }
    size_t row=0, i=0;
    while(!std::isfinite(x[i]) && i<n-1) i++;
    if (i==n-1) {
      O2SCL_ERR2("Entire vector not finite in ",
		 "function vector_lookup()",exc_einval);
      return 0;
    }
    double best=x[i], bdiff=fabs(x[i]-x0);
    for(;i<n;i++) {
      if (std::isfinite(x[i]) && fabs(x[i]-x0)<bdiff) {
	row=i;
	best=x[i];
	bdiff=fabs(x[i]-x0);
      }
    }
    return row;
  }

  /** \brief Lookup element \c x0 in vector \c x

      This function finds the element in vector \c x which is closest
      to \c x0. It ignores all elements in \c x which are not finite.
      If the vector is empty, or if all of the
      elements in \c x are not finite, then the error handler will be
      called.
      
      This function works for all classes \c vec_t with a
      <tt>size()</tt> method and where an operator[] is defined which
      returns a double (either as a value or a reference).
  */
  template<class vec_t>
    size_t vector_lookup(const vec_t &x, double x0) {
    return vector_lookup(x.size(),x,x0);
  }

  /** \brief Lookup an element in the first $(m,n)$ entries in a matrix

      Return the location <tt>(i,j)</tt> of the element closest to 
      \c x0. 
  */
  template<class mat_t>
    void matrix_lookup(size_t m, size_t n, const mat_t &A, 
		       double x0, size_t &i, size_t &j) {
    if (m==0 || n==0) {
      O2SCL_ERR("Empty matrix in matrix_lookup().",
		exc_einval);
    }
    double dist=0.0;
    bool found_one=false;
    for(size_t i2=0;i2<m;i2++) {
      for(size_t j2=0;j2<n;j2++) {
	if (std::isfinite(A(i,j))) {
	  if (found_one==false) {
	    dist=fabs(A(i,j)-x0);
	    found_one=true;
	    i=i2;
	    j=j2;
	  } else {
	    if (fabs(A(i,j)-x0)<dist) {
	      dist=fabs(A(i,j)-x0);
	      i=i2;
	      j=j2;
	    }
	  }
	}
      }
    }
    if (found_one==false) {
      O2SCL_ERR2("Entire matrix not finite in ",
		 "function matrix_lookup()",exc_einval);
    }
    return;
  }

  /** \brief Lookup an element in a matrix

      Return the location <tt>(i,j)</tt> of the element closest to 
      \c x0. 
  */
  template<class mat_t>
    void matrix_lookup(const mat_t &A, 
		       double x0, size_t &i, size_t &j) {
    matrix_lookup(A.size1(),A.size2(),x0,i,j);
    return;
  }

  /** \brief Binary search a part of an increasing vector for
      <tt>x0</tt>.
      
      This function performs a binary search of between
      <tt>x[lo]</tt> and <tt>x[hi]</tt>. It returns 
      - \c lo if \c x0 < <tt>x[lo]</tt>
      - \c i if <tt>x[i]</tt> <= \c x0 < <tt>x[i+2]</tt> 
      for \c lo <= \c i < \c hi
      - \c hi-1 if \c x0 >= \c <tt>x[hi-1]</tt>

      This function is designed to find the interval containing \c x0,
      not the index of the element closest to x0. To perform the
      latter operation, you can use \ref vector_lookup().

      The element at <tt>x[hi]</tt> is never referenced by this
      function. The parameter \c hi can be either the index of the
      last element (e.g. <tt>n-1</tt> for a vector of size <tt>n</tt>
      with starting index <tt>0</tt>), or the index of one element
      (e.g. <tt>n</tt> for a vector of size <tt>n</tt> and a starting
      index <tt>0</tt>) for a depending on whether or not the user
      wants to allow the function to return the index of the last
      element.

      This function operates in the same way as
      <tt>gsl_interp_bsearch()</tt>.

      The operation of this function is undefined if the data is
      not strictly monotonic, i.e. if some of the data elements are 
      equal.

      This function will call the error handler if \c lo is
      greater than \c hi.
  */
  template<class vec_t, class data_t> 
    size_t vector_bsearch_inc(const data_t x0, const vec_t &x,
			      size_t lo, size_t hi) {
    if (lo>hi) {
      O2SCL_ERR2("Low and high indexes backwards in ",
		 "function vector_bsearch_inc().",exc_einval);
    }
    while (hi>lo+1) {
      size_t i=(hi+lo)/2;
      if (x[i]>x0) {
	hi=i;
      } else {
	lo=i;
      }
    }
    
    return lo;
  }

  /** \brief Binary search a part of an decreasing vector for
      <tt>x0</tt>.
      
      This function performs a binary search of between
      <tt>x[lo]</tt> and <tt>x[hi]</tt> (inclusive). It returns 
      - \c lo if \c x0 > <tt>x[lo]</tt>
      - \c i if <tt>x[i]</tt> >= \c x0 > <tt>x[i+1]</tt> 
      for \c lo <= \c i < \c hi
      - \c hi-1 if \c x0 <= \c <tt>x[hi-1]</tt>
    
      The element at <tt>x[hi]</tt> is never referenced by this
      function. The parameter \c hi can be either the index of the
      last element (e.g. <tt>n-1</tt> for a vector of size <tt>n</tt>
      with starting index <tt>0</tt>), or the index of one element
      (e.g. <tt>n</tt> for a vector of size <tt>n</tt> and a starting
      index <tt>0</tt>) for a depending on whether or not the user
      wants to allow the function to return the index of the last
      element.

      The operation of this function is undefined if the data is
      not strictly monotonic, i.e. if some of the data elements are 
      equal.
      
      This function will call the error handler if \c lo is
      greater than \c hi.
  */
  template<class vec_t, class data_t> 
    size_t vector_bsearch_dec(const data_t x0, const vec_t &x,
			      size_t lo, size_t hi) {
    if (lo>hi) {
      O2SCL_ERR2("Low and high indexes backwards in ",
		 "function vector_bsearch_dec().",exc_einval);
    }
    while (hi>lo+1) {
      size_t i=(hi+lo)/2;
      if (x[i]<x0) {
	hi=i;
      } else {
	lo=i;
      }
    }
    
    return lo;
  }

  /** \brief Binary search a part of a monotonic vector for
      <tt>x0</tt>.

      This wrapper just calls \ref o2scl::vector_bsearch_inc() or
      \ref o2scl::vector_bsearch_dec() depending on the ordering
      of \c x. 
  */
  template<class vec_t, class data_t> 
    size_t vector_bsearch(const data_t x0, const vec_t &x,
			  size_t lo, size_t hi) {
    if (x[lo]<x[hi-1]) {
      return vector_bsearch_inc<vec_t,data_t>(x0,x,lo,hi);
    }
    return vector_bsearch_dec<vec_t,data_t>(x0,x,lo,hi);
  }
  
  /** \brief Binary search a monotonic vector for
      <tt>x0</tt>.
      
      This function calls \ref o2scl::vector_bsearch_inc() or
      \ref o2scl::vector_bsearch_dec() depending on the ordering
      of \c x. 
  */
  template<class vec_t, class data_t> 
    size_t vector_bsearch(const data_t x0, const vec_t &x) {
    size_t lo=0;
    size_t hi=x.size();
    if (x[lo]<x[hi-1]) {
      return vector_bsearch_inc<vec_t,data_t>(x0,x,lo,hi);
    }
    return vector_bsearch_dec<vec_t,data_t>(x0,x,lo,hi);
  }
  //@}

  /// \name Ordering and finite tests in src/base/vector.h
  //@{
  /** \brief Test if the first \c n elements of a vector are 
      monotonic and increasing or decreasing

      If \c n is zero or one, this function will return 0 without
      calling the error handler. If all the vector's elements are equal,
      this function will return 3. Otherwise, if the vector is not
      monotonic, then this function will return 0. Finally, if the
      vector is nondecreasing (increasing or equal intervals), this
      function will return 1, and if the vector is nonincreasing
      (decreasing or equal intervals), this function will return 2.
      This function assumes that simple comparison operators have been
      defined for the type of each vector element.
  */
  template<class vec_t>
    int vector_is_monotonic(size_t n, vec_t &data) {

    if (n<1) return 0;
    if (n<2) {
      if (data[0]==data[1]) {
	return 3;
      } else if (data[0]<data[1]) {
	return 1;
      } else {
	return 2;
      }
    }
    
    // Find first non-flat interval
    size_t start=0;
    bool done=false;
    for(size_t i=0;i<n-1 && done==false;i++) {
      if (data[i]!=data[i+1]) {
	done=true;
      } else {
	start++;
      }
    }

    // If all elements in the vector are equal, or if the only
    // distinct element is at the end, then return true.
    if (done==false) {
      return 3;
    } 
    
    if (start==n-2) {
      if (data[start]<data[start+1]) return 1;
      else return 2;
    }

    // Determine if the vector is increasing (inc=true) or decreasing
    // (inc=false)
    bool inc=true;
    if (data[start]>data[start+1]) inc=false;

    if (inc) {
      for(size_t i=start+1;i<n-1;i++) {
	// If there is one decreasing interval, return false
	if (data[i]>data[i+1]) return 0;
      }
      return 1;
    }
    
    // If there is one increasing interval, return false
    for(size_t i=start+1;i<n-1;i++) {
      if (data[i]<data[i+1]) return 0;
    }
    return 2;
  }

  /** \brief Test if the first \c n elements of a vector are 
      monotonic and increasing or decreasing

      If \c n is zero or one, this function will return 0 without
      calling the error handler. If all the vector's elements are equal,
      this function will return 3. Otherwise, if the vector is not
      monotonic, then this function will return 0. Finally, if the
      vector is nondecreasing (increasing or equal intervals), this
      function will return 1, and if the vector is nonincreasing
      (decreasing or equal intervals), this function will return 2.
      This function assumes that simple comparison operators have been
      defined for the type of each vector element.
  */
  template<class vec_t> int vector_is_monotonic(vec_t &data) {
    return vector_is_monotonic(data.size(),data);
  }

  /** \brief Test if the first \c n elements of a vector are 
      strictly monotonic and determine if they are increasing or decreasing

      If \c n is zero this function will return 0 without calling the
      error handler. Also, if the vector is not monotonic, this
      function will return 0. If the vector is strictly 
      monotonic, then this function will return 1 if it is 
      increasing and 2 if it is decreasing.
  */
  template<class vec_t>
    int vector_is_strictly_monotonic(size_t n, vec_t &data) {
    
    if (n<1) return 0;
    
    // Determine if the vector is increasing (inc=true) or decreasing
    // (inc=false)
    bool inc=true;
    if (data[0]==data[1]) {
      return 0;
    } else if (data[0]>data[1]) {
      inc=false;
    }

    if (inc) {
      for(size_t i=1;i<n-1;i++) {
	// If there is one nonincreasing interval, return 0
	if (data[i]>=data[i+1]) return 0;
      }
      return 1;
    } 

    // If there is one increasing interval, return 0
    for(size_t i=1;i<n-1;i++) {
      if (data[i]<=data[i+1]) return 0;
    }
    return 2;
  }

  /** \brief Test if the first \c n elements of a vector are 
      strictly monotonic and determine if they are increasing or decreasing

      If \c n is zero this function will return 0 without calling the
      error handler. Also, if the vector is not monotonic, this
      function will return 0. If the vector is strictly 
      monotonic, then this function will return 1 if it is 
      increasing and 2 if it is decreasing.
  */
  template<class vec_t>
    int vector_is_strictly_monotonic(vec_t &data) {
    return vector_is_strictly_monotonic(data.size(),data);
  }

  /** \brief Test if the first \c n elements of a vector are finite

      If \c n is zero, this will return true without throwing
      an exception.

      The corresponding tests for matrix functions are
      in clbas_base.h .
  */
  template<class vec_t>
    bool vector_is_finite(size_t n, vec_t &data) {
    for(size_t i=0;i<n;i++) {
      if (!std::isfinite(data[i])) return false;
    }
    return true;
  }

  /** \brief Test if a vector is finite

      If \c n is zero, this will return true without throwing
      an exception.

      The corresponding tests for matrix functions are
      in clbas_base.h .
  */
  template<class vec_t> bool vector_is_finite(vec_t &data) {
    return vector_is_finite(data.size(),data);
  }
  //@}

  /// \name Miscellaneous mathematical functions in src/base/vector.h
  //@{
  /** \brief Compute the sum of the first \c n elements of a vector

      If \c n is zero, this will return 0 without throwing
      an exception.
  */
  template<class vec_t, class data_t>
    data_t vector_sum(size_t n, vec_t &data) {
    
    data_t sum=0.0;
    for(size_t i=0;i<n;i++) {
      sum+=data[i];
    }
    return sum;
  }

  /** \brief Create a new vector containing the differences between
      adjacent entries
  */
  template<class vec_t, class rvec_t>
    void vector_diffs(const vec_t &v_data, rvec_t &v_diffs) {

    size_t n=v_data.size();
    v_diffs.resize(n-1);
    for(size_t i=0;i<n-1;i++) {
      v_diffs[i]=v_data[i+1]-v_data[i];
    }
    
    return;
  }
  
  /** \brief Compute the sum of all the elements of a vector

      If the vector has zero size, this will return 0 without
      calling the error handler.
  */
  template<class vec_t, class data_t> data_t vector_sum(vec_t &data) {
    data_t sum=0.0;
    for(size_t i=0;i<data.size();i++) {
      sum+=data[i];
    }
    return sum;
  }

  /** \brief Compute the sum of the first \c n elements of a vector
      of double-precision numbers

      If \c n is zero, this will return 0 without throwing
      an exception.
  */
  template<class vec_t>double vector_sum_double(size_t n, vec_t &data) {
    double sum=0.0;
    for(size_t i=0;i<n;i++) {
      sum+=data[i];
    }
    return sum;
  }

  /** \brief Compute the sum of all the elements of a vector
      of double-precision numbers

      If the vector has zero size, this will return 0 without
      calling the error handler.
  */
  template<class vec_t> double vector_sum_double(vec_t &data) {
    double sum=0.0;
    for(size_t i=0;i<data.size();i++) {
      sum+=data[i];
    }
    return sum;
  }

  /** \brief Compute the norm of the first \c n entries of a 
      vector of floating-point (single or double precision) numbers

      This function is a more generic version of 
      \ref o2scl_cblas::dnrm2 . 
  */
  template<class vec_t, class data_t>
    data_t vector_norm(size_t n, const vec_t &x) {
    
    data_t scale = 0.0;
    data_t ssq = 1.0;
    
    if (n <= 0) {
      return 0.0;
    } else if (n == 1) {
      return fabs(x[0]);
    }
      
    for (size_t i = 0; i < n; i++) {
      const data_t xx = x[i];

      if (xx != 0.0) {
	const data_t ax = fabs(xx);
          
	if (scale < ax) {
	  ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
	  scale = ax;
	} else {
	  ssq += (ax / scale) * (ax / scale);
	}
      }

    }
      
    return scale * sqrt(ssq);
  }

  /** \brief Compute the norm of a vector of floating-point 
      (single or double precision) numbers
  */
  template<class vec_t, class data_t> data_t vector_norm(const vec_t &x) {
    return vector_norm<vec_t,data_t>(x.size(),x);
  }

  /** \brief Compute the norm of the first \c n entries of a 
      vector of double precision numbers

      This function is a more generic version of 
      \ref o2scl_cblas::dnrm2 . 
  */
  template<class vec_t>
    double vector_norm_double(size_t n, const vec_t &x) {
    
    double scale = 0.0;
    double ssq = 1.0;
    
    if (n <= 0) {
      return 0.0;
    } else if (n == 1) {
      return fabs(x[0]);
    }
      
    for (size_t i = 0; i < n; i++) {
      const double xx = x[i];

      if (xx != 0.0) {
	const double ax = fabs(xx);
          
	if (scale < ax) {
	  ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
	  scale = ax;
	} else {
	  ssq += (ax / scale) * (ax / scale);
	}
      }

    }
      
    return scale * sqrt(ssq);
  }

  /** \brief Compute the norm of a vector of double precision numbers
   */
  template<class vec_t> double vector_norm_double(const vec_t &x) {
    return vector_norm_double<vec_t>(x.size(),x);
  }
  //@}

  /// \name Other vector and matrix functions in src/base/vector.h
  //@{
  /** \brief Set the first N entries in a vector to a particular value
   */
  template<class vec_t, class data_t> 
    void vector_set_all(size_t N, vec_t &src, data_t val) {
    for(size_t i=0;i<N;i++) {
      src[i]=val;
    }
    return;
  }
  
  /** \brief Set all entries in a vector to a particular value
   */
  template<class vec_t, class data_t> 
    void vector_set_all(vec_t &src, data_t val) {
    o2scl::vector_set_all(src.size(),src,val);
    return;
  }
  
  /** \brief Set the first (M,N) entries in a matrix to a particular value
   */
  template<class mat_t, class data_t> 
    void matrix_set_all(size_t M, size_t N, mat_t &src, data_t val) {
    for(size_t i=0;i<M;i++) {
      for(size_t j=0;j<N;j++) {
	src(i,j)=val;
      }
    }
    return;
  }
  
  /** \brief Set all entries in a matrix to a particular value
   */
  template<class mat_t, class data_t> 
    void matrix_set_all(mat_t &src, data_t val) {
    o2scl::matrix_set_all(src.size1(),src.size2(),src,val);
    return;
  }
  
  /** \brief From a given vector, create a new vector by removing a
      specified element

      This funciton is used in \ref o2scl::interp_krige_optim::qual_fun() .
  */
  template<class vec_t, class vec2_t> 
    void vector_copy_jackknife(const vec_t &src, size_t iout, vec2_t &dest) {
    if (src.size()==0) {
      O2SCL_ERR("Empty source vector.",o2scl::exc_einval);
    }
    if (iout>=src.size()) {
      O2SCL_ERR("Requested element beyond end.",o2scl::exc_einval);
    }
    dest.resize(src.size()-1);
    size_t j=0;
    for(size_t i=0;i<src.size();i++) {
      if (i!=iout) {
	dest[j]=src[i];
	j++;
      }
    }
    return;
  }

  /** \brief From a given vector, create a new vector by removing a
      specified element

      This funciton is used in \ref o2scl::interp_krige_optim::qual_fun() .
  */
  template<class vec_t, class vec2_t> 
    void vector_copy_jackknife(size_t sz, const vec_t &src,
			       size_t iout, vec2_t &dest) {
			       
    if (sz==0) {
      O2SCL_ERR("Empty source vector.",o2scl::exc_einval);
    }
    if (iout>=sz) {
      O2SCL_ERR("Requested element beyond end.",o2scl::exc_einval);
    }
    dest.resize(sz-1);
    size_t j=0;
    for(size_t i=0;i<sz;i++) {
      if (i!=iout) {
	dest[j]=src[i];
	j++;
      }
    }
    return;
  }

  /** \brief "Rotate" a vector so that the kth element is now the beginning

      This is a generic template function which will work for
      any types \c data_t and \c vec_t for which
      - \c data_t has an <tt>operator=</tt>
      - <tt>vec_t::operator[]</tt> returns a reference
      to an object of type \c data_t 

      This function is used, for example, in \ref o2scl::pinside.

      \note This function is not the same as a Givens rotation,
      which is typically referred to in BLAS routines as <tt>drot()</tt>.
  */
  template<class vec_t, class data_t>
    void vector_rotate(size_t n, vec_t &data, size_t k) {

    data_t *tmp=new data_t[n];
    for(size_t i=0;i<n;i++) {
      tmp[i]=data[(i+k)%n];
    }
    for(size_t i=0;i<n;i++) {
      data[i]=tmp[i];
    }
    delete[] tmp;
    
    return;
  }
  
  /** \brief Reverse the first \c n elements of a vector
      
      If \c n is zero, this function will silently do nothing.
  */
  template<class vec_t, class data_t>
    void vector_reverse(size_t n, vec_t &data) {
    data_t tmp;

    for(size_t i=0;i<n/2;i++) {
      tmp=data[n-1-i];
      data[n-1-i]=data[i];
      data[i]=tmp;
    }
    return;
  }

  /** \brief Reverse a vector

      If the <tt>size()</tt> method returns zero, this function will
      silently do nothing.
  */
  template<class vec_t, class data_t>
    void vector_reverse(vec_t &data) {
    data_t tmp;
    size_t n=data.size();

    for(size_t i=0;i<n/2;i++) {
      tmp=data[n-1-i];
      data[n-1-i]=data[i];
      data[i]=tmp;
    }
    return;
  }

  /** \brief Reverse the first n elements in a vector of double
      precision numbers

      If \c n is zero, this function will silently do nothing.
  */
  template<class vec_t>
    void vector_reverse_double(size_t n, vec_t &data) {
    double tmp;

    for(size_t i=0;i<n/2;i++) {
      tmp=data[n-1-i];
      data[n-1-i]=data[i];
      data[i]=tmp;
    }
    return;
  }

  /** \brief Reverse a vector of double precision numbers

      If the <tt>size()</tt> method returns zero, this function will
      silently do nothing.
  */
  template<class vec_t> void vector_reverse_double(vec_t &data) {
    double tmp;
    size_t n=data.size();
    
    for(size_t i=0;i<n/2;i++) {
      tmp=data[n-1-i];
      data[n-1-i]=data[i];
      data[i]=tmp;
    }
    return;
  }
  
  /** \brief Trivial index vector
      
      This object just returns the index whenever an object in the
      vector is requested, i.e. <tt>operator[](i)</tt> always returns
      \c i.
   */
  template<class data_t> class vector_index_vector {
  public:
    data_t operator[](size_t &i) const {
      return i;
    }
  };
  
  /** \brief Index vector with a size method

      This object just returns the index whenever an object in the
      vector is requested, i.e. <tt>operator[](i)</tt> always returns
      \c i.
   */
  template<class data_t> class vector_index_vector_size {

  protected:
    
    /// The vector size
    size_t n;
    
  public:

    /** \brief Create an index vector with size \c n_
     */
    vector_index_vector_size(size_t n_) {
      n=n_;
    }

    /** \brief Obtain the element with index \c i
     */
    data_t operator[](size_t &i) const {
      if (i>=n) {
	O2SCL_ERR("Out of bounds.",o2scl::exc_einval);
      }
      return i;
    }

    /** \brief Get the size of the vector
     */
    size_t size() const {
      return n;
    }

    /** \brief Resize the index vector
     */
    void resize(size_t n_) {
      n=n_;
    }
    
  };
  
  /** \brief Construct a row of a matrix

      This class template works with combinations of ublas
      <tt>matrix</tt> and <tt>matrix_row</tt> objects,
      <tt>arma::mat</tt> and <tt>arma::rowvec</tt>, and
      <tt>Eigen::MatrixXd</tt> and <tt>Eigen::VectorXd</tt>.

      \note When calling this function with ublas objects, the
      namespace prefix <tt>"o2scl::"</tt> must often be specified,
      otherwise some compilers will use argument dependent lookup and
      get (justifiably) confused with matrix_row in the ublas
      namespace.

      \note The template parameters must be explicitly specified
      when calling this template function. 
  */
  template<class mat_t, class mat_row_t>
    mat_row_t matrix_row(mat_t &M, size_t row) {
    return mat_row_t(M,row);
  }

  /** \brief Generic object which represents a row of a matrix

      \note This class is experimental.

      This class is used in <tt>o2scl::eos_sn_base::slice</tt>
      to construct a row of a matrix object of type
      \code 
      std::function<double &(size_t,size_t)>
      \endcode
  */
  template<class mat_t> class matrix_row_gen {

  protected:

    /// A reference to the original matrix
    mat_t &m_;

    /// The selected row
    size_t row_;

  public:

    /// Create a row object from row \c row of matrix \c m 
  matrix_row_gen(mat_t &m, size_t row) : m_(m), row_(row) {
    }
    
    /// Return a reference to the ith column of the selected row
    double &operator[](size_t i) {
      return m_(row_,i);
    }
    
    /// Return a const reference to the ith column of the selected row
    const double &operator[](size_t i) const {
      return m_(row_,i);
    }
  };

  /** \brief Matrix row object with a constructor and resize method

      This is used in \ref o2scl::ode_iv_solve_grid .
   */
  template<class mat_t> class matrix_row_gen_ctor {

  protected:

    /// A pointer to the matrix
    mat_t *mp;

    /// The selected row
    size_t row_;

    /// A matrix to point to
    mat_t mat;

  public:

    /// Create a row object from row \c row of matrix \c m 
  matrix_row_gen_ctor(mat_t &m, size_t row) : mp(&m), row_(row) {
    }

    /// Create a row object from row \c row of matrix \c m 
    matrix_row_gen_ctor(size_t n_cols=0) {
      if (n_cols==0) {
	mp=0;
      } else {
	mat.resize(1,n_cols);
	mp=&mat;
	row_=0;
      }
    }

    /** \brief Resize
     */
    void resize(size_t n_cols=0) {
      if (n_cols==0) {
	mp=0;
	mat.resize(0,0);
      } else {
	mat.resize(1,n_cols);
	mp=&mat;
	row_=0;
      }
      return;
    }
    
    /** \brief Return size
     */
    size_t size() const {
      if (mp==0) {
	return 0;
      }
      return mp->size2();
    }
    
    /// Return a reference to the ith column of the selected row
    double &operator[](size_t i) {
      if (mp==0) {
	O2SCL_ERR("No matrix in matrix_row_gen_ctor::operator[].",
		  o2scl::exc_efailed);
      }
      return (*mp)(row_,i);
    }
    
    /// Return a const reference to the ith column of the selected row
    const double &operator[](size_t i) const {
      if (mp==0) {
	O2SCL_ERR("No matrix in matrix_row_gen_ctor::operator[].",
		  o2scl::exc_efailed);
      }
      return (*mp)(row_,i);
    }
  };

  /** \brief Construct a view of the transpose of a matrix

      \note This class is experimental.
  */
  template<class mat_t> class matrix_view_transpose {

  protected:

    /// A reference to the original matrix
    mat_t &m_;

  public:

    /// Create a row object from row \c row of matrix \c m 
  matrix_view_transpose(mat_t &m) : m_(m) {
    }
    
    /// Return a reference to the ith column of the selected row
    double &operator()(size_t i, size_t j) {
      return m_(j,i);
    }
    
    /// Return a const reference to the ith column of the selected row
    const double &operator()(size_t i, size_t j) const {
      return m_(j,i);
    }

    /** \brief Return the number of rows
     */
    size_t size1() const {
      return m_.size2();
    }
    
    /** \brief Return the number of columns
     */
    size_t size2() const {
      return m_.size1();
    }
  
    
  };

  /** \brief Construct a view of a matrix omtting a specified row

      \note This class is experimental.
  */
  template<class mat_t> class matrix_view_omit_row {

  protected:

    /// A reference to the original matrix
    mat_t &m_;

    size_t ro;

  public:

    /// Create
  matrix_view_omit_row(mat_t &m, size_t row_omit) : m_(m) {
      ro=row_omit;
    }
    
    /// Return a reference
    double &operator()(size_t i, size_t j) {
      if (i>=ro) {
	return m_(i+1,j);
      }
      return m_(i,j);
    }
    
    /// Return a const reference
    const double &operator()(size_t i, size_t j) const {
      if (i>=ro) {
	return m_(i+1,j);
      }
      return m_(i,j);
    }

    /** \brief Return the number of rows
     */
    size_t size1() const {
      return m_.size1()-1;
    }
    
    /** \brief Return the number of columns
     */
    size_t size2() const {
      return m_.size2();
    }
  
    
  };

  /** \brief Construct a view of a matrix omitting one specified column

      \note This class is experimental.
  */
  template<class mat_t> class matrix_view_omit_column {

  protected:

    /// A reference to the original matrix
    mat_t &m_;

    size_t co;

  public:

    /// Create
  matrix_view_omit_column(mat_t &m, size_t column_omit) : m_(m) {
      co=column_omit;
    }
    
    /// Return a reference
    double &operator()(size_t i, size_t j) {
      if (j>=co) {
	return m_(i,j+1);
      }
      return m_(i,j);
    }
    
    /// Return a const reference
    const double &operator()(size_t i, size_t j) const {
      if (j>=co) {
	return m_(i,j+1);
      }
      return m_(i,j);
    }

    /** \brief Return the number of rows
     */
    size_t size1() const {
      return m_.size1();
    }
    
    /** \brief Return the number of columns
     */
    size_t size2() const {
      return m_.size2()-1;
    }
  
    
  };

  /** \brief Generic object which represents a row of a const matrix

      \note This class is experimental.

      This class is used in <tt>o2scl::eos_sn_base::slice</tt>
      to construct a row of a matrix object of type
      \code 
      std::function<double &(size_t,size_t)>
      \endcode
  */
  template<class mat_t> class const_matrix_row_gen {

  protected:

    /// A reference to the original matrix
    const mat_t &m_;

    /// The selected row
    size_t row_;

  public:

    /// Create a row object from row \c row of matrix \c m 
  const_matrix_row_gen(const mat_t &m, size_t row) : m_(m), row_(row) {
    }
    
    /// Return a const reference to the ith column of the selected row
    const double &operator[](size_t i) const {
      return m_(row_,i);
    }
  };

  /** \brief A simple matrix view object
   */
  class const_matrix_view {
  
  public:
  
    /** \brief Return a reference to the element at row \c row
	and column \c col
    */
    const double &operator()(size_t row, size_t col) const;
    
    /** \brief Return the number of rows
     */
    size_t size1() const;
    
    /** \brief Return the number of columns
     */
    size_t size2() const;
  
  };

  /** \brief A simple matrix view object
   */
  class matrix_view {
  
  public:
  
    /** \brief Return a reference to the element at row \c row
	and column \c col
    */
    const double &operator()(size_t row, size_t col) const;
    
    /** \brief Return a reference to the element at row \c row
	and column \c col
    */
    double &operator()(size_t row, size_t col);
    
    /** \brief Return the number of rows
     */
    size_t size1() const;
    
    /** \brief Return the number of columns
     */
    size_t size2() const;
  
  };

  /** \brief View a o2scl::table object as a matrix

      \note This stores a pointer to the table and the user must ensure
      that the pointer is valid with the matrix view is accessed.

      \future It would be nice to store a reference rather than a
      pointer, but this causes problems with \ref o2scl::interpm_idw .
  */
  template<class vec1_t, class vec2_t=std::vector<vec1_t> > 
    class matrix_view_vec_vec : public matrix_view {
  
  protected:
  
  /// Pointer to the table
  vec2_t *vvp;

  public:

  /** \brief Swap method
   */
  friend void swap(matrix_view_vec_vec &t1,
		   matrix_view_vec_vec &t2) {
    /// Just swap the pointer
    std::swap(t1.vvp,t2.vvp);
    return;
  }

  matrix_view_vec_vec() {
    vvp=0;
  }
  
  /** \brief Create a matrix view object from the specified 
      table and list of rows
  */
  matrix_view_vec_vec(vec2_t &vv) {
    vvp=&vv;
  }
  
  /** \brief Return the number of rows
   */
  size_t size1() const {
    if (vvp==0) return 0;
    return vvp->size();
  }
  
  /** \brief Return the number of columns
   */
  size_t size2() const {
    if (vvp==0) return 0;
    if (vvp->size()==0) return 0;
    return (*vvp)[0].size();
  }
  
  /** \brief Return a reference to the element at row \c row
      and column \c col
  */
  const double &operator()(size_t row, size_t col) const {
    if (vvp==0) {
      O2SCL_ERR2("Object empty in ",
		 "matrix_view_vec_vec::operator().",
		 o2scl::exc_einval);
    }
    if (row>=vvp->size()) {
      O2SCL_ERR2("Row exceeds max in ",
		 "matrix_view_vec_vec::operator().",
		 o2scl::exc_einval);
    }
    if (col>=(*vvp)[row].size()) {
      O2SCL_ERR2("Column exceeds max in ",
		 "matrix_view_vec_vec::operator().",
		 o2scl::exc_einval);
    }
    return (*vvp)[row][col];
  }
    
  /** \brief Return a reference to the element at row \c row
      and column \c col
  */
  double &operator()(size_t row, size_t col) {
      if (vvp==0) {
      O2SCL_ERR2("Object empty in ",
                "matrix_view_vec_vec::operator().",
                o2scl::exc_einval);
    }
    if (row>=vvp->size()) {
      O2SCL_ERR2("Row exceeds max in ",
		 "matrix_view_vec_vec::operator().",
		 o2scl::exc_einval);
    }
    if (col>=(*vvp)[row].size()) {
      std::cout << row << " " << col << " "
      << (*vvp)[row].size() << std::endl;
      O2SCL_ERR2("Column exceeds max in ",
		 "matrix_view_vec_vec::operator().",
		 o2scl::exc_einval);
    }
    return (*vvp)[row][col];
  }
    
  };

  /** \brief Construct a column of a matrix

      This class template works with combinations of ublas
      <tt>matrix</tt> and <tt>matrix_column</tt> objects,
      <tt>arma::mat</tt> and <tt>arma::colvec</tt>, and
      <tt>Eigen::MatrixXd</tt> and <tt>Eigen::VectorXd</tt>.

      \note When calling this function with ublas objects, the
      namespace prefix <tt>"o2scl::"</tt> must often be specified,
      otherwise some compilers will use argument dependent lookup and
      get (justifiably) confused with matrix_column in the ublas
      namespace.

      \note The template parameters must be explicitly specified
      when calling this template function.
  */
  template<class mat_t, class mat_column_t>
    mat_column_t matrix_column(mat_t &M, size_t column) {
    return mat_column_t(M,column);
  }

  /** \brief Generic object which represents a column of a matrix

      \note This class is experimental. 

      The only requirement on the type <tt>mat_t</tt> is that
      it must have an operator(size_t,size_t) method which
      accesses elements in the matrix.

      This class is used in <tt>o2scl::eos_sn_base::slice</tt>
      to construct a row of a matrix object of type
      \code 
      std::function<double &(size_t,size_t)>
      \endcode
  */
  template<class mat_t> class matrix_column_gen {
    
  protected:

    /// A reference to the original matrix
    mat_t &m_;

    /// The selected column
    size_t column_;
    
  public:
    
    /// Create a column object from column \c column of matrix \c m 
  matrix_column_gen(mat_t &m, size_t column) : m_(m), column_(column) {
    }
    
    /// Return a reference to the ith row of the selected column
    double &operator[](size_t i) {
      return m_(i,column_);
    }
    
    /// Return a const reference to the ith row of the selected column
    const double &operator[](size_t i) const {
      return m_(i,column_);
    }
    
  };

  /** \brief Generic object which represents a column of a const matrix

      \note This class is experimental. 
      
      The only requirement on the type <tt>mat_t</tt> is that
      it must have an operator(size_t,size_t) method which
      accesses elements in the matrix.

      This class is used in one of
      the \ref o2scl::prob_dens_mdim_gaussian constructors.
  */
  template<class mat_t> class const_matrix_column_gen {

  protected:
    
    /// A reference to the original matrix
    const mat_t &m_;
    
    /// The selected column
    size_t column_;
    
  public:
    
    /// Create a column object from column \c column of matrix \c m 
  const_matrix_column_gen(const mat_t &m, size_t column) :
    m_(m), column_(column) {
    }
    
    /// Return a const reference to the ith row of the selected column
    const double &operator[](size_t i) const {
      return m_(i,column_);
    }
    
  };
  
  /** \brief Output the first \c n elements of a vector to a stream,
      \c os
      
      No trailing space is output after the last element, and an
      endline is output only if \c endline is set to \c true.  If the
      parameter \c n is zero, this function silently does nothing.

      This works with any class <tt>vec_t</tt> which has an operator[]
      which returns either the value of or a reference to the ith
      element and the element type has its own output operator which
      has been defined.
  */
  template<class vec_t> 
    void vector_out(std::ostream &os, size_t n, const vec_t &v, 
		    bool endline=false) {
    
    // This next line is important since n-1 is not well-defined if n=0
    if (n==0) {
      if (endline) os << std::endl;
      return;
    }

    for(size_t i=0;i<n-1;i++) os << v[i] << " ";
    os << v[n-1];
    if (endline) os << std::endl;
    return;
  }

  /** \brief Output a vector to a stream
      
      No trailing space is output after the last element, and an
      endline is output only if \c endline is set to \c true.  If the
      parameter \c n is zero, this function silently does nothing.

      This works with any class <tt>vec_t</tt> which has an operator[]
      which returns either the value of or a reference to the ith
      element and the element type has its own output operator which
      has been defined.
  */
  template<class vec_t> 
    void vector_out(std::ostream &os, const vec_t &v, bool endline=false) {
    
    size_t n=v.size();
    
    // This next line is important since n-1 is not well-defined if n=0
    if (n==0) return;

    for(size_t i=0;i<n-1;i++) os << v[i] << " ";
    os << v[n-1];
    if (endline) os << std::endl;
    return;
  }
  
  /** \brief Fill a vector with a specified grid
   */
  template<class vec_t, class data_t>
    void vector_grid(uniform_grid<data_t> g, vec_t &v) {
    g.template vector<vec_t>(v);
    return;
  }

  /// Set a matrix to unity on the diagonal and zero otherwise
  template<class mat_t> 
    void matrix_set_identity(size_t M, size_t N, mat_t &m) {
    for(size_t i=0;i<M;i++) {
      for(size_t j=0;j<N;j++) {
	if (i==j) m(i,j)=1.0;
	else m(i,j)=0.0;
      }
    }
    return;
  }

  /// Set a matrix to unity on the diagonal and zero otherwise
  template<class mat_t> 
    void matrix_set_identity(mat_t &m) {
    matrix_set_identity(m.size1(),m.size2(),m);
    return;
  }
  //@}

  /// \name Vector range classes and functions in src/base/vector.h
  //@{
  /** \brief Vector range function for pointers

      \note In this case, the return type is the same as the
      type of the first parameter. 
  */
  template<class dat_t> dat_t *vector_range
    (dat_t *v, size_t start, size_t last) {
    return v+start;
  }
  
  /** \brief Vector range function for const pointers

      \note In this case, the return type is the same as the
      type of the first parameter. 
  */
  template<class dat_t> const dat_t *const_vector_range
    (const dat_t *v, size_t start, size_t last) {
    return v+start;
  }
  
  /** \brief Vector range function template for ublas vectors

      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is not the same as the
      type of the first parameter. 
  */
  template<class dat_t> boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<dat_t> >
    vector_range(boost::numeric::ublas::vector<dat_t> &v,
		 size_t start, size_t last) {
    return boost::numeric::ublas::vector_range
      <boost::numeric::ublas::vector<dat_t> >
      (v,boost::numeric::ublas::range(start,last));
  }

  /** \brief Const vector range function template for ublas vectors

      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is not the same as the
      type of the first parameter. 
  */
  template<class dat_t> const boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<dat_t> >
    const_vector_range(boost::numeric::ublas::vector<dat_t> &v,
		       size_t start, size_t last) {
    return boost::numeric::ublas::vector_range
      <boost::numeric::ublas::vector<dat_t> >
      (v,boost::numeric::ublas::range(start,last));
  }
  
  /** \brief Const vector range function template for const ublas 
      vectors
      
      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is not the same as the
      type of the first parameter. 
  */
  template<class dat_t> const boost::numeric::ublas::vector_range
    <const boost::numeric::ublas::vector<dat_t> >
    const_vector_range(const boost::numeric::ublas::vector<dat_t> &v,
		       size_t start, size_t last) {
    return boost::numeric::ublas::vector_range
      <const boost::numeric::ublas::vector<dat_t> >
      (v,boost::numeric::ublas::range(start,last));
  }
  
  /** \brief Vector range function template for ublas vector
      ranges of ublas vectors

      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is not the same as the
      type of the first parameter. 
  */
  template<class dat_t>
    boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<dat_t> > >
    vector_range
    (boost::numeric::ublas::vector_range
     <boost::numeric::ublas::vector<dat_t> > &v,
     size_t start, size_t last) {
    return boost::numeric::ublas::vector_range
      <boost::numeric::ublas::vector_range
      <boost::numeric::ublas::vector<dat_t> > >
      (v,boost::numeric::ublas::range(start,last));
  }
  
  /** \brief Const vector range function template for ublas vector
      ranges of ublas vectors

      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is not the same as the
      type of the first parameter. 
  */
  template<class dat_t>
    const boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<dat_t> > >
    const_vector_range
    (boost::numeric::ublas::vector_range
     <boost::numeric::ublas::vector<dat_t> > &v,
     size_t start, size_t last) {
    return boost::numeric::ublas::vector_range
      <boost::numeric::ublas::vector_range
      <boost::numeric::ublas::vector<dat_t> > >
      (v,boost::numeric::ublas::range(start,last));
  }

  /** \brief Const vector range function template for const 
      ublas vector ranges of ublas vectors

      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is not the same as the
      type of the first parameter. 
  */
  template<class dat_t>
    const boost::numeric::ublas::vector_range
    <const boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<dat_t> > >
    const_vector_range
    (const boost::numeric::ublas::vector_range
     <boost::numeric::ublas::vector<dat_t> > &v,
     size_t start, size_t last) {
    return boost::numeric::ublas::vector_range
      <const boost::numeric::ublas::vector_range
      <boost::numeric::ublas::vector<dat_t> > >
      (v,boost::numeric::ublas::range(start,last));
  }

  /** \brief Const vector range function template for const 
      ublas vector ranges of const ublas vectors

      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is not the same as the
      type of the first parameter. 
  */
  template<class dat_t>
    const boost::numeric::ublas::vector_range
    <const boost::numeric::ublas::vector_range
    <const boost::numeric::ublas::vector<dat_t> > >
    const_vector_range
    (const boost::numeric::ublas::vector_range
     <const boost::numeric::ublas::vector<dat_t> > &v,
     size_t start, size_t last) {
    return boost::numeric::ublas::vector_range
      <const boost::numeric::ublas::vector_range
      <const boost::numeric::ublas::vector<dat_t> > >
      (v,boost::numeric::ublas::range(start,last));
  }

  // Forward definition for friendship
  template<class vec_t> class const_vector_range_gen;
  
  /** \brief Experimental vector range object
   */
  template<class vec_t> class vector_range_gen {
    
  protected:

    friend class const_vector_range_gen<vec_t>;

    /// A reference to the original vector
    vec_t &v_;

    /// The index offset
    size_t start_;

    /// The end() iterator
    size_t last_;
    
  public:
    
    /// Create an object starting with index \c start in vector \c v
  vector_range_gen(vec_t &v, size_t start, size_t last) : v_(v), 
      start_(start), last_(last) {
#if !O2SCL_NO_RANGE_CHECK
      if (last<start) {
	O2SCL_ERR2("End before beginning in vector_range_gen::",
		   "vector_range_gen(vec_t,size_t,size_t)",
		   o2scl::exc_einval);
      }
#endif
    }
    
    /// Create an object from a previously constructed range object
  vector_range_gen(const vector_range_gen &v2, size_t start,
		   size_t last) : v_(v2.v_), 
      start_(start+v2.start_), last_(last+v2.start_) {
#if !O2SCL_NO_RANGE_CHECK
      if (last<start) {
	O2SCL_ERR2("End before beginning in vector_range_gen::",
		   "vector_range_gen(vector_range_gen,size_t,size_t)",
		   o2scl::exc_einval);
      }
      if (last>v2.last_) {
	O2SCL_ERR2("End beyond end of previous vector in vector_range_gen::",
		   "vector_range_gen(vector_range_gen,size_t,size_t)",
		   o2scl::exc_einval);
      }
#endif
    }
      
    /// Return the vector size
    size_t size() const {
      return last_-start_;
    }
    
    /// Return a reference ith element
    double &operator[](size_t i) {
#if !O2SCL_NO_RANGE_CHECK
      if (i+start_>=last_) {
	O2SCL_ERR("Index out of range in vector_range_gen::operator[].",
		  o2scl::exc_einval);
      }
#endif
      return v_[i+start_];
    }
    
    /// Return a const reference ith element
    const double &operator[](size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
      if (i+start_>=last_) {
	O2SCL_ERR2("Index out of range in ",
		   "vector_range_gen::operator[] const.",o2scl::exc_einval);
      }
#endif
      return v_[i+start_];
    }
  };

  /** \brief Experimental const vector range object
   */
  template<class vec_t> class const_vector_range_gen {
    
  protected:

    /// A reference to the original vector
    const vec_t &v_;

    /// The index offset
    size_t start_;

    /// The end() iterator
    size_t last_;
    
  public:
    
    /// Create an object starting with index \c start in vector \c v
  const_vector_range_gen(const vec_t &v, size_t start, size_t last) : v_(v), 
      start_(start), last_(last) {
#if !O2SCL_NO_RANGE_CHECK
      if (last<start) {
	O2SCL_ERR2("End before beginning in vector_range_gen::",
		   "vector_range_gen(vec_t,size_t,size_t)",
		   o2scl::exc_einval);
      }
#endif
    }
    
    /// Create an object from a previously constructed range object
  const_vector_range_gen(const const_vector_range_gen &v2, size_t start,
			 size_t last) : v_(v2.v_), 
      start_(start+v2.start_), last_(last+v2.start_) {
#if !O2SCL_NO_RANGE_CHECK
      if (last<start) {
	O2SCL_ERR2("End before beginning in vector_range_gen::",
		   "vector_range_gen(vector_range_gen,size_t,size_t)",
		   o2scl::exc_einval);
      }
      if (last>v2.last_) {
	O2SCL_ERR2("End beyond end of previous vector in vector_range_gen::",
		   "vector_range_gen(vector_range_gen,size_t,size_t)",
		   o2scl::exc_einval);
      }
#endif
    }
      
    /// Create an object from a previously constructed range object
  const_vector_range_gen(const vector_range_gen<vec_t> &v2, size_t start,
			 size_t last) : v_(v2.v_), 
      start_(start+v2.start_), last_(last+v2.start_) {
#if !O2SCL_NO_RANGE_CHECK
      if (last<start) {
	O2SCL_ERR2("End before beginning in vector_range_gen::",
		   "vector_range_gen(vector_range_gen,size_t,size_t)",
		   o2scl::exc_einval);
      }
      if (last>v2.last_) {
	O2SCL_ERR2("End beyond end of previous vector in vector_range_gen::",
		   "vector_range_gen(vector_range_gen,size_t,size_t)",
		   o2scl::exc_einval);
      }
#endif
    }
      
    /// Return the vector size
    size_t size() const {
      return last_-start_;
    }
    
    /// Return a const reference ith element
    const double &operator[](size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
      if (i+start_>=last_) {
	O2SCL_ERR2("Index out of range in ",
		   "vector_range_gen::operator[] const.",o2scl::exc_einval);
      }
#endif
      return v_[i+start_];
    }
  };

  /** \brief Create a \ref o2scl::vector_range_gen object 
      from a <tt>std::vector</tt>
  */
  template<class data_t> vector_range_gen<std::vector<data_t> >
    vector_range(std::vector<data_t> &v, size_t start, size_t last) {
    return vector_range_gen<std::vector<data_t> >(v,start,last);
  }

  /** \brief Create a \ref o2scl::vector_range_gen object 
      from a <tt>std::vector</tt>
  */
  template<class data_t> const const_vector_range_gen<std::vector<data_t> >
    const_vector_range(const std::vector<data_t> &v, size_t start,
		       size_t last) {
    return const_vector_range_gen<std::vector<data_t> >(v,start,last);
  }
      
  /** \brief Create a \ref o2scl::vector_range_gen object 
      from a <tt>std::vector</tt>
  */
  template<class data_t> const const_vector_range_gen<std::vector<data_t> >
    const_vector_range(std::vector<data_t> &v, size_t start,
		       size_t last) {
    return const_vector_range_gen<std::vector<data_t> >(v,start,last);
  }
      
  /** \brief Recursively create a \ref o2scl::vector_range_gen object 
      from a vector range
  */
  template<class vec_t> vector_range_gen<vec_t>
    vector_range(vector_range_gen<vec_t> &v, size_t start, size_t last) {
    return vector_range_gen<vec_t>(v,start,last);
  }

  /** \brief Recursively create a const \ref o2scl::vector_range_gen
      object from a vector range
  */
  template<class vec_t> const const_vector_range_gen<vec_t>
    const_vector_range(vector_range_gen<vec_t> &v,
		       size_t start, size_t last) {
    return const_vector_range_gen<vec_t>(v,start,last);
  }

  /** \brief Recursively create a const \ref o2scl::vector_range_gen
      object from a const vector range
  */
  template<class vec_t> const const_vector_range_gen<vec_t>
    const_vector_range(const vector_range_gen<vec_t> &v,
		       size_t start, size_t last) {
    return const_vector_range_gen<vec_t>(v,start,last);
  }
  
  /** \brief Recursively create a const \ref o2scl::vector_range_gen
      object from a const vector range
  */
  template<class vec_t> const const_vector_range_gen<vec_t>
    const_vector_range(const const_vector_range_gen<vec_t> &v,
		       size_t start, size_t last) {
    return const_vector_range_gen<vec_t>(v,start,last);
  }
  
  /** \brief Vector range function template for <tt>std::vector</tt>
      
      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is the same as the
      type of the first parameter. 
      \note Unlike the ublas and pointer cases, this forces
      a copy. 
  */
  template<class dat_t> std::vector<dat_t>
    vector_range_copy(const std::vector<dat_t> &v, size_t start, size_t last) {
    return std::vector<dat_t> (v.begin()+start,v.begin()+last);
  }

  /** \brief Const vector range function template for <tt>std::vector</tt>
      
      The element with index \c start in the original vector
      will become the first argument in the new vector, and
      the new vector will have size <tt>last-start</tt> .

      \note In this case, the return type is the same as the
      type of the first parameter. 
      \note Unlike the ublas and pointer cases, this forces
      a copy. 
  */
  template<class dat_t> const std::vector<dat_t>
    vector_range_copy(const std::vector<dat_t> &v, size_t start,
		      size_t last) {
    return std::vector<dat_t> (v.begin()+start,v.begin()+last);
  }

  //@}
  
}

#if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)

#if defined (O2SCL_ARMA) || defined (DOXYGEN)
#include <armadillo>
namespace o2scl {

  /// \name Armadillo specializations in src/base/vector.h
  //@{
  /// Armadillo version of \ref matrix_max()
  double matrix_max(const arma::mat &data);

  /// Armadillo version of \ref matrix_min()
  double matrix_min(const arma::mat &data);

  /// Armadillo version of \ref matrix_row()
  template<> arma::subview_row<double>  
    matrix_row<arma::mat,arma::subview_row<double> >
    (arma::mat &M, size_t row);

  /// Armadillo version of \ref matrix_column()
  template<> arma::subview_col<double>  
    matrix_column<arma::mat,arma::subview_col<double> >
    (arma::mat &M, size_t column);
  //@}

}

#endif

#if defined (O2SCL_EIGEN) || defined (DOXYGEN)
#include <eigen3/Eigen/Dense>

namespace o2scl {

  /// \name Eigen specializations in src/base/vector.h
  //@{
  /// Eigen version of \ref matrix_max()
  double matrix_max(const Eigen::MatrixXd &data);

  /// Eigen version of \ref matrix_min()
  double matrix_min(const Eigen::MatrixXd &data);

  /// Eigen version of \ref matrix_row()
  template<> Eigen::MatrixXd::RowXpr 
    matrix_row<Eigen::MatrixXd,Eigen::MatrixXd::RowXpr>
    (Eigen::MatrixXd &M, size_t row);

  /// Eigen version of \ref matrix_column()
  template<> Eigen::MatrixXd::ColXpr 
    matrix_column<Eigen::MatrixXd,Eigen::MatrixXd::ColXpr>
    (Eigen::MatrixXd &M, size_t column);
  //@}

}

#endif

#else

#include <o2scl/vector_special.h>

// End of "#if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)"
#endif

#ifdef DOXYGEN
/** \brief Placeholder documentation of some related Boost objects
 */
namespace boost {
  /** \brief Documentation of Boost::numeric objects
   */
  namespace numeric {
    /** \brief Documentation of uBlas objects
     */
    namespace ublas {
      /** \brief The default vector type from uBlas 

	  The uBlas types aren't documented here, but the full documentation 
	  is available at
	  http://www.boost.org/doc/libs/release/libs/numeric/ublas/doc/index.htm

	  Internally in \o2, this is often typedef'd using
	  \code
	  typedef boost::numeric::ublas::vector<double> ubvector;
	  typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
	  typedef boost::numeric::ublas::vector<int> ubvector_int;
	  \endcode

	  This is documented in \ref vector.h .
      */
      template<class T, class A> class vector {
      };
      /** \brief The default matrix type from uBlas 
	  
	  The uBlas types aren't documented here, but the full documentation 
	  is available at
	  http://www.boost.org/doc/libs/release/libs/numeric/ublas/doc/index.htm

	  Internally in \o2, this is often typedef'd using
	  \code
	  typedef boost::numeric::ublas::matrix<double> ubmatrix;
	  typedef boost::numeric::ublas::matrix<size_t> ubmatrix_size_t;
	  typedef boost::numeric::ublas::matrix<int> ubmatrix_int;
	  \endcode

	  This is documented in \ref vector.h .
      */
      template<class T, class F, class A> class matrix {
      };
    }
  }
}
// End of "#ifdef DOXYGEN"
#endif

// End of "#ifndef O2SCL_VECTOR_H"
#endif
