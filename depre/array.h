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
#ifndef O2SCL_ARRAY_H
#define O2SCL_ARRAY_H

/** \file array.h
    \brief Various array classes

    For a more general discussion of vectors and matrices in \o2,
    see the \ref vecmat_section
    of the User's Guide.

    This file contains classes and functions for operating with
    C-style 1- or 2-dimensional arrays and pointers to double. For an
    example of the usage of the array allocation classes, see the \ref
    ex_mroot_sect . For more generic operations on generic vector
    objects (including in some cases C-style arrays), see also the
    file \ref vector.h . 

    This file contains the allocation classes
    - array_alloc
    - array_2d_alloc
    - pointer_alloc
    - pointer_2d_alloc

    the classes for the manipulation of arrays in smart_interp
    - array_reverse
    - array_subvector
    - array_subvector_reverse
    - array_const_reverse
    - array_const_subvector
    - array_const_subvector_reverse 

    the array equivalent of omatrix_row and omatrix_col (see 
    usage in <tt>src/ode/ode_it_solve_ts.cpp</tt>)
    - array_2d_row
    - array_2d_col

    \note The classes 
    - array_reverse
    - array_subvector
    - array_subvector_reverse
    - array_const_reverse
    - array_const_subvector
    - array_const_subvector_reverse 

    can be used with pointers or arrays, but array_alloc and
    pointer_alloc are \e not interchangable.

    \future Create a class which views a C-style array
    as a matrix and offers an <tt>operator(,)</tt>

*/
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <o2scl/err_hnd.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_sort.h>

#ifndef DOXYGENP
namespace o2scl
{
#endif

  /** \brief A simple class to provide an \c allocate() function
      for arrays
      
      The functions here are blank, as fixed-length arrays are
      automatically allocated and destroyed by the compiler. This
      class is present to provide an analog to \ref pointer_alloc and
      \ref ovector_alloc . 

      \comment
      This statement is commented out since interp_sma is deprecated
      at the moment.
      This class is used, for example, for \ref interp_sma .
      \endcomment

      \future Might it be possible to rework this so that 
      it does range checking and ensures that the user
      doesn't try to allocate more or less space? I.e.
      array_alloc<double[2]> complains if you try an allocate(x,3)?
  */
  template<class vec_t> class array_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(vec_t &v, size_t i) {}
    /// Free memory
    void free(vec_t &v) {}
  };

  /** \brief A simple class to provide an \c allocate() function
      for 2-dimensional arrays
      
      The functions here are blank, as fixed-length arrays are
      automatically allocated and destroyed by the compiler. This
      class is present to provide an analog to \ref pointer_2d_alloc
      and \ref omatrix_alloc. This class is used in
      <tt>gsl_mroot_hybrids_ts.cpp</tt> .
  */
  template<class mat_t> class array_2d_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(mat_t &v, size_t i, size_t j) {}
    /// Free memory
    void free(mat_t &v, size_t i) {}
  };

  /** \brief A simple class to provide an \c allocate() function
      for pointers
      
      This class uses \c new and \c delete and may throw a C++
      exception in the usual way. It contains a simple garbage
      collection mechanism to assist in exception safety.
      
      \future There may be a slight issue here if the pointer
      allocation succeeds and the list insertion fails, and the user
      catches the exception thrown by the list insertion failure, then
      a memory leak will ensue. This might be partially ameliorated by
      counting allocations and insertions. This failure mechanism is
      rarely encountered in practice.
  */
  template<class base_t> class pointer_alloc {
    
  protected:

    /// Store allocated pointers
    std::vector<base_t *> list;

  public:
    
    ~pointer_alloc() {
      typename std::vector<base_t *>::iterator it;
      for (it=list.begin();it!=list.end();it++) {
	base_t *v=*it;
	delete[] v;
      }
    }

    /// Allocate \c v for \c i elements
    void allocate(base_t *&v, size_t i) { 
      v=new base_t[i]; 
      list.push_back(v);
    }

    /// Free memory
    void free(base_t *&v) { 
      typename std::vector<base_t *>::iterator it;
      for (it=list.begin();it!=list.end();) {
	if (*it==v) {
	  list.erase(it);
	  it=list.end();
	} else {
	  it++;
	}
      }
      delete[] v; 
    }

  };

  /** \brief A simple class to provide an \c allocate() function
      for pointers

      This class uses \c new and \c delete and may throw a C++
      exception in the usual way. It contains a simple garbage
      collection mechanism to assist in exception safety.

      This class is used in \ref columnify and \ref cern_mroot.

      \future There may be a slight issue here if the pointer
      allocation succeeds and the list insertion fails, and the user
      catches the exception thrown by the list insertion failure, then
      a memory leak will ensue. This might be partially ameliorated by
      counting allocations and insertions. This failure mechanism is
      rarely encountered in practice.
  */
  template<class base_t> class pointer_2d_alloc {

  protected:

    typedef struct {
      /// Return \c p1<p2
      bool operator()(base_t ** const &p1, base_t ** const &p2) const {
	return p1<p2;
      }
    } pointer_comp;

    /// Store allocated pointers
    std::map<base_t **,size_t,pointer_comp> list;

    /// Iterator type
    typedef typename std::map<base_t **,size_t,pointer_comp>::iterator iter;

  public:

    ~pointer_2d_alloc() {
      for(iter it=list.begin();it!=list.end();it++) {
	base_t **v=it->first;
	for(size_t i=0;i<it->second;i++) {
	  delete[] v[i];
	}
	delete[] v;
      }
    }

    /// Allocate \c v for \c i elements
    void allocate(base_t **&v, size_t i, size_t j) { 
      v=new base_t *[i]; 
      size_t cnt=0;
      for(size_t m=0;m<i;m++) {
	v[m]=new base_t[j];
	cnt++;
      }
      list.insert(std::make_pair(v,cnt));
    }

    /// Free memory
    void free(base_t **&v, size_t i) { 
      iter it=list.find(v);
      size_t cnt=it->second;
      if (it!=list.end()) {
	list.erase(it);
      }
      for(size_t m=0;m<cnt;m++) delete[] v[m];
      delete[] v; 
    }
  };

  /** \brief A simple class which reverses the order of an array

      See an example of usage in <tt>interp_ts.cpp</tt> and
      <tt>smart_interp_ts.cpp</tt>.
   */
  template<size_t sz> class array_reverse {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The array pointer
    double *a;

#endif

  public:

    /// Create a reversed array from \c arr of size \c sz
    array_reverse(double *arr) {
      a=arr;
    }

    /** \brief Array-like indexing 
    */
    double &operator[](size_t i) {
      return a[sz-1-i];
    }
    
    /** \brief Array-like indexing 
    */
    const double &operator[](size_t i) const {
      return a[sz-1-i];
    }
    
  };

  /** \brief A simple class which reverses the order of an array

      See an example of usage in <tt>interp_ts.cpp</tt> and
      <tt>smart_interp_ts.cpp</tt>.
   */
  template<size_t sz> class array_const_reverse {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The array pointer
    double *a;

#endif

  public:

    /// Create a reversed array from \c arr of size \c sz
    array_const_reverse(const double *arr) {
      a=(double *)arr;
    }

    /** \brief Array-like indexing 
    */
    const double &operator[](size_t i) const {
      return a[sz-1-i];
    }
    
  };

  /** \brief A simple subvector class for an array (without error checking)
   */
  class array_subvector {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The array pointer
    double *a;

    /// The offset
    size_t off;

    /// The subvector length
    size_t len;

#endif

  public:

    /// Create a reversed array from \c arr of size \c sz
    array_subvector(double *arr, size_t offset, size_t n) {
      a=arr+offset;
      off=offset;
      len=n;
    }

    /** \brief Array-like indexing 
    */
    double &operator[](size_t i) {
      return a[i];
    }
    
    /** \brief Array-like indexing 
    */
    const double &operator[](size_t i) const {
      return a[i];
    }
    
  };

  /** \brief Column of a 2d array

      This works because two-dimensional arrays are always
      continguous (as indicated in appendix C of Soustroup's book).
  */
  template<size_t R, size_t C, class data_t=double> class array_2d_col {
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The array pointer
    data_t *a;
    
#endif
    
    public:
    
    /// Create an object as the <tt>i</tt>th column of \c mat
    array_2d_col(data_t mat[R][C], size_t i) {
      a=&(mat[0][i]);
    }
    
    /** \brief Array-like indexing 
    */
    data_t &operator[](size_t i) {
      return a[R*i];
    }
    
    /** \brief Array-like indexing 
    */
    const data_t &operator[](size_t i) const {
      return a[R*i];
    }
    
    };
  
  /** \brief Row of a 2d array
  */
  template<class array_2d_t, class data_t=double> class array_2d_row {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The array pointer
    data_t *a;

#endif

  public:

    /// Create an object as the <tt>i</tt>th row of \c mat
    array_2d_row(array_2d_t &mat, size_t i) {
      a=&(mat[i][0]);
    }

    /** \brief Array-like indexing, returns the element in the <tt>i</tt>th
	column of the chosen row
    */
    data_t &operator[](size_t i) {
      return a[i];
    }
    
    /** \brief Array-like indexing, returns the element in the <tt>i</tt>th
	column of the chosen row
    */
    const data_t &operator[](size_t i) const {
      return a[i];
    }
    
  };

  /** \brief A simple subvector class for a const array 
      (without error checking)

      \future Make the member data truly const and remove
      the extra typecast.
   */
  class array_const_subvector {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The array pointer
    double *a;

    /// The offset
    size_t off;

    /// The subvector length
    size_t len;

#endif

  public:

    /// Create a reversed array from \c arr of size \c sz
    array_const_subvector(const double *arr, size_t offset, size_t n) {
      a=(double *)arr+offset;
      off=offset;
      len=n;
    }

    /** \brief Array-like indexing 
    */
    const double &operator[](size_t i) const {
      return a[i];
    }
    
  };

  /** \brief Reverse a subvector of an array
   */
  class array_subvector_reverse {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The array pointer
    double *a;

    /// The offset
    size_t off;
    
    /// The subvector length
    size_t len;

#endif

  public:

    /// Create a reversed array from \c arr of size \c sz
    array_subvector_reverse(double *arr, size_t offset, size_t n) {
      a=arr+offset;
      off=offset;
      len=n;
    }

    /** \brief Array-like indexing 
    */
    double &operator[](size_t i) {
      return a[len-1-i];
    }
    
    /** \brief Array-like indexing 
    */
    const double &operator[](size_t i) const {
      return a[len-1-i];
    }
    
  };

  /** \brief Reverse a subvector of a const array
  */
  class array_const_subvector_reverse {
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The array pointer
    double *a;

    /// The offset
    size_t off;

    /// The subvector length
    size_t len;

#endif

  public:
    
    /// Create a reversed array from \c arr of size \c sz
    array_const_subvector_reverse(const double *arr, size_t offset, 
				  size_t n) {
      a=(double *)arr+offset;
      off=offset;
      len=n;
    }

    /** \brief Array-like indexing 
    */
    const double &operator[](size_t i) const {
      return a[len-1-i];
    }
    
  };

#ifndef DOXYGENP
}
#endif

#endif
