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
#ifndef O2SCL_OMATRIX_CX_TLATE_H
#define O2SCL_OMATRIX_CX_TLATE_H

/** \file omatrix_cx_tlate.h
    \brief File for definitions of complex matrices
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>

#include <o2scl/err_hnd.h>
#include <o2scl/ovector_tlate.h>
#include <o2scl/ovector_cx_tlate.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief A matrix view of double-precision numbers
  */
  template<class data_t, class parent_t, class block_t, class complex_t> 
    class omatrix_cx_view_tlate : 
  public parent_t {
    
  public:

    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same matrix
    omatrix_cx_view_tlate(const omatrix_cx_view_tlate &v) {
      parent_t::block=0;
      parent_t::data=v.data;
      parent_t::size1=v.size1;
      parent_t::size2=v.size2;
      parent_t::tda=v.tda;
      parent_t::owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    omatrix_cx_view_tlate& operator=(const omatrix_cx_view_tlate &v) {
      parent_t::block=0;
      parent_t::data=v.data;
      parent_t::size1=v.size1;
      parent_t::size2=v.size2;
      parent_t::tda=v.tda;
      parent_t::owner=0;      

      return *this;
    }
    //@}

    ~omatrix_cx_view_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    complex_t *operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1) {
	O2SCL_ERR((((std::string)"Array index ")+itos_nothrow(i)+
		   " out of bounds in omatrix_cx_view_tlate::"+
		   "operator[]. Size: "+itos_nothrow(parent_t::size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (complex_t *)(parent_t::data);
      }
#endif
      return (complex_t *)(parent_t::data+i*parent_t::tda*2);
    }
    
    /** \brief Array-like indexing 
    */
    const complex_t *operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1) {
	O2SCL_ERR((((std::string)"Array index ")+itos_nothrow(i)+
		   " out of bounds in omatrix_cx_view_tlate::oper"+
		   "ator[] const. Size: "+itos_nothrow(parent_t::size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (const complex_t *)(parent_t::data);
      }
#endif
      return (const complex_t *)(parent_t::data+i*parent_t::tda*2);
    }
    
    /** \brief Array-like indexing 
    */
    complex_t &operator()(size_t i, size_t j) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		   ", "+itos_nothrow(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::operator(). Sizes: "+
		   itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return *(complex_t *)(parent_t::data);
      }
#endif
      return *(complex_t *)(parent_t::data+(i*parent_t::tda+j)*2);
    }
    
    /** \brief Array-like indexing 
    */
    const complex_t &operator()(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		   ", "+itos_nothrow(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::operator() const. Sizes: "+
		   itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return *(const complex_t *)(parent_t::data);
      }
#endif
      return *(const complex_t *)(parent_t::data+(i*parent_t::tda+j)*2);
    }
    
    /** \brief Get (with optional range-checking) */
    complex_t get(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		 ", "+itos_nothrow(j)+" out of bounds"
		 +" in omatrix_cx_view_tlate::get(). Sizes: "+
		 itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		 " (indices should be less than sizes).").c_str(),
		gsl_eindex);
	complex_t zero={{0,0}};
	return zero;
      }
#endif
      return *(complex_t *)(parent_t::data+(i*parent_t::tda+j)*2);
    }
    
    /** \brief Get STL-like complex number (with optional range-checking) */
    std::complex<data_t> get_stl(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		 ", "+itos_nothrow(j)+" out of bounds"
		 +" in omatrix_cx_view_tlate::get_stl(). Sizes: "+
		 itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		 " (indices should be less than sizes).").c_str(),
		gsl_eindex);
	std::complex<data_t> zero=0;
	return zero;
      }
#endif
      return *(std::complex<data_t> *)(parent_t::data+(i*parent_t::tda+j)*2);
    }
    
    /** \brief Get real part (with optional range-checking) */
    data_t real(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		 ", "+itos_nothrow(j)+" out of bounds"
		 +" in omatrix_cx_view_tlate::real(). Sizes: "+
		 itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		 " (indices should be less than sizes).").c_str(),
		gsl_eindex);
	return 0;
      }
#endif
      return *(parent_t::data+(i*parent_t::tda+j)*2);
    }
    
    /** \brief Get imaginary part(with optional range-checking) */
    data_t imag(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		 ", "+itos_nothrow(j)+" out of bounds"
		 +" in omatrix_cx_view_tlate::imag(). Sizes: "+
		 itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		 " (indices should be less than sizes).").c_str(),
		gsl_eindex);
	return 0;
      }
#endif
      return *(parent_t::data+(i*parent_t::tda+j)*2+1);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    complex_t *get_ptr(size_t i, size_t j) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		 ", "+itos_nothrow(j)+" out of bounds"
		 +" in omatrix_cx_view_tlate::get_ptr(). Sizes: "+
		 itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		 " (indices should be less than sizes).").c_str(),
		gsl_eindex);
	return (complex_t *)(parent_t::data);
      }
#endif
      return (complex_t *)(parent_t::data+(i*parent_t::tda+j)*2);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const complex_t *get_const_ptr(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		 ", "+itos_nothrow(j)+" out of bounds"
		 +" in omatrix_cx_view_tlate::get_const_ptr(). Sizes: "+
		 itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		 " (indices should be less than sizes).").c_str(),
		gsl_eindex);
	return (const complex_t *)(parent_t::data);
      }
#endif
      return (const complex_t *)(parent_t::data+(i*parent_t::tda+j)*2);
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, size_t j, complex_t &val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		   ", "+itos_nothrow(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::set(). Sizes: "+
		   itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
      }
#endif
      parent_t::data[(i*parent_t::tda+j)*2]=GSL_REAL(val);
      parent_t::data[(i*parent_t::tda+j)*2+1]=GSL_IMAG(val);
      return;
    }

    /** \brief Set (with optional range-checking) */
    void set(size_t i, size_t j, data_t vr, data_t vi) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		   ", "+itos_nothrow(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::set(). Sizes: "+
		   itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
      }
#endif
      parent_t::data[(i*parent_t::tda+j)*2]=vr;
      parent_t::data[(i*parent_t::tda+j)*2+1]=vi;
      return;
    }

    /** \brief Set (with optional range-checking) */
    void set_real(size_t i, size_t j, data_t vr) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		   ", "+itos_nothrow(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::set_real(). Sizes: "+
		   itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
      }
#endif
      parent_t::data[(i*parent_t::tda+j)*2]=vr;
      return;
    }

    /** \brief Set (with optional range-checking) */
    void set_imag(size_t i, size_t j, data_t vi) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=parent_t::size1 || j>=parent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+itos_nothrow(i)+
		   ", "+itos_nothrow(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::set_imag(). Sizes: "+
		   itos_nothrow(parent_t::size1)+","+
		   itos_nothrow(parent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
      }
#endif
      parent_t::data[(i*parent_t::tda+j)*2+1]=vi;
      return;
    }

    /** \brief Set all*/
    void set_all(complex_t &val) {
      for(size_t i=0;i<parent_t::size1;i++) {
	for(size_t j=0;j<parent_t::size2;j++) {
	  parent_t::data[(i*parent_t::tda+j)*2]=GSL_REAL(val);
	  parent_t::data[(i*parent_t::tda+j)*2+1]=GSL_IMAG(val);
	}
      }
      return;
    }

    /** \brief Method to return number of rows
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t rows() const {
      return parent_t::size1;
    }

    /** \brief Method to return number of columns
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t cols() const {
      return parent_t::size2;
    }

    /** \brief Method to return matrix tda 
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t tda() const {
      return parent_t::tda;
    }
    //@}
    
    /// \name Other methods
    //@{
    /// Return true if this object owns the data it refers to
    bool is_owner() const {
      if (parent_t::owner==1) return true;
      return false;
    }
    //@}
    
#ifndef DOXYGENP

  protected:

    /** \brief Empty constructor provided for use by 
	omatrix_cx_tlate(const omatrix_cx_tlate &v)
    */
    omatrix_cx_view_tlate() {};

#endif
    
  };

  /** \brief A matrix of double-precision numbers
  
  */
  template<class data_t, class parent_t, class block_t, class complex_t> 
    class omatrix_cx_tlate : 
  public omatrix_cx_view_tlate<data_t,parent_t,block_t,complex_t> {

  public:

    /// \name Standard constructor 
    //@{
    /** \brief Create an omatrix of size \c n with owner as 'true' */
    omatrix_cx_tlate(size_t r=0, size_t c=0) {

      parent_t::data=0;

      // This must be set to 1 even if n=0 so that future
      // calls to operator= work properly
      parent_t::owner=1;
      
      if (r>0 && c>0) {
	parent_t::block=(block_t *)malloc(sizeof(block_t));
	if (parent_t::block) {
	  parent_t::block->data=(data_t *)malloc(2*r*c*sizeof(data_t));
	  if (parent_t::block->data) {
	    parent_t::block->size=r*c;
	    parent_t::data=parent_t::block->data;
	    parent_t::size1=r;
	    parent_t::size2=c;
	    parent_t::tda=c;
	  } else {
	    std::free(parent_t::block);
	    O2SCL_ERR("No memory for data in omatrix_cx_tlate constructor",
		    gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in omatrix_cx_tlate contructor",
		  gsl_enomem);
	}
      } else {
	parent_t::size1=0;
	parent_t::size2=0;
	parent_t::tda=0;
      }
    }
    //@}
    
    /// \name Copy constructors
    //@{
    /// Deep copy constructor, allocate new space and make a copy
    omatrix_cx_tlate(const omatrix_cx_tlate &v) : 
      omatrix_cx_view_tlate<data_t,parent_t,block_t,complex_t>() {
      size_t n=v.size1;
      size_t n2=v.size2;
      if (n>0 && n2>0) {
	parent_t::block=(block_t *)malloc(sizeof(block_t));
	if (parent_t::block) {
	  parent_t::block->data=(data_t *)malloc(2*n*n2*sizeof(data_t));
	  if (parent_t::block->data) {
	    parent_t::block->size=n*n2;
	    parent_t::data=parent_t::block->data;
	    parent_t::size1=n;
	    parent_t::size2=n2;
	    parent_t::tda=n2;
	    parent_t::owner=1;
	    for(size_t i=0;i<n;i++) {
	      for(size_t j=0;j<n2;j++) {
		*(parent_t::data+i*parent_t::tda+j)=*(v.data+i*v.tda+j);
		*(parent_t::data+i*parent_t::tda+j+1)=*(v.data+i*v.tda+j+1);
	      }
	    }
	  } else {
	    std::free(parent_t::block);
	    O2SCL_ERR("No memory for data in omatrix_cx_tlate constructor",
		    gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in omatrix_cx_tlate contructor",
		  gsl_enomem);
	}
      } else {
	parent_t::size1=0;
	parent_t::size2=0;
	parent_t::tda=0;
      }
    }

#ifdef O2SCL_NEVER_DEFINED
  }
  {
#endif
    
    /// Deep copy constructor, allocate new space and make a copy
    omatrix_cx_tlate
      (const omatrix_cx_view_tlate<data_t,parent_t,block_t,complex_t> &v) : 
	omatrix_cx_view_tlate<data_t,parent_t,block_t,complex_t>() {
	size_t r=v.rows();
	size_t c=v.cols();
	if (r>0 && c>0) {
	  parent_t::block=(block_t *)malloc(sizeof(block_t));
	  if (parent_t::block) {
	    parent_t::block->data=(data_t *)malloc(2*r*c*sizeof(data_t));
	    if (parent_t::block->data) {
	      parent_t::block->size=r*c;
	      parent_t::data=parent_t::block->data;
	      parent_t::size1=v.size1;
	      parent_t::size2=v.size2;
	      parent_t::tda=v.tda;
	      parent_t::owner=1;
	      for(size_t i=0;i<r;i++) {
		for(size_t j=0;i<c;j++) {
		  *(parent_t::data+i*parent_t::tda+j)=*(v.data+i*v.tda+j);
		  *(parent_t::data+i*parent_t::tda+j+1)=
		    *(v.data+i*v.tda+j+1);
		}
	      }
	    } else {
	      std::free(parent_t::block);
	      O2SCL_ERR("No memory for data in omatrix_cx_tlate constructor",
		      gsl_enomem);
	    }
	  } else {
	    O2SCL_ERR("No memory for block in omatrix_cx_tlate contructor",
		    gsl_enomem);
	  }
	} else {
	  parent_t::size1=0;
	  parent_t::size2=0;
	  parent_t::tda=0;
	}
      }
      //@}
    
      ~omatrix_cx_tlate() {
	if (parent_t::size1>0) {
	  if (parent_t::owner==1) {
	    if (parent_t::block->size>0) {
	      std::free(parent_t::block->data);
	    }
	    std::free(parent_t::block);
	    parent_t::size1=0;
	    parent_t::size2=0;
	  }
	}
      }
    
      /// \name Memory allocation
      //@{
      /** \brief Allocate memory for size \c n after freeing any memory
	  presently in use
      */
      void allocate(size_t nrows, size_t ncols) {
	if (parent_t::size1>0 || parent_t::size2>0) free();
      
	if (nrows>0 && ncols>0) {
	  parent_t::block=(block_t *)malloc(sizeof(block_t));
	  if (parent_t::block) {
	    parent_t::block->data=(data_t *)
	      malloc(2*nrows*ncols*sizeof(data_t));
	    if (parent_t::block->data) {
	      parent_t::block->size=nrows*ncols;
	      parent_t::data=parent_t::block->data;
	      parent_t::size1=nrows;
	      parent_t::size2=ncols;
	      parent_t::tda=ncols;
	      parent_t::owner=1;
	    } else {
	      std::free(parent_t::block);
	      O2SCL_ERR("No memory for data in omatrix_cx_tlate::allocate()",
			gsl_enomem);
	    }
	  } else {
	    O2SCL_ERR("No memory for block in omatrix_cx_tlate::allocate()",
		      gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("Zero size in omatrix::allocate()",gsl_einval);
	}
	return;
      }

      /** \brief Free the memory 
    
	  This function will safely do nothing if used without first
	  allocating memory or if called multiple times in succession.
      */
      void free() {
	if (parent_t::size1>0) {
	  if (parent_t::owner==1) {
	    if (parent_t::block->size>0) {
	      std::free(parent_t::block->data);
	    }
	    std::free(parent_t::block);
	  }
	  parent_t::size1=0;
	  parent_t::size2=0;
	  parent_t::tda=0;
	}
	return;
      }
      //@}
    
      /// \name Other methods
      //@{
      /// \brief Compute the transpose
      omatrix_cx_tlate<data_t,parent_t,block_t,complex_t> transpose() {
	omatrix_cx_tlate<data_t,parent_t,block_t,complex_t> 
	  result(parent_t::size2,parent_t::size1);
	for(size_t i=0;i<parent_t::size1;i++) {
	  for(size_t j=0;j<parent_t::size2;j++) {
	    result[j][i]=(*this)[i][j];
	  }
	}
      }
    
      /// Compute the conjugate transpose
      omatrix_cx_tlate<data_t,parent_t,block_t,complex_t> htranspose() {
	omatrix_cx_tlate<data_t,parent_t,block_t,complex_t> 
	  result(parent_t::size2,parent_t::size1);
	for(size_t i=0;i<parent_t::size1;i++) {
	  for(size_t j=0;j<parent_t::size2;j++) {
	    result[j][i]=gsl_complex_conjugate((*this)[i][j]);
	  }
	}
      }
      //@}

  };

  /** \brief Create a vector from a row of a matrix 
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t,
    class complex_t> class omatrix_cx_row_tlate : 
  public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {
  public:
    /** \brief Create a vector from row \c i of matrix \c m */
    omatrix_cx_row_tlate
      (omatrix_cx_view_tlate<data_t,mparent_t,block_t,complex_t> &m, 
       size_t i) {
      if (i<m.size1) {
	vparent_t::size=m.size2;
	vparent_t::stride=1;
	vparent_t::data=m.data+m.tda()*i;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };
    
  /** \brief Create a vector from a row of a matrix 
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t,
    class complex_t> class omatrix_cx_const_row_tlate :
  public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {
  public:
    /** \brief Create a vector from row \c i of matrix \c m */
    omatrix_cx_const_row_tlate
      (const omatrix_cx_view_tlate<data_t,mparent_t,block_t,complex_t> &m, 
       size_t i) {
      if (i<m.size1) {
	vparent_t::size=m.size2;
	vparent_t::stride=1;
	vparent_t::data=m.data+m.tda*i;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };
    
  /** \brief Create a vector from a column of a matrix
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t, 
    class complex_t> class omatrix_cx_col_tlate :
  public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {
  public:
    /** \brief Create a vector from col \c i of matrix \c m */
    omatrix_cx_col_tlate
      (omatrix_cx_view_tlate<data_t,mparent_t,block_t,complex_t> &m, 
       size_t i) {
      if (i<m.size2) {
	vparent_t::size=m.size1;
	vparent_t::stride=m.tda();
	vparent_t::data=m.data+2*i;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };
    
  /** \brief Create a vector from a column of a matrix
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t,
    class complex_t> class omatrix_cx_const_col_tlate :
  public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {
  public:
    /** \brief Create a vector from col \c i of matrix \c m */
    omatrix_cx_const_col_tlate
      (omatrix_cx_view_tlate<data_t,mparent_t,block_t,complex_t> &m,
       size_t i) {
      if (i<m.size2) {
	vparent_t::size=m.size1;
	vparent_t::stride=m.tda();
	vparent_t::data=m.data+2*i;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };

  /// omatrix_cx typedef
  typedef omatrix_cx_tlate<double,gsl_matrix_complex,gsl_block_complex,
    gsl_complex> omatrix_cx;
  /// omatrix_cx_view typedef
  typedef omatrix_cx_view_tlate<double,gsl_matrix_complex,gsl_block_complex,
    gsl_complex> omatrix_cx_view;
  /// omatrix_cx_row typedef
  typedef omatrix_cx_row_tlate<double,gsl_matrix_complex,gsl_vector_complex,
    gsl_block_complex,gsl_complex> omatrix_cx_row;
  /// omatrix_cx_col typedef
  typedef omatrix_cx_col_tlate<double,gsl_matrix_complex,gsl_vector_complex,
    gsl_block_complex,gsl_complex> omatrix_cx_col;
  /// omatrix_cx_const_row typedef
  typedef omatrix_cx_const_row_tlate<double,gsl_matrix_complex,
    gsl_vector_complex,gsl_block_complex,gsl_complex> omatrix_cx_const_row;
  /// omatrix_cx_const_col typedef
  typedef omatrix_cx_const_col_tlate<double,gsl_matrix_complex,
    gsl_vector_complex,gsl_block_complex,gsl_complex> omatrix_cx_const_col;

#ifndef DOXYGENP
}
#endif

#endif

