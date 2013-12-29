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
#ifndef O2SCL_OMATRIX_TLATE_H
#define O2SCL_OMATRIX_TLATE_H

/** \file omatrix_tlate.h
    \brief Matrix classes

    \future The \c xmatrix class demonstrates how operator[] could
    return an ovector_array object and thus provide more
    bounds-checking. This would demand including a new parameter in
    omatrix_view_tlate which contains the vector type.
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_ieee_utils.h>

#include <o2scl/err_hnd.h>
#include <o2scl/ovector_tlate.h>
#include <o2scl/umatrix_tlate.h>
#include <o2scl/array.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief A const matrix view of omatrix objects
  */
  template<class data_t, class mparent_t, class block_t> 
    class omatrix_const_view_tlate : 
  public mparent_t {

  public:
    
    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same matrix
    omatrix_const_view_tlate(const omatrix_const_view_tlate &v) {
      mparent_t::block=0;
      mparent_t::data=v.data;
      mparent_t::size1=v.size1;
      mparent_t::size2=v.size2;
      mparent_t::tda=v.tda();
      mparent_t::owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    omatrix_const_view_tlate& operator=(const omatrix_const_view_tlate &v) {
      mparent_t::block=0;
      mparent_t::data=v.data;
      mparent_t::size1=v.size1;
      mparent_t::size2=v.size2;
      mparent_t::tda=v.tda;
      mparent_t::owner=0;      

      return *this;
    }
    //@}

    ~omatrix_const_view_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    const data_t *operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in omatrix_const_view_tlate::operator[] const. Size: "+
		   szttos(mparent_t::size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (mparent_t::data);
      }
#endif
      return mparent_t::data+i*mparent_t::tda;
    }
    
    /** \brief Array-like indexing 
    */
    const data_t &operator()(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::operator() const. Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return *(mparent_t::data);
      }
#endif
      return *(mparent_t::data+i*mparent_t::tda+j);
    }
    
    /** \brief Get (with optional range-checking) */
    data_t get(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::get() const. Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return *(mparent_t::data);
      }
#endif
      return *(mparent_t::data+i*mparent_t::tda+j);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i, size_t j) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::get_ptr(). Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return (mparent_t::data);
      }
#endif
      return mparent_t::data+i*mparent_t::tda+j;
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const data_t *get_const_ptr(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::get_const_ptr(). Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return (mparent_t::data);
      }
#endif
      return (const data_t *)(mparent_t::data+i*mparent_t::tda+j);
    }
    
    /** \brief Method to return number of rows
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t rows() const {
      return mparent_t::size1;
    }

    /** \brief Method to return number of columns
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t cols() const {
      return mparent_t::size2;
    }

    /** \brief Method to return matrix tda 
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t tda() const {
      return mparent_t::tda;
    }
    //@}

    /// \name Other methods
    //@{
    /** \brief Return true if this object owns the data it refers to
	
	This can be used to determine if an object is a "matrix_view",
	or a "matrix". If is_owner() is true, then it is an \ref
	omatrix_tlate object.
    */
    bool is_owner() const {
      if (mparent_t::owner==1) return true;
      return false;
    }

    /// The largest matrix element
    data_t max() const {
      data_t mx=*mparent_t::data;
      for(size_t i=0;i<mparent_t::size1;i++) {
	for(size_t j=0;j<mparent_t::size2;j++) {
	  if (*(mparent_t::data+i*mparent_t::tda+j)>mx) {
	    mx=*(mparent_t::data+i*mparent_t::tda+j);
	  }
	}
      }
      return mx;
    }
    
    /// The smallest matrix element
    data_t min() const {
      data_t mx=*mparent_t::data;
      for(size_t i=0;i<mparent_t::size1;i++) {
	for(size_t j=0;j<mparent_t::size2;j++) {
	  if (*(mparent_t::data+i*mparent_t::tda+j)<mx) {
	    mx=*(mparent_t::data+i*mparent_t::tda+j);
	  }
	}
      }
      return mx;
    }
    
    /** \brief Return a \c const gsl matrix */
    const mparent_t *get_gsl_matrix_const() const { return this; };
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief The default constructor
     */
    omatrix_const_view_tlate() {};

#endif
    
  };

  /** \brief A base class for omatrix and omatrix_view
  */
  template<class data_t, class mparent_t, class block_t> 
    class omatrix_base_tlate : 
  public omatrix_const_view_tlate<data_t,mparent_t,block_t> {
    
  public:
    
    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same matrix
    omatrix_base_tlate(const omatrix_base_tlate &v) {
      mparent_t::block=0;
      mparent_t::data=v.data;
      mparent_t::size1=v.size1;
      mparent_t::size2=v.size2;
      mparent_t::tda=v.tda();
      mparent_t::owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    omatrix_base_tlate& operator=(const omatrix_base_tlate &v) {
      mparent_t::block=0;
      mparent_t::data=v.data;
      mparent_t::size1=v.size1;
      mparent_t::size2=v.size2;
      mparent_t::tda=v.tda;
      mparent_t::owner=0;      

      return *this;
    }
    //@}

    ~omatrix_base_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    data_t *operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in omatrix_base_tlate::operator[]. Size: "+
		   szttos(mparent_t::size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (mparent_t::data);
      }
#endif
      return mparent_t::data+i*mparent_t::tda;
    }
    
    /** \brief Array-like indexing 
    */
    const data_t *operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in omatrix_base_tlate::operator[] const. Size: "+
		   szttos(mparent_t::size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (mparent_t::data);
      }
#endif
      return mparent_t::data+i*mparent_t::tda;
    }
    
    /** \brief Array-like indexing 
    */
    data_t &operator()(size_t i, size_t j) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::operator(). Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return *(mparent_t::data);
      }
#endif
      return *(mparent_t::data+i*mparent_t::tda+j);
    }
    
    /** \brief Array-like indexing 
    */
    const data_t &operator()(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::operator() const. Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return *(mparent_t::data);
      }
#endif
      return *(mparent_t::data+i*mparent_t::tda+j);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i, size_t j) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::get_ptr(). Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return (mparent_t::data);
      }
#endif
      return mparent_t::data+i*mparent_t::tda+j;
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, size_t j, data_t val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::set(). Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
      }
#endif
      *(mparent_t::data+i*mparent_t::tda+j)=val;
      return;
    }

    /** \brief Set all of the value to be the value \c val */
    void set_all(data_t val) {
      for(size_t i=0;i<mparent_t::size1;i++) {
	for(size_t j=0;j<mparent_t::size2;j++) {
	  *(mparent_t::data+i*mparent_t::tda+j)=val;
	}
      }
      return;
    }
    //@}

    /// \name Other methods
    //@{
    /** \brief Return a gsl matrix */
    mparent_t *get_gsl_matrix() { return this; };
    //@}

    /// \name Arithmetic 
    //@{
    /** \brief operator+= */
    omatrix_base_tlate<data_t,mparent_t,block_t> &operator+=
      (const omatrix_base_tlate<data_t,mparent_t,block_t> &x) {
      size_t lsize=x.size1;
      if (lsize>mparent_t::size1) lsize=mparent_t::size1;
      size_t lsize2=x.size2;
      if (lsize2>mparent_t::size2) lsize2=mparent_t::size2;
      for(size_t i=0;i<lsize;i++) {
	for(size_t j=0;j<lsize2;j++) {
	  (*this)[i][j]+=x[i][j];
	}
      }
      
      return *this;
    }
    
    /** \brief operator-= */
    omatrix_base_tlate<data_t,mparent_t,block_t> &operator-=
      (const omatrix_base_tlate<data_t,mparent_t,block_t> &x) {
      size_t lsize=x.size1;
      if (lsize>mparent_t::size1) lsize=mparent_t::size1;
      size_t lsize2=x.size2;
      if (lsize2>mparent_t::size2) lsize2=mparent_t::size2;
      for(size_t i=0;i<lsize;i++) {
	for(size_t j=0;j<lsize2;j++) {
	  (*this)[i][j]+=x[i][j];
	}
      }

      return *this;
    }
    
    /** \brief operator+= */
    omatrix_base_tlate<data_t,mparent_t,block_t> 
      &operator+=(const data_t &y) {
      for(size_t i=0;i<mparent_t::size1;i++) {
	for(size_t j=0;j<mparent_t::size2;j++) {
	  (*this)[i][j]+=y;
	}
      }
      
      return *this;
    }

    /** \brief operator-= */
    omatrix_base_tlate<data_t,mparent_t,block_t> 
      &operator-=(const data_t &y) {
      for(size_t i=0;i<mparent_t::size1;i++) {
	for(size_t j=0;j<mparent_t::size2;j++) {
	  (*this)[i][j]-=y;
	}
      }
      
      return *this;
    }

    /** \brief operator*= */
    omatrix_base_tlate<data_t,mparent_t,block_t> 
      &operator*=(const data_t &y) {
      for(size_t i=0;i<mparent_t::size1;i++) {
	for(size_t j=0;j<mparent_t::size2;j++) {
	  (*this)[i][j]*=y;
	}
      }
      
      return *this;
    }
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief The default constructor
     */
    omatrix_base_tlate() {};

#endif
    
  };

  /** \brief A matrix view of double-precision numbers

      This is a matrix view, which views a matrix stored somewhere
      else. For a full matrix template class with its own memory, see
      \ref omatrix_tlate . The basic matrix classes are built upon
      these templates. A matrix of double-precision numbers is an
      object of type \ref omatrix . See \ref vecmat_section in the
      User's Guide for more information.
  */
  template<class data_t, class mparent_t, class block_t> 
    class omatrix_view_tlate : 
  public omatrix_base_tlate<data_t,mparent_t,block_t> {

  public:

    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same matrix
  omatrix_view_tlate(const omatrix_view_tlate &v) : 
    omatrix_base_tlate<data_t,mparent_t,block_t>() {
      mparent_t::block=0;
      mparent_t::data=v.data;
      mparent_t::size1=v.size1;
      mparent_t::size2=v.size2;
      mparent_t::tda=v.tda();
      mparent_t::owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    omatrix_view_tlate& operator=(const omatrix_view_tlate &v) {
      mparent_t::block=0;
      mparent_t::data=v.data;
      mparent_t::size1=v.size1;
      mparent_t::size2=v.size2;
      mparent_t::tda=v.tda;
      mparent_t::owner=0;      

      return *this;
    }

    /// Shallow copy constructor - create a new view of the same matrix
  omatrix_view_tlate(omatrix_base_tlate<data_t,mparent_t,block_t> &v) :
    omatrix_base_tlate<data_t,mparent_t,block_t>() {
      mparent_t::block=0;
      mparent_t::data=v.data;
      mparent_t::size1=v.size1;
      mparent_t::size2=v.size2;
      mparent_t::tda=v.tda();
      mparent_t::owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    omatrix_view_tlate& operator=
      (omatrix_base_tlate<data_t,mparent_t,block_t> &v) {
      mparent_t::block=0;
      mparent_t::data=v.data;
      mparent_t::size1=v.size1;
      mparent_t::size2=v.size2;
      mparent_t::tda=v.tda;
      mparent_t::owner=0;      

      return *this;
    }
    //@}

    ~omatrix_view_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    data_t *operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in omatrix_view_tlate::operator[]. Size: "+
		   szttos(mparent_t::size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (mparent_t::data);
      }
#endif
      return mparent_t::data+i*mparent_t::tda;
    }
    
    /** \brief Array-like indexing 
    */
    data_t &operator()(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::operator(). Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return *(mparent_t::data);
      }
#endif
      return *(mparent_t::data+i*mparent_t::tda+j);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::get_ptr(). Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
	return (mparent_t::data);
      }
#endif
      return mparent_t::data+i*mparent_t::tda+j;
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, size_t j, data_t val) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=mparent_t::size1 || j>=mparent_t::size2) {
	O2SCL_ERR((((std::string)"Array indices ")+szttos(i)+
		   ", "+szttos(j)+" out of bounds"
		   +" in omatrix_cx_view_tlate::set(). Sizes: "+
		   szttos(mparent_t::size1)+","+szttos(mparent_t::size2)+
		   " (indices should be less than sizes).").c_str(),
		  gsl_eindex);
      }
#endif
      *(mparent_t::data+i*mparent_t::tda+j)=val;
      return;
    }

    /** \brief Set all of the value to be the value \c val */
    void set_all(double val) const {
      for(size_t i=0;i<mparent_t::size1;i++) {
	for(size_t j=0;j<mparent_t::size2;j++) {
	  *(mparent_t::data+i*mparent_t::tda+j)=val;
	}
      }
      return;
    }
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief The default constructor
     */
    omatrix_view_tlate() {};

#endif

  };

  /** \brief A matrix of double-precision numbers
  
      The basic matrix classes are built upon this template. A matrix
      of double-precision numbers is an object of type \ref omatrix ,
      which is just a <tt>typedef</tt> defined using this class
      template. See \ref vecmat_section in the User's Guide for more
      information.
  */
  template<class data_t, class mparent_t, class vparent_t, class block_t> 
    class omatrix_tlate : 
  public omatrix_base_tlate<data_t,mparent_t,block_t> {
  public:
    
    /// \name Standard constructor
    //@{
    /** \brief Create an omatrix of size \c n with owner as \c true.
     */
    omatrix_tlate(size_t r=0, size_t c=0) {

      mparent_t::block=0;
      mparent_t::data=0;
      mparent_t::size1=0;
      mparent_t::size2=0;
      mparent_t::tda=0;

      // This must be set to 1 even if n=0 so that future
      // calls to operator= work properly
      mparent_t::owner=1;

      if (r>0 && c>0) {
	mparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (mparent_t::block) {
	  mparent_t::block->data=(data_t *)malloc(r*c*sizeof(data_t));
	  if (mparent_t::block->data) {
	    mparent_t::block->size=r*c;
	    mparent_t::data=mparent_t::block->data;
	    mparent_t::size1=r;
	    mparent_t::size2=c;
	    mparent_t::tda=c;
	  } else {
	    std::free(mparent_t::block);
	    O2SCL_ERR("No memory for data in omatrix_tlate constructor",
		      gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in omatrix_tlate contructor",
		    gsl_enomem);
	}
      }
    }
    //@}
    
    /// \name Copy constructors
    //@{
    /// Deep copy constructor, allocate new space and make a copy
  omatrix_tlate(const omatrix_tlate &v) : 
    omatrix_base_tlate<data_t,mparent_t,block_t>() {

      mparent_t::block=0;
      mparent_t::data=0;
      mparent_t::size1=0;
      mparent_t::size2=0;
      mparent_t::tda=0;
      mparent_t::owner=1;

      size_t n=v.size1;
      size_t n2=v.size2;

      if (n>0 && n2>0) {
	mparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (mparent_t::block) {
	  mparent_t::block->data=(data_t *)malloc(n*n2*sizeof(data_t));
	  if (mparent_t::block->data) {
	    mparent_t::block->size=n*n2;
	    mparent_t::data=mparent_t::block->data;
	    mparent_t::size1=n;
	    mparent_t::size2=n2;
	    mparent_t::tda=n2;
	    mparent_t::owner=1;
	    for(size_t i=0;i<n;i++) {
	      for(size_t j=0;j<n2;j++) {
		*(mparent_t::data+i*mparent_t::tda+j)=v[i][j];
	      }
	    }
	  } else {
	    std::free(mparent_t::block);
	    O2SCL_ERR("No memory for data in omatrix_tlate constructor",
		      gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in omatrix_tlate contructor",
		    gsl_enomem);
	}
      } else {
	mparent_t::size1=0;
	mparent_t::size2=0;
	mparent_t::tda=0;
      }
    }
    
    /// Deep copy constructor, allocate new space and make a copy
    omatrix_tlate
      (const omatrix_const_view_tlate<data_t,mparent_t,block_t> &v) : 
    omatrix_view_tlate<data_t,mparent_t,block_t>() {

      mparent_t::block=0;
      mparent_t::data=0;
      mparent_t::size1=0;
      mparent_t::size2=0;
      mparent_t::tda=0;
      mparent_t::owner=1;

      size_t r=v.rows();
      size_t c=v.cols();
      
      if (r>0 && c>0) {
	mparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (mparent_t::block) {
	  mparent_t::block->data=(data_t *)malloc(r*c*sizeof(data_t));
	  if (mparent_t::block->data) {
	    mparent_t::block->size=r*c;
	    mparent_t::data=mparent_t::block->data;
	    mparent_t::size1=v.size1;
	    mparent_t::size2=v.size2;
	    mparent_t::tda=v.tda;
	    mparent_t::owner=1;
	    for(size_t i=0;i<r;i++) {
	      for(size_t j=0;i<c;j++) {
		*(mparent_t::data+i*mparent_t::tda+j)=v[i][j];
	      }
	    }
	  } else {
	    std::free(mparent_t::block);
	    O2SCL_ERR("No memory for data in omatrix_tlate constructor",
		      gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in omatrix_tlate contructor",
		    gsl_enomem);
	}
      }
    }
    
    /** \brief Deep copy constructor, if owner is true, allocate
	space and make a new copy, otherwise, just copy into the
	view
    */
    omatrix_tlate& operator=(const omatrix_tlate &v) {
	  
      // Check for self-assignment
      if (this==&v) return *this;
      
      mparent_t::block=0;
      mparent_t::data=0;
      mparent_t::size1=0;
      mparent_t::size2=0;
      mparent_t::tda=0;
      mparent_t::owner=1;

      size_t sz=v.rows();
      size_t sz2=v.cols();

      // Changed here
      if (sz>0 && sz2>0) {
	if (mparent_t::owner) {
	  allocate(sz,sz2);
	} else {
	  if (mparent_t::size1!=sz || mparent_t::size2!=sz2) {
	    O2SCL_ERR("Sizes don't match in omatrix_tlate::operator=()",
		      gsl_ebadlen);
	    return *this;
	  }
	}
	for(size_t i=0;i<sz;i++) {
	  for(size_t j=0;j<sz2;j++) {
	    *(mparent_t::data+i*mparent_t::tda+j)=v[i][j];
	  }
	}
      }
      return *this;
    }

    /** \brief Deep copy constructor, if owner is true, allocate
	space and make a new copy, otherwise, just copy into the
	view
    */
    omatrix_tlate& operator=
      (const omatrix_const_view_tlate<data_t,mparent_t,block_t> &v) {

      mparent_t::block=0;
      mparent_t::data=0;
      mparent_t::size1=0;
      mparent_t::size2=0;
      mparent_t::tda=0;
      mparent_t::owner=1;
      
      // Check for self-assignment
      if (this==&v) return *this;

      size_t r=v.rows();
      size_t c=v.cols();

      mparent_t::data=0;
      if (mparent_t::owner) {
	allocate(r,c);
      } else {
	if (mparent_t::size1!=r || mparent_t::size2!=c) {
	  O2SCL_ERR("Sizes don't match in omatrix_tlate::operator=()",
		    gsl_ebadlen);
	  return *this;
	}
      }
      for(size_t i=0;i<r;i++) {
	for(size_t j=0;j<c;j++) {
	  *(mparent_t::data+i*mparent_t::tda+j)=v[i][j];
	}
      }
      return *this;
    }
    
    /** \brief Deep copy from an array of ovectors
     */
    omatrix_tlate(size_t n,
		  ovector_base_tlate<data_t,vparent_t,block_t> ova[]) {
      if (n>0) {
	size_t n2=ova[0];
	if (n2>0) {
	  allocate(n,n2);
	  for(size_t i=0;i<n;i++) {
	    for(size_t j=0;j<n2;j++) {
	      (*this)[i][j]=ova[i][j];
	    }
	  }
	}
      }
    }

    /** \brief Deep copy from an array of uvectors
     */
    omatrix_tlate(size_t n, uvector_base_tlate<data_t> uva[]) {
      if (n>0) {
	size_t n2=uva[0];
	if (n2>0) {
	  allocate(n,n2);
	  for(size_t i=0;i<n;i++) {
	    for(size_t j=0;j<n2;j++) {
	      (*this)[i][j]=uva[i][j];
	    }
	  }
	}
      }
    }

    /** \brief Deep copy from a C-style 2-d array
     */
    omatrix_tlate(size_t n, size_t n2, data_t **csa) {
      if (n>0 && n2>0) {
	allocate(n,n2);
	for(size_t i=0;i<n;i++) {
	  for(size_t j=0;j<n2;j++) {
	    (*this)[i][j]=csa[i][j];
	  }
	}
      }
    }

    //@}
    
    ~omatrix_tlate() {
      if (mparent_t::size1>0) {
	if (mparent_t::owner==1) {
	  if (mparent_t::block->size>0) {
	    std::free(mparent_t::block->data);
	  }
	  std::free(mparent_t::block);
	  mparent_t::size1=0;
	  mparent_t::size2=0;
	}
      }
    }

    /// \name Memory allocation
    //@{
    void resize(size_t new_rows, size_t new_cols) {
      // New size is zero so just free()
      if (new_rows==0 && new_cols==0) {
	free();
	return;
      }
      // Old size was zero so just allocate()
      if (mparent_t::size1==0 && mparent_t::size2==0) {
	allocate(new_rows,new_cols);
	return;
      }
      // If the sizes are the same, do nothing
      if (new_rows==mparent_t::size1 && new_cols==mparent_t::size2) {
	return;
      }
      // Neither size was zero, so allocate
      block_t *new_block=(block_t *)malloc(sizeof(block_t));
      if (new_block==0) {
	O2SCL_ERR2("Ran out of memory for block in ",
		   "omatrix_tlate::resize().",gsl_enomem);
      }
      new_block->data=(data_t *)malloc(new_rows*new_cols*sizeof(data_t));
      if (new_block->data==0) {
	std::free(new_block);
	O2SCL_ERR2("Ran out of memory for data in ",
		   "omatrix_tlate::resize().",gsl_enomem);
      }
      new_block->size=new_rows*new_cols;
      // Perform copy
      size_t min_rows, min_cols;
      if (new_rows<mparent_t::size1) min_rows=new_rows;
      else min_rows=mparent_t::size1;
      if (new_cols<mparent_t::size2) min_cols=new_cols;
      else min_cols=mparent_t::size2;
      for(size_t i=0;i<min_rows;i++) {
	for(size_t j=0;j<min_cols;j++) {
	  *(new_block->data+i*new_cols+j)=(*this)[i][j];
	}
      }
      // Free
      free();
      // Set new array
      mparent_t::block=new_block;
      mparent_t::data=mparent_t::block->data;
      mparent_t::size1=new_rows;
      mparent_t::size2=new_cols;
      mparent_t::tda=new_cols;
      mparent_t::owner=1;
      return;
    }

    /** \brief Allocate memory after freeing any memory
	presently in use
    */
    void allocate(size_t nrows, size_t ncols) {
      if (mparent_t::size1>0 || mparent_t::size2>0) free();
      
      if (nrows>0 && ncols>0) {
	mparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (mparent_t::block) {
	  mparent_t::block->data=(data_t *)
	    malloc(nrows*ncols*sizeof(data_t));
	  if (mparent_t::block->data) {
	    mparent_t::block->size=nrows*ncols;
	    mparent_t::data=mparent_t::block->data;
	    mparent_t::size1=nrows;
	    mparent_t::size2=ncols;
	    mparent_t::tda=ncols;
	    mparent_t::owner=1;
	  } else {
	    std::free(mparent_t::block);
	    O2SCL_ERR("No memory for data in omatrix_tlate::allocate()",
		      gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in omatrix_tlate::allocate()",
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
      if (mparent_t::size1>0) {
	if (mparent_t::owner==1) {
	  if (mparent_t::block->size>0) {
	    std::free(mparent_t::block->data);
	  }
	  std::free(mparent_t::block);
	}
	mparent_t::size1=0;
	mparent_t::size2=0;
	mparent_t::tda=0;
      }
      return;
    }
    //@}

    /// \name Other methods
    //@{
    /// \brief Compute the transpose (even if matrix is not square)
    omatrix_tlate<data_t,mparent_t,vparent_t,block_t> transpose() {
      omatrix_tlate<data_t,mparent_t,vparent_t,block_t> 
	result(mparent_t::size2,mparent_t::size1);
      for(size_t i=0;i<mparent_t::size1;i++) {
	for(size_t j=0;j<mparent_t::size2;j++) {
	  result[j][i]=(*this)[i][j];
	}
      }
    }
    //@}
    
  };

#ifdef O2SCL_NEVER_DEFINED
}
{
#endif
  
  /** \brief Create a matrix from an array 
  */
  template<class data_t, class mparent_t, class block_t> 
    class omatrix_array_tlate : 
  public omatrix_view_tlate<data_t,mparent_t,block_t> {

  public:

    /** \brief Create a vector from \c dat with size \c n */
    omatrix_array_tlate(data_t *dat, size_t ltda,
			size_t sz1, size_t sz2)  {
      if (sz1==0 || sz2==0) {
	O2SCL_ERR("Failed to create omatrix from array, insufficient size.",
		  gsl_einval);
	mparent_t::block=0;
	mparent_t::data=0;
	mparent_t::size1=0;
	mparent_t::size2=0;
	mparent_t::tda=0;
	mparent_t::owner=0;
      }
      mparent_t::block=0;
      mparent_t::data=dat;
      mparent_t::size1=sz1;
      mparent_t::size2=sz2;
      mparent_t::tda=ltda;
      mparent_t::owner=0;
    }
  };
  
  /** \brief Create a vector from a row of a matrix 
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t> 
    class omatrix_row_tlate : 
  public ovector_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector from row \c i of matrix \c m */
    omatrix_row_tlate
      (omatrix_base_tlate<data_t,mparent_t,block_t> &m, 
       size_t i) {
      if (i<m.size1) {
	vparent_t::size=m.size2;
	vparent_t::stride=1;
	vparent_t::data=m.data+m.tda()*i;
	vparent_t::owner=0;
	vparent_t::block=0;
      } else {
	O2SCL_ERR("Invalid row in omatrix_row_tlate constructor.",
		  gsl_einval);
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::data=0;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };
    
  /** \brief Create a const vector from a row of a matrix 
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t> 
    class omatrix_const_row_tlate :
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector from row \c i of matrix \c m */
    omatrix_const_row_tlate
      (const omatrix_base_tlate<data_t,mparent_t,block_t> &m, 
       size_t i) {
      if (i<m.size1) {
	vparent_t::size=m.size2;
	vparent_t::stride=1;
	vparent_t::data=m.data+m.tda*i;
	vparent_t::owner=0;
	vparent_t::block=0;
      } else {
	O2SCL_ERR("Invalid row in omatrix_const_row_tlate constructor.",
		  gsl_einval);
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::data=0;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };
    
  /** \brief Create a vector from a column of a matrix
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t> 
    class omatrix_col_tlate :
  public ovector_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector from col \c i of matrix \c m */
    omatrix_col_tlate(omatrix_base_tlate<data_t,mparent_t,block_t> &m, 
		      size_t i) {
      if (i<m.size2) {
	vparent_t::size=m.size1;
	vparent_t::stride=m.tda();
	vparent_t::data=m.data+i;
	vparent_t::owner=0;
	vparent_t::block=0;
      } else {
	O2SCL_ERR("Invalid column in omatrix_col_tlate constructor.",
		  gsl_einval);
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::data=0;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };
    
  /** \brief Create a vector from a column of a matrix
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t> 
    class umatrix_col_tlate :
  public ovector_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector from col \c i of matrix \c m */
    umatrix_col_tlate(umatrix_base_tlate<data_t> &m, size_t i) {
      if (i<m.cols()) {
	vparent_t::size=m.rows();
	vparent_t::stride=m.cols();
	vparent_t::data=(&m.get(0,0))+i;
	vparent_t::owner=0;
	vparent_t::block=0;
      } else {
	O2SCL_ERR("Invalid column in umatrix_col_tlate constructor.",
		  gsl_einval);
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::data=0;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };
    
  /** \brief Create a const vector from a column of a matrix
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t> 
    class omatrix_const_col_tlate :
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector from col \c i of matrix \c m */
    omatrix_const_col_tlate
      (const omatrix_base_tlate<data_t,mparent_t,block_t> &m, size_t i) {
      if (i<m.size2) {
	vparent_t::size=m.size1;
	vparent_t::stride=m.tda();
	vparent_t::data=m.data+i;
	vparent_t::owner=0;
	vparent_t::block=0;
      } else {
	O2SCL_ERR("Invalid column in omatrix_const_col_tlate constructor.",
		  gsl_einval);
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::data=0;
	vparent_t::owner=0;
	vparent_t::block=0;
      }
    }
  };
  
  /** \brief Create a vector from the main diagonal
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t> 
    class omatrix_diag_tlate :
  public ovector_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector of the diagonal of matrix \c m */
    omatrix_diag_tlate(omatrix_view_tlate<data_t,mparent_t,block_t> &m) {
      
      if (m.size2<m.size1) vparent_t::size=m.size2;
      else vparent_t::size=m.size1;
      vparent_t::stride=m.tda()+1;
      vparent_t::data=m.data;
      vparent_t::owner=0;
      vparent_t::block=0;
    }
  };
  
  /** \brief Create a vector from the main diagonal
   */
  template<class data_t, class mparent_t, class vparent_t, class block_t> 
    class omatrix_const_diag_tlate :
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector of the diagonal of matrix \c m */
    omatrix_const_diag_tlate
      (const omatrix_view_tlate<data_t,mparent_t,block_t> &m) {
      
      if (m.size2<m.size1) vparent_t::size=m.size2;
      else vparent_t::size=m.size1;
      vparent_t::stride=m.tda()+1;
      vparent_t::data=m.data;
      vparent_t::owner=0;
      vparent_t::block=0;
    }
  };
  
    
  /// omatrix typedef
  typedef omatrix_tlate<double,gsl_matrix,gsl_vector_norm,gsl_block> omatrix;
  /// omatrix_view typedef
  typedef omatrix_view_tlate<double,gsl_matrix,gsl_block> omatrix_view;
  /// omatrix_const_view typedef
  typedef omatrix_const_view_tlate<double,gsl_matrix,gsl_block> 
    omatrix_const_view;
  /// omatrix_base typedef
  typedef omatrix_base_tlate<double,gsl_matrix,gsl_block> omatrix_base;
  /// omatrix_row typedef
  typedef omatrix_row_tlate<double,gsl_matrix,gsl_vector_norm,gsl_block> 
    omatrix_row;
  /// omatrix_col typedef
  typedef omatrix_col_tlate<double,gsl_matrix,gsl_vector_norm,gsl_block> 
    omatrix_col;
  /// umatrix_col typedef
  typedef umatrix_col_tlate<double,gsl_matrix,gsl_vector_norm,gsl_block> 
    umatrix_col;
  /// omatrix_const_row typedef
  typedef omatrix_const_row_tlate<double,gsl_matrix,gsl_vector_norm,gsl_block> 
    omatrix_const_row;
  /// omatrix_const_col typedef
  typedef omatrix_const_col_tlate<double,gsl_matrix,gsl_vector_norm,gsl_block> 
    omatrix_const_col;
  /// omatrix_diag typedef
  typedef omatrix_diag_tlate<double,gsl_matrix,gsl_vector_norm,gsl_block> 
    omatrix_diag;
  /// omatrix_const_diag typedef
  typedef omatrix_const_diag_tlate<double,gsl_matrix,gsl_vector_norm,
    gsl_block> omatrix_const_diag;
  /// omatrix_array typedef
  typedef omatrix_array_tlate<double,gsl_matrix,gsl_block> 
    omatrix_array;

  /// omatrix_int typedef
  typedef omatrix_tlate<int,gsl_matrix_int,gsl_vector_int,gsl_block_int> 
    omatrix_int;
  /// omatrix_int_view typedef
  typedef omatrix_view_tlate<int,gsl_matrix_int,gsl_block_int> 
    omatrix_int_view;
  /// omatrix_int_base typedef
  typedef omatrix_base_tlate<int,gsl_matrix_int,gsl_block_int> 
    omatrix_int_base;
  /// omatrix_int_const_view typedef
  typedef omatrix_const_view_tlate<int,gsl_matrix_int,gsl_block_int> 
    omatrix_int_const_view;
  /// omatrix_int_row typedef
  typedef omatrix_row_tlate<int,gsl_matrix_int,gsl_vector_int,
    gsl_block_int> omatrix_int_row;
  /// omatrix_int_col typedef
  typedef omatrix_col_tlate<int,gsl_matrix_int,gsl_vector_int,
    gsl_block_int> omatrix_int_col;
  /// umatrix_int_col typedef
  typedef umatrix_col_tlate<int,gsl_matrix_int,gsl_vector_int,
    gsl_block_int> umatrix_int_col;
  /// omatrix_int_const_row typedef
  typedef omatrix_const_row_tlate<int,gsl_matrix_int,gsl_vector_int,
    gsl_block_int> omatrix_int_const_row;
  /// omatrix_int_const_col typedef
  typedef omatrix_const_col_tlate<int,gsl_matrix_int,gsl_vector_int,
    gsl_block_int> omatrix_int_const_col;
  /// omatrix_int_diag typedef
  typedef omatrix_diag_tlate<int,gsl_matrix_int,gsl_vector_int,
    gsl_block_int> omatrix_int_diag;
  /// omatrix_int_const_diag typedef
  typedef omatrix_const_diag_tlate<int,gsl_matrix_int,gsl_vector_int,
    gsl_block_int> omatrix_int_const_diag;
  /// omatrix_int_array typedef
  typedef omatrix_array_tlate<int,gsl_matrix_int,gsl_block_int> 
    omatrix_int_array;

  /** \brief A operator for output of omatrix objects

      This outputs all of the matrix elements. Each row is output with
      an endline character at the end of each row. Positive values are
      preceeded by an extra space. A 2x2 example:
      \verbatim
      -3.751935e-05 -6.785864e-04
      -6.785864e-04  1.631984e-02
      \endverbatim

      The function \c gsl_ieee_double_to_rep() is used to determine
      the sign of a number, so that "-0.0" as distinct from "+0.0" 
      is handled correctly.

      \comment
      I had originally thought about removing this function, since
      it's superceded by matrix_out, but it's simpler to be able
      to just use the << operator so I keep it for now.
      \endcomment

      \future This assumes that scientific mode is on and showpos 
      is off. It'd be nice to fix this.

  */
  template<class data_t, class parent_t, class block_t>
    std::ostream &operator<<
    (std::ostream &os, 
     const omatrix_const_view_tlate<data_t,parent_t,block_t> &v) {
    size_t i;
    gsl_ieee_double_rep r;
    for(i=0;i<v.rows()-1;i++) {
      for(size_t j=0;j<v.cols();j++) {
	gsl_ieee_double_to_rep(&(v[i][j]), &r);
	if (r.sign==1) os << v[i][j] << ' ';
	else os << ' ' << v[i][j] << ' ';
      }
      os << '\n';
    }
    i=v.rows()-1;
    if (i>0) {
      for(size_t j=0;j<v.cols();j++) {
	gsl_ieee_double_to_rep(&(v[i][j]), &r);
	if (r.sign==1) os << v[i][j] << ' ';
	else os << ' ' << v[i][j] << ' ';
      }
    }
    return os;
  }
  
  /** \brief A simple class to provide an \c allocate() function
      for \ref omatrix
      
  */
  class omatrix_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(omatrix &o, size_t i, size_t j) { o.allocate(i,j); }
    /// Free memory
    void free(omatrix &o, size_t i) { o.free(); }
  };

#ifdef DOXGYENP2

  /** \brief A matrix where the memory allocation is performed in 
      the constructor
  */
  template<size_t N, size_t M> class ofmatrix : public omatrix_tlate {

  public:

    /// Create a matrix
    ofmatrix();
  };

#else 

  template<size_t N, size_t M> class ofmatrix : 
  public omatrix_tlate<double,gsl_matrix,gsl_vector,gsl_block>
    {
    public:
    ofmatrix() : omatrix_tlate<double,gsl_matrix,gsl_vector,
	gsl_block>(N,M) {
      }
    };
  
#endif

  /** \brief A version of \ref omatrix with better error checking

      This demonstrates how operator[] could return an ovector_array
      object and thus provide more bounds-checking. This would demand
      including a new parameter in omatrix_view_tlate which contains
      the vector type.
  */
  class xmatrix : public omatrix {
  public:
  xmatrix(size_t r=0, size_t c=0) : omatrix(r,c) {
    }
    
    /** \brief Array-like indexing 
    */
    ovector_array operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=gsl_matrix::size1) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in omatrix_view_tlate::operator[]. Size: "+
		   szttos(gsl_matrix::size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
      }
#endif
      return ovector_array(gsl_matrix::size2,gsl_matrix::data+
			   i*gsl_matrix::tda);
    }

    /** \brief Array-like indexing 
    */
    const ovector_const_array operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=gsl_matrix::size1) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in omatrix_view_tlate::operator[]. Size: "+
		   szttos(gsl_matrix::size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
      }
#endif
      return ovector_const_array(gsl_matrix::size2,gsl_matrix::data+
				 i*gsl_matrix::tda);
    }
  };

#ifndef DOXYGENP
}
#endif

#endif

