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
#ifndef O2SCL_UMATRIX_TLATE_H
#define O2SCL_UMATRIX_TLATE_H

/** \file umatrix_tlate.h
    \brief Unit-stride matrix classes
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_ieee_utils.h>

#include <o2scl/err_hnd.h>
#include <o2scl/uvector_tlate.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief A matrix view of double-precision numbers
   */
  template<class data_t> class umatrix_const_view_tlate {
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /// The data
    data_t *data;
    /// The number of rows
    size_t size1;
    /// The number of columns
    size_t size2;
    /// Zero if memory is owned elsewhere, 1 otherwise
    int owner;
    
#endif
    
  public:

    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same matrix
    umatrix_const_view_tlate(const umatrix_const_view_tlate &v) {
      data=v.data;
      size1=v.size1;
      size2=v.size2;
      owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    umatrix_const_view_tlate& operator=(const umatrix_const_view_tlate &v) {
      data=v.data;
      size1=v.size1;
      size2=v.size2;
      owner=0;      

      return *this;
    }
    //@}

    ~umatrix_const_view_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
     */
    const data_t *operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=size1) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in umatrix_const_view_tlate::operator[] const. Size: "+
		   szttos(size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return data;
      }
#endif
      return data+i*size2;
    }
    
    /** \brief Array-like indexing 
     */
    const data_t &operator()(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=size1 || j>=size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_const_view_tlate::operator() const. Sizes: ("+
		   szttos(size1)+","+szttos(size2)+
		   ") (index should be less than size).").c_str(),gsl_eindex);
	return *data;
      }
#endif
      return *(data+i*size2+j);
    }
    
    /** \brief Get (with optional range-checking) */
    data_t &get(size_t i, size_t j) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=size1 || j>=size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_const_view_tlate::get(). Sizes: ("+
		   szttos(size1)+","+szttos(size2)+
		   ") (index should be less than size).").c_str(),gsl_eindex);
	return *data;
      }
#endif
      return *(data+i*size2+j);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const data_t *get_const_ptr(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=size1 || j>=size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_const_view_tlate::get_const_ptr(). Sizes: ("+
		   szttos(size1)+","+szttos(size2)+
		   ") (index should be less than size).").c_str(),gsl_eindex);
	return (const data_t *)data;
      }
#endif
      return (const data_t *)(data+i*size2+j);
    }
    
    /** \brief Method to return number of rows
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t rows() const {
      return size1;
    }

    /** \brief Method to return number of columns
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t cols() const {
      return size2;
    }
    //@}

    /// \name Other methods
    //@{
    /// Return true if this object owns the data it refers to
    bool is_owner() const {
      if (owner==1) return true;
      return false;
    }
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor provided for use by 
	umatrix_tlate(const umatrix_tlate &v)
    */
    umatrix_const_view_tlate() {};

#endif
    
  };

  /** \brief A matrix view of double-precision numbers
   */
  template<class data_t> class umatrix_base_tlate :
  public umatrix_const_view_tlate<data_t> {
    
  public:

    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same matrix
    umatrix_base_tlate(umatrix_base_tlate &v) {
      this->data=v.data;
      this->size1=v.size1;
      this->size2=v.size2;
      this->owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    umatrix_base_tlate& operator=(umatrix_base_tlate &v) {
      this->data=v.data;
      this->size1=v.size1;
      this->size2=v.size2;
      this->owner=0;      

      return *this;
    }
    //@}

    ~umatrix_base_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
     */
    data_t *operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1) {
	O2SCL_ERR((((std::string)"Row index ")+szttos(i)+" out of bounds"
		   +" in umatrix_base_tlate::operator[]. Size: "+
		   szttos(this->size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data;
      }
#endif
      return this->data+i*this->size2;
    }
    
    /** \brief Array-like indexing 
     */
    const data_t *operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1) {
	O2SCL_ERR((((std::string)"Row index ")+szttos(i)+" out of bounds"
		   +" in umatrix_base_tlate::operator[] const. Size: "+
		   szttos(this->size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data;
      }
#endif
      return this->data+i*this->size2;
    }
    
    /** \brief Array-like indexing 
     */
    data_t &operator()(size_t i, size_t j) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1 || j>=this->size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_base_tlate::operator(). Sizes: ("+
		   szttos(this->size1)+","+szttos(this->size2)+
		   ") (index should be less than size).").c_str(),gsl_eindex);
	return *this->data;
      }
#endif
      return *(this->data+i*this->size2+j);
    }
    
    /** \brief Array-like indexing 
     */
    const data_t &operator()(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1 || j>=this->size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_base_tlate::operator() const. Sizes: ("+
		   szttos(this->size1)+","+szttos(this->size2)+
		   ") (index should be less than size).").c_str(),gsl_eindex);
	return *this->data;
      }
#endif
      return *(this->data+i*this->size2+j);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i, size_t j) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1 || j>=this->size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_base_tlate::get_ptr(). Sizes: ("+
		   szttos(this->size1)+","+szttos(this->size2)+
		   ") (index should be less than size).").c_str(),gsl_eindex);
	return this->data;
      }
#endif
      return this->data+i*this->size2+j;
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, size_t j, data_t val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1 || j>=this->size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_base_tlate::set(). Sizes: ("+
		   szttos(this->size1)+","+szttos(this->size2)+
		   ") (index should be less than size).").c_str(),
		  gsl_eindex);
      }
#endif
      *(this->data+i*this->size2+j)=val;
      return;
    }

    /** \brief Set all of the value to be the value \c val */
    void set_all(data_t val) {
      for(size_t i=0;i<this->size1;i++) {
	for(size_t j=0;j<this->size2;j++) {
	  *(this->data+i*this->size2+j)=val;
	}
      }
      return;
    }
    //@}

    /// \name Arithmetic 
    //@{
    /** \brief operator+= */
    umatrix_base_tlate<data_t> &operator+=
      (const umatrix_base_tlate<data_t> &x) {
      size_t lsize=x.size1;
      if (lsize>this->size1) lsize=this->size1;
      size_t lsize2=x.size2;
      if (lsize2>this->size2) lsize2=this->size2;
      for(size_t i=0;i<lsize;i++) {
	for(size_t j=0;j<lsize2;j++) {
	  (*this)[i][j]+=x[i][j];
	}
      }
      
      return *this;
    }
    
    /** \brief operator-= */
    umatrix_base_tlate<data_t> &operator-=
      (const umatrix_base_tlate<data_t> &x) {
      size_t lsize=x.size1;
      if (lsize>this->size1) lsize=this->size1;
      size_t lsize2=x.size2;
      if (lsize2>this->size2) lsize2=this->size2;
      for(size_t i=0;i<lsize;i++) {
	for(size_t j=0;j<lsize2;j++) {
	  (*this)[i][j]+=x[i][j];
	}
      }

      return *this;
    }
    
    /** \brief operator+= */
    umatrix_base_tlate<data_t> &operator+=(const data_t &y) {
      for(size_t i=0;i<this->size1;i++) {
	for(size_t j=0;j<this->size2;j++) {
	  (*this)[i][j]+=y;
	}
      }
      
      return *this;
    }

    /** \brief operator-= */
    umatrix_base_tlate<data_t> &operator-=(const data_t &y) {
      for(size_t i=0;i<this->size1;i++) {
	for(size_t j=0;j<this->size2;j++) {
	  (*this)[i][j]-=y;
	}
      }
      
      return *this;
    }

    /** \brief operator*= */
    umatrix_base_tlate<data_t> &operator*=(const data_t &y) {
      for(size_t i=0;i<this->size1;i++) {
	for(size_t j=0;j<this->size2;j++) {
	  (*this)[i][j]*=y;
	}
      }
      
      return *this;
    }
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor
     */
    umatrix_base_tlate() {};

#endif
    
  };

  /** \brief A matrix view of double-precision numbers
   */
  template<class data_t> class umatrix_view_tlate :
  public umatrix_base_tlate<data_t> {
    
  public:

    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same matrix
    umatrix_view_tlate(const umatrix_view_tlate &v) {
      this->data=v.data;
      this->size1=v.size1;
      this->size2=v.size2;
      this->owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    umatrix_view_tlate& operator=(const umatrix_view_tlate &v) {
      this->data=v.data;
      this->size1=v.size1;
      this->size2=v.size2;
      this->owner=0;      

      return *this;
    }

    /// Shallow copy constructor - create a new view of the same matrix
    umatrix_view_tlate(umatrix_base_tlate<data_t> &v) {
      this->data=v.data;
      this->size1=v.size1;
      this->size2=v.size2;
      this->owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same matrix
    umatrix_view_tlate& operator=(umatrix_base_tlate<data_t> &v) {
      this->data=v.data;
      this->size1=v.size1;
      this->size2=v.size2;
      this->owner=0;      

      return *this;
    }
    //@}

    ~umatrix_view_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
     */
    data_t *operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1) {
	O2SCL_ERR((((std::string)"Row index ")+szttos(i)+" out of bounds"
		   +" in umatrix_view_tlate::operator[]. Size: "+
		   szttos(this->size1)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data;
      }
#endif
      return this->data+i*this->size2;
    }
    
    /** \brief Array-like indexing 
     */
    data_t &operator()(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1 || j>=this->size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_view_tlate::operator(). Sizes: ("+
		   szttos(this->size1)+","+szttos(this->size2)+
		   ") (index should be less than size).").c_str(),gsl_eindex);
	return *this->data;
      }
#endif
      return *(this->data+i*this->size2+j);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i, size_t j) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1 || j>=this->size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_view_tlate::get_ptr(). Sizes: ("+
		   szttos(this->size1)+","+szttos(this->size2)+
		   ") (index should be less than size).").c_str(),gsl_eindex);
	return this->data;
      }
#endif
      return this->data+i*this->size2+j;
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, size_t j, data_t val) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->size1 || j>=this->size2) {
	O2SCL_ERR((((std::string)"Indices (")+szttos(i)+","+szttos(j)+
		   ") out of bounds"
		   +" in umatrix_view_tlate::set(). Sizes: ("+
		   szttos(this->size1)+","+szttos(this->size2)+
		   ") (index should be less than size).").c_str(),
		  gsl_eindex);
      }
#endif
      *(this->data+i*this->size2+j)=val;
      return;
    }

    /** \brief Set all of the value to be the value \c val */
    void set_all(double val) const {
      for(size_t i=0;i<this->size1;i++) {
	for(size_t j=0;j<this->size2;j++) {
	  *(this->data+i*this->size2+j)=val;
	}
      }
      return;
    }
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor
     */
    umatrix_view_tlate() {};

#endif
    
  };

  /** \brief A matrix of double-precision numbers
   */
  template<class data_t> class umatrix_tlate : 
  public umatrix_base_tlate<data_t> {
  public:
    
    /// \name Standard constructor
    //@{
    /** \brief Create an umatrix of size \c n with owner as 'true' 
     */
    umatrix_tlate(size_t r=0, size_t c=0) {

      this->data=0;
      this->size1=0;
      this->size2=0;

      // This must be set to 1 even if n=0 so that future
      // calls to operator= work properly
      this->owner=1;
      
      if (r>0 && c>0) {
	this->data=new data_t[r*c];
	this->size1=r;
	this->size2=c;
      }
    }
    //@}
    
    /// \name Copy constructors
    //@{
    /// Deep copy constructor, allocate new space and make a copy
  umatrix_tlate(const umatrix_tlate &v) : umatrix_base_tlate<data_t>() {
      size_t n=v.size1;
      size_t n2=v.size2;
      if (n>0 && n2>0) {
	this->data=new data_t[n*n2];
	this->size1=n;
	this->size2=n2;
	this->owner=1;
	for(size_t i=0;i<n;i++) {
	  for(size_t j=0;j<n2;j++) {
	    *(this->data+i*this->size2+j)=v[i][j];
	  }
	}
      } else {
	this->size1=0;
	this->size2=0;
      }
    }
    
    /// Deep copy constructor, allocate new space and make a copy
    umatrix_tlate
      (const umatrix_view_tlate<data_t> &v) : 
    umatrix_base_tlate<data_t>() {
      size_t r=v.rows();
      size_t c=v.cols();
      if (r>0 && c>0) {
	this->data=new data_t[r*c];
	this->size1=v.size1;
	this->size2=v.size2;
	this->owner=1;
	for(size_t i=0;i<r;i++) {
	  for(size_t j=0;i<c;j++) {
	    *(this->data+i*this->size2+j)=v[i][j];
	  }
	}
      } else {
	this->size1=0;
	this->size2=0;
      }
    }
    
    /** \brief Deep copy constructor, if owner is true, allocate space and
	make a new copy, otherwise, just copy into the view
    */
    umatrix_tlate& operator=(const umatrix_tlate &v) {
      
      // Check for self-assignment
      if (this==&v) return *this;

      size_t sze=v.size1;
      size_t sze2=v.size2;
      
      if (sze==0 && sze2==0) return *this;

      if (this->owner) {
	allocate(sze,sze2);
      } else {
	if (this->size1!=sze || this->size2!=sze2) {
	  O2SCL_ERR("Sizes don't match in umatrix_tlate::operator=()",
		    gsl_ebadlen);
	  return *this;
	}
      }
      for(size_t i=0;i<sze;i++) {
	for(size_t j=0;j<sze2;j++) {
	  *(this->data+i*this->size2+j)=v[i][j];
	}
      }

      return *this;
    }

    /** \brief Deep copy constructor, if owner is true, allocate space and
	make a new copy, otherwise, just copy into the view
    */
    umatrix_tlate& operator=
      (const umatrix_view_tlate<data_t> &v) {

      // Check for self-assignment
      if (this==&v) return *this;

      size_t sze=v.rows();
      size_t sze2=v.cols();

      if (sze==0 && sze2==0) return *this;

      if (this->owner) {
	allocate(sze,sze2);
      } else {
	if (this->size1!=sze || this->size2!=sze2) {
	  O2SCL_ERR("Sizes don't match in umatrix_tlate::operator=()",
		    gsl_ebadlen);
	  return *this;
	}
      }
      for(size_t i=0;i<sze;i++) {
	for(size_t j=0;j<sze2;j++) {
	  *(this->data+i*this->size2+j)=v[i][j];
	}
      }

      return *this;
    }
    
#ifdef O2SCL_NEVER_DEFINED

    /*
      AWS, 11/30/12: It's not clear what this was intended
      for, and the 'n2=uva[0]' line doesn't look right 
      anyway. We comment it out for now. 
    */
    
    /** \brief Deep copy from an array of uvectors
     */
    umatrix_tlate(size_t n, uvector_view_tlate<data_t> uva[]) {
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

#endif

    /** \brief Deep copy from a C-style 2-d array
     */
    umatrix_tlate(size_t n, size_t n2, data_t **csa) {
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
    
    ~umatrix_tlate() {
      if (this->size1>0) {
	if (this->owner==1) {
	  delete[] this->data;
	  this->size1=0;
	  this->size2=0;
	}
      }
    }

    /// \name Memory allocation
    //@{
    /** \brief Reallocate memory for new size
	    
	If \c nrows and \c ncols are both zero, then all memory is
	deallocated. Otherwise, new memory is allocated and the
	old contents are copied into the new memory (up to the
	lesser of the old and new size).
    */
    void resize(size_t nrows, size_t ncols) {
      // New size is zero
      if (nrows==0 && ncols==0) {
	free();
	return;
      }
      // Old size was zero
      if (this->size1==0 && this->size2==0) {
	this->data=new data_t[nrows*ncols];
	this->size1=nrows;
	this->size2=ncols;
	this->owner=1;
	return;
      }

      // Neither size was zero
      data_t *newp=new data_t[nrows*ncols];
      size_t r=nrows, c=ncols;
      if (this->size1<nrows) r=this->size1;
      if (this->size2<ncols) c=this->size2;
      vector_copy(r*c,this->data,newp);
      if (this->owner==1) {
	delete[] this->data;
      }
      this->data=newp;
      this->size1=nrows;
      this->size2=ncols;
      this->owner=1;
      return;
    }

    /** \brief Allocate memory after freeing any memory presently in use
     */
    void allocate(size_t nrows, size_t ncols) {
    
      if (this->size1>0 || this->size2>0) free();
      
      if (nrows>0 && ncols>0) {
	this->data=new data_t[nrows*ncols];
	this->size1=nrows;
	this->size2=ncols;
	this->owner=1;
      } else {
	O2SCL_ERR("Zero size in umatrix::allocate()",gsl_einval);
      }
      return;
    }

    /** \brief Free the memory 
    
	This function will safely do nothing if used without first
	allocating memory or if called multiple times in succession.
    */
    void free() {
      if (this->size1>0) {
	if (this->owner==1) {
	  delete[] this->data;
	}
	this->size1=0;
	this->size2=0;
      }
      return;
    }
    //@}

    /// \name Other methods
    //@{
    /** \brief Compute the transpose (even if matrix is not square)
     */
    umatrix_tlate<data_t> transpose() {
      umatrix_tlate<data_t> result(this->size2,this->size1);
      for(size_t i=0;i<this->size1;i++) {
	for(size_t j=0;j<this->size2;j++) {
	  result[j][i]=(*this)[i][j];
	}
      }
    }
    //@}
    
  };

  /** \brief Create a vector from a row of a matrix 
   */
  template<class data_t> class umatrix_row_tlate : 
  public uvector_view_tlate<data_t> {
  public:
    /** \brief Create a vector from row \c i of matrix \c m */
    umatrix_row_tlate(umatrix_base_tlate<data_t> &m, 
		      size_t i) {
      this->sz=0;
      this->data=0;
      this->owner=0;
      if (i<m.rows()) {
	this->sz=m.cols();
	this->data=m[0]+m.cols()*i;
      }
    }
  };
    
  /** \brief Create a const vector from a row of a matrix 
   */
  template<class data_t> class umatrix_const_row_tlate :
  public uvector_const_view_tlate<data_t> {
  public:
    /** \brief Create a vector from row \c i of matrix \c m */
    umatrix_const_row_tlate
      (const umatrix_base_tlate<data_t> &m, 
       size_t i) {
      if (i<m.size1) {
	this->sz=m.cols();
	this->data=m[0]+m.cols()*i;
	this->owner=0;
      }
    }
  };
    
  /// umatrix typedef
  typedef umatrix_tlate<double> umatrix;
  /// umatrix_view typedef
  typedef umatrix_view_tlate<double> umatrix_view;
  /// umatrix_base typedef
  typedef umatrix_base_tlate<double> umatrix_base;
  /// umatrix_const_view typedef
  typedef umatrix_const_view_tlate<double> umatrix_const_view;
  /// umatrix_row typedef
  typedef umatrix_row_tlate<double> umatrix_row;
  /// umatrix_const_row typedef
  typedef umatrix_const_row_tlate<double> umatrix_const_row;

  /// umatrix_int typedef
  typedef umatrix_tlate<int> umatrix_int;
  /// umatrix_int_view typedef
  typedef umatrix_view_tlate<int> umatrix_int_view;
  /// umatrix_int_const_view typedef
  typedef umatrix_const_view_tlate<int> umatrix_int_const_view;
  /// umatrix_int_base typedef
  typedef umatrix_base_tlate<int> umatrix_int_base;
  /// umatrix_int_row typedef
  typedef umatrix_row_tlate<int> umatrix_int_row;
  /// umatrix_int_const_row typedef
  typedef umatrix_const_row_tlate<int> umatrix_int_const_row;

  /// umatrix_size_t typedef
  typedef umatrix_tlate<size_t> umatrix_size_t;
  /// umatrix_size_t_view typedef
  typedef umatrix_view_tlate<size_t> umatrix_size_t_view;
  /// umatrix_size_t_const_view typedef
  typedef umatrix_const_view_tlate<size_t> umatrix_size_t_const_view;
  /// umatrix_size_t_base typedef
  typedef umatrix_base_tlate<size_t> umatrix_size_t_base;
  /// umatrix_size_t_row typedef
  typedef umatrix_row_tlate<size_t> umatrix_size_t_row;
  /// umatrix_size_t_const_row typedef
  typedef umatrix_const_row_tlate<size_t> umatrix_size_t_const_row;

  /** \brief A operator for simple matrix output
      
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

      \future This assumes that scientific mode is on and showpos 
      is off. It'd be nice to fix this.

  */
  template<class data_t> std::ostream &operator<<
    (std::ostream &os, const umatrix_const_view_tlate<data_t> &v) {
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
      for \ref umatrix
  */
  class umatrix_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(umatrix &o, int i, int j) { o.allocate(i,j); }
    /// Free memory
    void free(umatrix &o) { o.free(); }
  };

  /** \brief A matrix where the memory allocation is performed in 
      the constructor
  */
#ifdef DOXYGENP
  template<size_t N, size_t M> class ufmatrix : 
  public umatrix_tlate<data_t>
#else  
    template<size_t N, size_t M> class ufmatrix :   
  public umatrix_tlate<double>
#endif
  {
  public:
  ufmatrix() : umatrix_tlate<double>(N,M) {
    }
  };
  
  
#ifndef DOXYGENP
}
#endif

#endif



