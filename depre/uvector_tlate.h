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
#ifndef O2SCL_UVECTOR_TLATE_H
#define O2SCL_UVECTOR_TLATE_H

/** \file uvector_tlate.h
    \brief File for definitions of unit-stride vectors
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

// For gsl_finite
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>
#include <o2scl/array.h>
#include <o2scl/vector.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  // Forward declaration
  template<class data_t> class uvector_base_tlate;
  // Forward declaration
  template<class data_t> class uvector_view_tlate;
  // Forward declaration
  template<class data_t> class uvector_tlate;
  // Forward declarations so that we can friend these classes in
  // uvector_const_view_tlate below
  template<class data_t> class uvector_subvector_tlate;
  template<class data_t> class uvector_const_subvector_tlate;

  /** \brief A const vector view with unit stride

      See the \ref vecmat_section section of the User's guide
      for general information about vector classes in \o2 .
      
      See the explicit instantations of this class, \ref
      uvector_const_view, \ref uvector_int_const_view , and \ref
      uvector_size_t_const_view .

      \future Could allow user-defined specification of 
      restrict keyword
  */
  template<class data_t> class uvector_const_view_tlate {
    
#ifndef DOXYGEN_INTERNAL

  protected:

    // Make sure everyone can access the base data
    friend class uvector_base_tlate<data_t>;
    friend class uvector_view_tlate<data_t>;
    friend class uvector_tlate<data_t>;
    friend class uvector_subvector_tlate<data_t>;
    friend class uvector_const_subvector_tlate<data_t>;

    /// The data
    data_t *data;
    /// The vector sz
    size_t sz;
    /// Zero if memory is owned elsewhere, 1 otherwise
    int owner;
    
#endif

  public:

    /// \name Copy constructors
    //@{
    /// Copy constructor - create a new view of the same vector
    uvector_const_view_tlate(const uvector_const_view_tlate &v) {
      data=v.data;
      sz=v.sz;
      owner=0;      
    }
    
    /// Copy constructor - create a new view of the same vector
    uvector_const_view_tlate& operator=(const uvector_const_view_tlate &v) {
      data=v.data;
      sz=v.sz;
      owner=0;      

      return *this;
    }
    //@}

    ~uvector_const_view_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
     */
    const data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::operator[] const. Size: "+
		   szttos(sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return data[0];
      }
#endif
      return data[i];
    }
    
    /** \brief Array-like indexing 
     */
    const data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::operator() const. Size: "+
		   szttos(sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return data[0];
      }
#endif
      return data[i];
    }
    
    /** \brief Get (with optional range-checking) */
    data_t get(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::get(). Size: "+
		   szttos(sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return data[0];
      }
#endif
      return data[i];
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const data_t *get_const_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR2("Index out of range in ",
		   "uvector_const_view::get_const_ptr().",1);
	return 0;
      }
#endif
      return (const data_t *)(data+i);
    }
    
    /** \brief Method to return vector size 
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t size() const {
      return sz;
    }
    //@}

    /// \name Other methods
    //@{

    /// Return true if this object owns the data it refers to
    bool is_owner() const {
      if (owner==1) return true;
      return false;
    }

    /** \brief Exhaustively look through the array for a
	particular value 

	This can only fail if \em all of the entries in the array are
	not finite, in which case it calls O2SCL_ERR() and returns
	0. 

	If more than one entry is the same distance from <tt>x0</tt>,
	this function returns the entry with smallest index.
    */
    size_t lookup(const data_t x0) const {

      const uvector_const_view_tlate<data_t> *a=this;
      size_t row=0, i=0, nvar=size();
      while(!gsl_finite((*a)[i]) && i<nvar-1) i++;
      if (i==nvar-1) {
	O2SCL_ERR("Array not finite in intp_base::lookup()",1);
	return 0;
      }
      data_t best=(*a)[i], bdiff=fabs((*a)[i]-x0);
      for(;i<nvar;i++) {
	if (gsl_finite((*a)[i]) && fabs((*a)[i]-x0)<bdiff) {
	  row=i;
	  best=(*a)[i];
	  bdiff=fabs((*a)[i]-x0);
	}
      }
      return row;
    }

    /** \brief Find the maximum element */
    data_t max() const {
      data_t maxval;
      if (sz>0) {
	maxval=data[0];
	for(size_t i=1;i<sz;i++) {
	  if (data[i]>maxval) {
	    maxval=data[i];
	  }
	}
      } else {
	return 0;
      }
      return maxval;
    }

    /** \brief Find the minimum element */
    data_t min() const {
      data_t minval;
      if (sz>0) {
	minval=data[0];
	for(size_t i=1;i<sz;i++) {
	  if (data[i]<minval) {
	    minval=data[i];
	  }
	}
      } else {
	return 0.0;
      }
      return minval;
    }
    //@}

    /** \brief Norm 

	\todo Fix this so that norm() is computed as 
	in ovector and so that integer norms are performed
	separately.
    */
    data_t norm() const {
      data_t result=0;
      for(size_t i=0;i<sz;i++) {
	result+=(*this)[i]*(*this)[i];
      }
      return std::sqrt(result);
    }
        
#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor provided for use by 
	uvector_tlate(const uvector_tlate &v)
    */
    uvector_const_view_tlate() {};

#endif
    
  };
  
  /** \brief A base class for uvector and uvector_view
      
      This class provides a base class for uvector and uvector_view,
      mostly useful for creating function arguments which accept
      either uvector or uvector_view. 

      See the \ref vecmat_section section of the User's guide
      for general information about vector classes in \o2 .
      
      See the explicit instantations of this class, \ref
      uvector_base, \ref uvector_int_base , and \ref
      uvector_size_t_base .

  */
  template<class data_t> class uvector_base_tlate : 
  public uvector_const_view_tlate<data_t> {
    
  public:
    
    /// \name Copy constructors
    //@{
    /// Copy constructor - create a new view of the same vector
  uvector_base_tlate(uvector_base_tlate &v) : 
    uvector_const_view_tlate<data_t>() {
      this->data=v.data;
      this->sz=v.sz;
      this->owner=0;      
    }
    
    /// Copy constructor - create a new view of the same vector
    uvector_base_tlate& operator=(uvector_base_tlate &v) {
      this->data=v.data;
      this->sz=v.sz;
      this->owner=0;      

      return *this;
    }
    //@}

    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
     */
    data_t &operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::operator[]. Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data[0];
      }
#endif
      return this->data[i];
    }
    
    /** \brief Array-like indexing 
     */
    const data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::operator[] const. Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data[0];
      }
#endif
      return this->data[i];
    }
    
    /** \brief Array-like indexing 
     */
    data_t &operator()(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::operator(). Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data[0];
      }
#endif
      return this->data[i];
    }
    
    /** \brief Array-like indexing 
     */
    const data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::operator() const. Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data[0];
      }
#endif
      return this->data[i];
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::get_ptr(). Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data;
      }
#endif
      return this->data+i;
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, data_t val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_const_view_tlate::set(). Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),
		  gsl_eindex);
      }
#endif
      this->data[i]=val;
      return;
    }

    /** \brief Set all of the value to be the value \c val */
    void set_all(data_t val) {
      for(size_t i=0;i<this->sz;i++) {
	this->data[i]=val;
      }
      return;
    }

    //@}

    /// \name Arithmetic
    //@{
    /** \brief operator+= */
    uvector_base_tlate<data_t> &operator+=
      (const uvector_base_tlate<data_t> &x) {
      size_t lsz=x.size();
      if (lsz>this->sz) lsz=this->sz;
      for(size_t i=0;i<lsz;i++) (*this)[i]+=x[i];
      
      return *this;
    }
    
    /** \brief operator-= */
    uvector_base_tlate<data_t> &operator-=
      (const uvector_base_tlate<data_t> &x) {
      size_t lsz=x.size();
      if (lsz>this->sz) lsz=this->sz;
      for(size_t i=0;i<lsz;i++) (*this)[i]-=x[i];
      
      return *this;
    }
    
    /** \brief operator+= */
    uvector_base_tlate<data_t> &operator+=(const data_t &y) {
      for(size_t i=0;i<this->sz;i++) (*this)[i]+=y;

      return *this;
    }

    /** \brief operator-= */
    uvector_base_tlate<data_t> &operator-=(const data_t &y) {
      for(size_t i=0;i<this->sz;i++) (*this)[i]-=y;

      return *this;
    }

    /** \brief operator*= */
    uvector_base_tlate<data_t> &operator*=(const data_t &y) {
      for(size_t i=0;i<this->sz;i++) (*this)[i]*=y;

      return *this;
    }
    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor for use by children

	\comment
	This is private because I want to prevent the casual end-user
	from creating uvector_base objects. End-users should really be
	focused on making normal vectors and views.
	\endcomment
    */
    uvector_base_tlate() {}

#endif

  };
    
#ifdef O2SCL_NEVER_DEFINED
}{
#endif

  /** \brief A base class for uvectors

      See the \ref vecmat_section section of the User's guide
      for general information about vector classes in \o2 .
      
      See the explicit instantations of this class, \ref
      uvector_view, \ref uvector_int_view , and \ref
      uvector_size_t_view .
   */
  template<class data_t> class uvector_view_tlate : 
  public uvector_base_tlate<data_t> {

  public:

    /// \name Copy constructors
    //@{
    /// Copy constructor - create a new view of the same vector
  uvector_view_tlate(const uvector_view_tlate &v) : 
    uvector_base_tlate<data_t>() {
      this->data=v.data;
      this->sz=v.sz;
      this->owner=0;      
    }
    
    /// Copy constructor - create a new view of the same vector
    uvector_view_tlate& operator=(const uvector_view_tlate &v) {
      this->data=v.data;
      this->sz=v.sz;
      this->owner=0;      

      return *this;
    }

    /// Copy constructor - create a new view of the same vector
  uvector_view_tlate(uvector_base_tlate<data_t> &v) :
    uvector_base_tlate<data_t>() {
      this->data=v.data;
      this->sz=v.sz;
      this->owner=0;      
    }
    
    /// Copy constructor - create a new view of the same vector
    uvector_view_tlate& operator=(uvector_base_tlate<data_t> &v) {
      this->data=v.data;
      this->sz=v.sz;
      this->owner=0;      

      return *this;
    }
    //@}

    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
     */
    data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_view_tlate::operator[]. Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data[0];
      }
#endif
      return this->data[i];
    }
    
    /** \brief Array-like indexing 
     */
    data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_view_tlate::operator(). Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data[0];
      }
#endif
      return this->data[i];
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_view_tlate::get_ptr(). Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return this->data;
      }
#endif
      return this->data+i;
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, data_t val) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=this->sz) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in uvector_view_tlate::set(). Size: "+
		   szttos(this->sz)+
		   " (index should be less than size).").c_str(),
		  gsl_eindex);
      }
#endif
      this->data[i]=val;
      return;
    }

    /** \brief Set all of the value to be the value \c val */
    void set_all(data_t val) const {
      for(size_t i=0;i<this->sz;i++) {
	this->data[i]=val;
      }
      return;
    }

    //@}

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor provided for use by 
	uvector_view_tlate(const uvector_view_tlate &v)
    */
    uvector_view_tlate() {}

#endif

  };
  
  /** \brief A vector with unit stride
      
      See the \ref vecmat_section section of the User's guide
      for general information about vector classes in \o2 .
      
      See the explicit instantations of this class, \ref
      uvector , \ref uvector_int , and \ref
      uvector_size_t .

      There are several useful methods which are defined in the parent
      class, \ref uvector_const_view_tlate . There is also an <<
      operator for this class documented "Functions" section of \ref
      uvector_tlate.h.
  */
  template<class data_t> class uvector_tlate : 
  public uvector_base_tlate<data_t> {

  public:

    // This is required so that the uvector_view_tlate constructor
    // can access uvector_tlate data.
    friend class uvector_view_tlate<data_t>;

    /// \name Standard constructor
    //@{
    /** \brief Create an uvector of size \c n with owner as 'true' */
    uvector_tlate(size_t n=0) {
      
      this->sz=0;
      this->data=0;

      // This must be set to 1 even if n=0 so that future
      // calls to operator= work properly
      this->owner=1;

      if (n>0) {
	this->data=new data_t[n];
	if (this->data) {
	  this->sz=n;
	} else {
	  O2SCL_ERR("No memory for data in uvector_tlate constructor",
		    gsl_enomem);
	}
      }
    }
    //@}
    
    /// \name Copy constructors
    //@{
    /// Deep copy constructor - allocate new space and make a copy
  uvector_tlate(const uvector_tlate &v) : 
    uvector_base_tlate<data_t>() {

      this->sz=0;
      this->data=0;
      this->owner=1;

      size_t n=v.sz;
      if (n>0) {
	this->data=new data_t[n];
	if (this->data) {
	  this->sz=n;
	  this->owner=1;
	  for(size_t i=0;i<n;i++) {
	    this->data[i]=v[i];
	  }
	} else {
	  O2SCL_ERR("No memory for data in uvector_tlate constructor",
		    gsl_enomem);
	}
      } else {
	this->sz=0;
      }
    }
    
    /// Deep copy constructor - allocate new space and make a copy
  uvector_tlate(const uvector_const_view_tlate<data_t> &v) : 
    uvector_base_tlate<data_t>() {

      this->sz=0;
      this->data=0;
      this->owner=1;

      size_t n=v.size();
      if (n>0) {
	this->data=new data_t[n];
	if (this->data) {
	  this->sz=n;
	  this->owner=1;
	  for(size_t i=0;i<n;i++) {
	    this->data[i]=v[i];
	  }
	} else {
	  O2SCL_ERR("No memory for data in uvector_tlate constructor",
		    gsl_enomem);
	}
      } else {
	this->sz=0;
      }
    }

    /** \brief Deep copy constructor - if owner is true, allocate space and
	make a new copy, otherwise, just copy into the view
    */
    uvector_tlate& operator=(const uvector_tlate &v) {

      // Check for self-assignment
      if (this==&v) return *this;

      size_t sz2=v.sz;
      this->sz=0;
      this->data=0;
      if (this->owner) {
	allocate(sz2);
      } else {
	if (this->sz!=sz2) {
	  O2SCL_ERR("Sizes don't match in uvector_tlate::operator=()",
		    gsl_ebadlen);
	  return *this;
	}
      }
      for(size_t i=0;i<sz2;i++) {
	this->data[i]=v[i];
      }
      return *this;
    }

    /** \brief Deep copy constructor - if owner is true, allocate space and
	make a new copy, otherwise, just copy into the view
    */
    uvector_tlate& operator=(const uvector_const_view_tlate<data_t> &v) {

      // Check for self-assignment
      if (this==&v) return *this;

      size_t sz2=v.size();
      this->sz=0;
      this->data=0;
      if (this->owner) {
	allocate(sz2);
      } else {
	if (this->sz!=sz2) {
	  O2SCL_ERR("Sizes don't match in uvector_tlate::operator=()",
		    gsl_ebadlen);
	  return *this;
	}
      }
      for(size_t i=0;i<sz2;i++) {
	this->data[i]=v[i];
      }
      return *this;
    }
    //@}

    ~uvector_tlate() {
      if (this->sz>0) {
	if (this->owner==1) {
	  delete[] this->data;
	}
	this->sz=0;
      }
    }

    /// \name Memory allocation
    //@{
    /** \brief Reallocate memory for size \c nsize

	If \c nsize is zero, then all memory is deallocated. Otherwise,
	new memory is allocated and the old contents are copied into 
	the new memory (up to the lesser of the old and new size).
    */
    void resize(size_t nsize) {
      // New size is zero
      if (nsize==0) {
	free();
	return;
      }
      // Old size was zero
      if (this->sz==0) {
	this->data=new data_t[nsize];
	this->sz=nsize;
	this->owner=1;
	return;
      }
      // New and old size are the same
      if (nsize==this->sz) return;
      // Neither size was zero
      data_t *newp=new data_t[nsize];
      if (nsize<this->sz) {
	vector_copy(nsize,this->data,newp);
      } else {
	vector_copy(this->sz,this->data,newp);
      }
      if (this->owner==1) {
	delete[] this->data;
      }
      this->data=newp;
      this->sz=nsize;
      this->owner=1;
      return;
    }
    
    /** \brief Allocate memory for size \c nsize after freeing any memory 
	presently in use
	
	If \c nsize is zero, this only frees the memory and allocates
	no additional memory. 
	\comment
	9/16/09 - I used to set an error here for empty allocations,
	but that was causing problems with std::vector<uvector> type
	objects, and also we don't want the user to have to make an
	extra nsize>0 check every time they call allocate(), so we
	don't set an error here.
	\endcomment
    */
    void allocate(size_t nsize) {
      if (this->sz>0) free();

      if (nsize>0) {
	this->data=new data_t[nsize];
	this->sz=nsize;
	this->owner=1;
      }
      return;
    }

    /** \brief Free the memory 
    
	This function will safely do nothing if used without first
	allocating memory or if called multiple times in succession.
    */
    void free() {
      if (this->sz>0) {
	if (this->owner==1) {
	  delete[] this->data;
	}
	this->sz=0;
      }
      return;
    }
    //@}

    /// \name Other methods
    //@{
    /// Erase an element from the array. 
    void erase(size_t ix) {

      // Only proceed if the user gave an element inside the array
      if (this->sz>ix) {

	// Decrement the size 
	this->sz--;

	// Allocate new space
	data_t *newdat=new data_t[this->sz];

	// Copy into the new space
	for(size_t i=0;i<this->sz;i++) {
	  if (i<ix) newdat[i]=this->data[i];
	  else newdat[i]=this->data[i+1];
	}

	// Free the old space and reset the pointer
	  delete[] this->data;
	this->data=newdat;

      } else {
	O2SCL_ERR("Cannot erase() beyond end of uvector_tlate.",
		  gsl_einval);
      }

      return;
    }
    //@}

    /** \brief Sort the vector and ensure all elements are unique by 
	removing duplicates

	\note This just sorts the vector and deletes duplicate 
	elements afterwards. 
    */
    void sort_unique() {
      o2scl::vector_sort<uvector_tlate<data_t>,data_t >(this->sz,*this);
      for(size_t i=1;i<this->sz;i++) {
	if ((*this)[i]==(*this)[i-1]) {
	  erase(i-1);
	  i--;
	}
      }
      return;
    }
    //@}

  };

  /** \brief Create a vector from an array 
   */
  template<class data_t> class uvector_array_tlate : 
  public uvector_view_tlate<data_t> {
  public:
    /** \brief Create a vector from \c dat with size \c n */
    uvector_array_tlate(size_t n, data_t *dat) {
      if (n>0) {
	this->data=dat;
	this->sz=n;
	this->owner=0;
      }
    }
  };
    
  /** \brief Create a vector from a subvector of another
   */
  template<class data_t> class uvector_subvector_tlate : 
  public uvector_view_tlate<data_t> {
  public:
    /** \brief Create a vector from \c orig */
    uvector_subvector_tlate(uvector_base_tlate<data_t> &orig, 
			    size_t offset, size_t n) {
      if (offset+n-1<orig.size()) {
	this->data=orig.data+offset;
	this->sz=n;
	this->owner=0;
      } else {
	this->sz=0;
	O2SCL_ERR("Subvector failed in uvector_sub_view().",1);
      }
    }
  };

  /** \brief Create a vector from an const array 
   */
  template<class data_t> class uvector_const_array_tlate :
  public uvector_const_view_tlate<data_t> {
  public:
    /** \brief Create a vector from \c dat with size \c n */
    uvector_const_array_tlate(size_t n, const data_t *dat) {
      if (n>0) {
	// We have to do an explicit cast here, but we prevent the
	// user from changing the data.
	this->data=(data_t *)dat;
	this->sz=n;
	this->owner=0;
      }
    }
    
    ~uvector_const_array_tlate() {};
    
  };

  /** \brief Create a const vector from a subvector of another vector
   */
  template<class data_t> class uvector_const_subvector_tlate :
  public uvector_const_view_tlate<data_t> {
  public:
    /** \brief Create a vector from \c orig 
     */
    uvector_const_subvector_tlate
      (const uvector_base_tlate<data_t> &orig, size_t offset, size_t n) {
      if (offset+n-1<orig.sz) {
	this->data=orig.data+offset;
	this->sz=n;
	this->owner=0;
      } else {
	this->sz=0;
	O2SCL_ERR("Subvector failed in uvector_subvector().",1);
      }
    }

  };

  /// uvector typedef
  typedef uvector_tlate<double> uvector;
  /// uvector_view typedef
  typedef uvector_view_tlate<double> uvector_view;
  /// uvector_base typedef
  typedef uvector_base_tlate<double> uvector_base;
  /// uvector_const_view typedef
  typedef uvector_const_view_tlate<double> uvector_const_view;
  /// uvector_array typedef
  typedef uvector_array_tlate<double> uvector_array;
  /// uvector_subvector typedef
  typedef uvector_subvector_tlate<double> uvector_subvector;
  /// uvector_const_array typedef
  typedef uvector_const_array_tlate<double> uvector_const_array;
  /// uvector_const_subvector typedef
  typedef uvector_const_subvector_tlate<double> uvector_const_subvector;
  
  /// uvector_int typedef
  typedef uvector_tlate<int> uvector_int;
  /// uvector_int_view typedef
  typedef uvector_view_tlate<int> uvector_int_view;
  /// uvector_int_base typedef
  typedef uvector_base_tlate<int> uvector_int_base;
  /// uvector_int_const_view typedef
  typedef uvector_const_view_tlate<int> uvector_int_const_view;
  /// uvector_int_array typedef
  typedef uvector_array_tlate<int> uvector_int_array;
  /// uvector_int_subvector typedef
  typedef uvector_subvector_tlate<int> uvector_int_subvector;
  /// uvector_int_const_array typedef
  typedef uvector_const_array_tlate<int> uvector_int_const_array;
  /// uvector_int_const_subvector typedef
  typedef uvector_const_subvector_tlate<int> uvector_int_const_subvector;

  /// uvector_size_t typedef
  typedef uvector_tlate<size_t> uvector_size_t;
  /// uvector_size_t_view typedef
  typedef uvector_view_tlate<size_t> uvector_size_t_view;
  /// uvector_size_t_base typedef
  typedef uvector_base_tlate<size_t> uvector_size_t_base;
  /// uvector_size_t_const_view typedef
  typedef uvector_const_view_tlate<size_t> uvector_size_t_const_view;
  /// uvector_size_t_array typedef
  typedef uvector_array_tlate<size_t> uvector_size_t_array;
  /// uvector_size_t_subvector typedef
  typedef uvector_subvector_tlate<size_t> uvector_size_t_subvector;
  /// uvector_size_t_const_array typedef
  typedef uvector_const_array_tlate<size_t> uvector_size_t_const_array;
  /// uvector_size_t_const_subvector typedef
  typedef uvector_const_subvector_tlate<size_t> 
    uvector_size_t_const_subvector;

  /** \brief A operator for simple vector output

      This outputs all of the vector elements. All of these are
      separated by one space character, though no trailing space or \c
      endl is sent to the output. If the vector is empty, nothing
      is done.
  */
  template<class data_t>
    std::ostream &operator<<
    (std::ostream &os, 
     const uvector_const_view_tlate<data_t> &v) {
    if (v.size()>0) {
      for(size_t i=0;i<v.size()-1;i++) {
	os << v[i] << ' ';
      }
      os << v[v.size()-1];
    }
    return os;
  }
  
  /** \brief A simple class to provide an \c allocate() function
      for \ref uvector

      See the \ref vec_for_tlate_subsect section of the \o2
      User's guide for more information.
  */
  class uvector_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(uvector &o, size_t i) { o.allocate(i); }
    /// Free memory
    void free(uvector &o) { o.free(); }
  };

  /** \brief A simple class to provide an \c allocate() function
      for \ref uvector_int

      See the \ref vec_for_tlate_subsect section of the \o2
      User's guide for more information.
  */
  class uvector_int_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(uvector_int &o, size_t i) { o.allocate(i); }
    /// Free memory
    void free(uvector_int &o) { o.free(); }
  };

  /** \brief A simple class to provide an \c allocate() function
      for \ref uvector_size_t

      See the \ref vec_for_tlate_subsect section of the \o2
      User's guide for more information.
  */
  class uvector_size_t_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(uvector_size_t &o, size_t i) { o.allocate(i); }
    /// Free memory
    void free(uvector_size_t &o) { o.free(); }
  };

  /** \brief A vector with unit-stride where the memory allocation is 
      performed in the constructor
  */
  template<size_t N=0> class ufvector : public uvector_tlate<double> {

  public:

  ufvector() : uvector_tlate<double>(N) {
    }

  };
  
#ifndef DOXYGENP
}
#endif

#endif

