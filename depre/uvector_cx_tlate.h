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
#ifndef O2SCL_UVECTOR_CX_TLATE_H
#define O2SCL_UVECTOR_CX_TLATE_H

/** \file uvector_cx_tlate.h
    \brief File for definitions of complex unit-stride vectors
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include <gsl/gsl_vector.h>

#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>

#ifndef DOXYGENP
namespace o2scl {
#endif
  
  // Forward declarations so that we can friend these classes in
  // uvector_const_view_tlate below
  template<class data_t, class complex_t> class uvector_cx_subvector_tlate;

  /** \brief A vector view of complex numbers with unit stride

      \future Write lookup() method, and possibly an erase() method.
  */
  template<class data_t, class complex_t> class uvector_cx_view_tlate {
    
#ifndef DOXYGEN_INTERNAL

  protected:

    friend class uvector_cx_subvector_tlate<data_t,complex_t>;

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
    uvector_cx_view_tlate(const uvector_cx_view_tlate &v) {
      data=v.data;
      sz=v.sz;
      owner=0;      
    }
    
    /// Copy constructor - create a new view of the same vector
    uvector_cx_view_tlate& operator=(const uvector_cx_view_tlate &v) {
      data=v.data;
      sz=v.sz;
      owner=0;      

      return *this;
    }
    //@}

    ~uvector_cx_view_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    complex_t &operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::operator[]. Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(complex_t *)(data);
      }
#endif
      return *(complex_t *)(&data[i*2]);
    }
    
    /** \brief Array-like indexing 
    */
    const complex_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::operator[] const. Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(const complex_t *)(data);
      }
#endif
      return *(const complex_t *)(&data[i*2]);
    }
    
    /** \brief Array-like indexing 
    */
    complex_t &operator()(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::operator(). Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(complex_t *)(data);
      }
#endif
      return *(complex_t *)(&data[i*2]);
    }
    
    /** \brief Array-like indexing 
    */
    const complex_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::operator() const. Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(const complex_t *)(data);
      }
#endif
      return *(const complex_t *)(&data[i*2]);
    }
    
    /** \brief Get (with optional range-checking) */
    complex_t get(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in uvector_cx_view_tlate::get(). Size: "+itos(sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return *(complex_t *)(data);
      }
#endif
      return *(complex_t *)(&data[i*2]);
    }

    /** \brief Get real part (with optional range-checking) */
    data_t real(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::get(). Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(data_t *)(data);
      }
#endif
      return *(data_t *)(&data[i*2]);
    }
    
    /** \brief Get imaginary part (with optional range-checking) */
    data_t imag(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::get(). Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(data_t *)(data);
      }
#endif
      return *(data_t *)(&data[i*2+1]);
    }
    
    /** \brief Get (with optional range-checking) */
    std::complex<data_t> get_stl(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::get_stl(). Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	std::complex<data_t> zero(0,0);
	return zero;
      }
#endif
      return *(std::complex<data_t> *)(&data[i*2]);
    }

    /** \brief Get pointer (with optional range-checking) */
    complex_t *get_ptr(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::get_ptr(). Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return (complex_t *)(data);
      }
#endif
      return (complex_t *)(&data[i*2]);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const complex_t *get_const_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in uvector_cx_view_tlate::get_const_ptr(). Size: "+
		 itos(sz)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return (const complex_t *)(data);
      }
#endif
      return (const complex_t *)(&data[i*2]);
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, const complex_t &val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in uvector_cx_view_tlate::set(). Size: "+
		   itos(sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
      }
#endif
      data[2*i]=val.dat[0];
      data[2*i+1]=val.dat[1];
      return;
    }

    /** \brief Set (with optional range-checking) */
    void set(size_t i, double r, double c) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=sz) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in uvector_cx_view_tlate::set(). Size: "+
		   itos(sz)+
		   " (index should be less than size).").c_str(),gsl_eindex);
      }
#endif
      data[2*i]=r;
      data[2*i+1]=c;
      return;
    }

    /** \brief Set all of the value to be the value \c val */
    void set_all(const complex_t &val) {
      double re=val.dat[0];
      double im=val.dat[1];
      for(size_t i=0;i<sz;i++) {
	data[2*i]=re;
	data[2*i+1]=im;
      }
      return;
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

    /// \name Arithmetic
    //@{
    /** \brief operator+= */
    uvector_cx_view_tlate<data_t,complex_t> &operator+=
      (const uvector_cx_view_tlate<data_t,complex_t> &x) {
      size_t lsz=x.sz;
      if (lsz>sz) lsz=sz;
      for(size_t i=0;i<lsz;i++) (*this)[i]+=x[i];
      
      return *this;
    }
    
    /** \brief operator-= */
    uvector_cx_view_tlate<data_t,complex_t> &operator-=
      (const uvector_cx_view_tlate<data_t,complex_t> &x) {
      size_t lsz=x.sz;
      if (lsz>sz) lsz=sz;
      for(size_t i=0;i<lsz;i++) (*this)[i]-=x[i];
      
      return *this;
    }
    
    /** \brief operator+= */
    uvector_cx_view_tlate<data_t,complex_t> &operator+=(const data_t &y) {
      for(size_t i=0;i<sz;i++) (*this)[i]+=y;

      return *this;
    }

    /** \brief operator-= */
    uvector_cx_view_tlate<data_t,complex_t> &operator-=(const data_t &y) {
      for(size_t i=0;i<sz;i++) (*this)[i]-=y;

      return *this;
    }

    /** \brief operator*= */
    uvector_cx_view_tlate<data_t,complex_t> &operator*=(const data_t &y) {
      for(size_t i=0;i<sz;i++) (*this)[i]*=y;

      return *this;
    }
    
    /** \brief Norm */
    data_t norm() const {
      data_t result=0;
      for(size_t i=0;i<sz;i++) {
	result+=(*this)[i]*(*this)[i];
      }
      return sqrt(result);
    }
    //@}
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor provided for use by 
	uvector_cx_tlate(const uvector_cx_tlate &v)
    */
    uvector_cx_view_tlate() {};

#endif
    
  };

  /** \brief A vector of double-precision numbers with unit stride
  
      There is also an << operator for this class documented
      "Functions" section of \ref uvector_tlate.h.
  */
  template<class data_t, class complex_t> class uvector_cx_tlate : 
  public uvector_cx_view_tlate<data_t,complex_t> {
  public:

    /// \name Standard constructor
    //@{
    /** \brief Create an uvector_cx of size \c n with owner as 'true' */
    uvector_cx_tlate(size_t n=0) {
      // This must be set to 1 even if n=0 so that future
      // calls to operator= work properly
      this->owner=1;
      if (n>0) {
	this->data=new data_t[2*n];
	if (this->data) {
	  this->sz=n;
	} else {
	  O2SCL_ERR("No memory for data in uvector_cx_tlate constructor",
		  gsl_enomem);
	}
      } else {
	this->sz=0;
      }
    }
    //@}
    
    /// \name Copy constructors
    //@{
    /// Deep copy constructor - allocate new space and make a copy
    uvector_cx_tlate(const uvector_cx_tlate &v) : 
      uvector_cx_view_tlate<data_t,complex_t>() {
      size_t n=v.sz;
      if (n>0) {
	this->data=new data_t[2*n];
	if (this->data) {
	  this->sz=n;
	  this->owner=1;
	  for(size_t i=0;i<n;i++) {
	    this->data[2*i]=v[i].dat[0];
	    this->data[2*i+1]=v[i].dat[1];
	  }
	} else {
	  O2SCL_ERR("No memory for data in uvector_cx_tlate constructor",
		  gsl_enomem);
	}
      } else {
	this->sz=0;
      }
    }
      
      /// Deep copy constructor - allocate new space and make a copy
      uvector_cx_tlate(const uvector_cx_view_tlate<data_t,complex_t> &v) : 
	uvector_cx_view_tlate<data_t,complex_t>() {
	size_t n=v.sz;
	if (n>0) {
	  this->data=new data_t[2*n];
	  if (this->data) {
	    this->sz=n;
	    this->owner=1;
	    for(size_t i=0;i<n;i++) {
	      this->data[2*i]=v[i].dat[0];
	      this->data[2*i+1]=v[i].dat[1];
	    }
	  } else {
	    O2SCL_ERR("No memory for data in uvector_cx_tlate constructor",
		    gsl_enomem);
	  }
	} else {
	  this->sz=0;
	}
      }

	/** \brief Deep copy constructor - if owner is true, allocate space and
	    make a new copy, otherwise, just copy into the view
	*/
	uvector_cx_tlate& operator=(const uvector_cx_tlate &v) {

      // Check for self-assignment
      if (this==&v) return *this;

	  size_t sz2=v.sz;
	  if (this->owner) {
	    allocate(sz2);
	  } else {
	    if (this->sz!=sz2) {
	      O2SCL_ERR("Sizes don't match in uvector_cx_tlate::operator=()",
		      gsl_ebadlen);
	      return *this;
	    }
	  }
	  for(size_t i=0;i<sz2;i++) {
	    this->data[2*i]=v[i].dat[0];
	    this->data[2*i+1]=v[i].dat[1];
	  }
	  return *this;
	}

	/** \brief Deep copy constructor - if owner is true, allocate space and
	    make a new copy, otherwise, just copy into the view
	*/
	uvector_cx_tlate& operator=
	  (const uvector_cx_view_tlate<data_t,complex_t> &v) {

      // Check for self-assignment
      if (this==&v) return *this;

	  size_t sz2=v.size();
	  if (this->owner) {
	    allocate(sz2);
	  } else {
	    if (this->sz!=sz2) {
	      O2SCL_ERR("Sizes don't match in uvector_cx_tlate::operator=()",
		      gsl_ebadlen);
	      return *this;
	    }
	  }
	  for(size_t i=0;i<sz2;i++) {
	    this->data[2*i]=v[i].dat[0];
	    this->data[2*i+1]=v[i].dat[1];
	  }
	  return *this;
	}
	//@}

	~uvector_cx_tlate() {
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
	    this->data=new data_t[2*nsize];
	    this->sz=nsize;
	    this->owner=1;
	    return;
	  }
	  // Neither size was zero
	  data_t *newp=new data_t[2*nsize];
	  if (nsize<this->sz) {
	    vector_copy(2*nsize,this->data,newp);
	  } else {
	    vector_copy(2*this->sz,this->data,newp);
	  }
	  if (this->owner==1) {
	    delete[] this->data;
	  }
	  this->data=newp;
	  this->sz=nsize;
	  this->owner=1;
	  return;
	}

	/** \brief Allocate memory for size \c n after freeing any memory 
	    presently in use
	*/
	void allocate(size_t nsize) {
	  if (this->sz>0) free();

	  if (nsize>0) {
	    this->data=new data_t[2*nsize];
	    this->sz=nsize;
	    this->owner=1;
	  } else {
	    O2SCL_ERR_RET("Zero size in uvector_cx::allocate()",gsl_einval);
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

  };

  /** \brief Create a vector from an array 
   */
  template<class data_t, class complex_t> class uvector_cx_array_tlate : 
  public uvector_cx_view_tlate<data_t,complex_t> {
  public:
    /** \brief Create a vector from \c dat with size \c n */
    uvector_cx_array_tlate(size_t n, data_t *dat) {
      if (n>0) {
	this->data=dat;
	this->sz=n;
	this->owner=0;
      }
    }
  };

  /** \brief Create a vector from a subvector of another
   */
  template<class data_t, class complex_t> class uvector_cx_subvector_tlate : 
  public uvector_cx_view_tlate<data_t,complex_t> {
  public:
    /** \brief Create a vector from \c orig      */
    uvector_cx_subvector_tlate(uvector_cx_view_tlate<data_t,complex_t> &orig, 
			       size_t offset, size_t n) {
      if (offset+n-1<orig.size()) {
	this->data=orig.data+2*offset;
	this->sz=n;
	this->owner=0;
      } else {
	this->sz=0;
	this->data=0;
	O2SCL_ERR("Subvector failed in uvector_cx_sub_view().",1);
      }
    }
  };

  /** \brief Create a vector from an array 
   */
  template<class data_t, class complex_t> class uvector_cx_const_array_tlate :
  public uvector_cx_view_tlate<data_t,complex_t> {
  public:
    /** \brief Create a vector from \c dat with size \c n */
    uvector_cx_const_array_tlate(size_t n, const data_t *dat) {
      if (n>0) {
	// We have to do an explicit cast here, but we prevent the
	// user from changing the data.
	this->data=(data_t *)dat;
	this->sz=n;
	this->owner=0;
      }
    }
    
    ~uvector_cx_const_array_tlate() {};
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /** \name These are inaccessible to ensure the vector is \c const.
     */
    //@{
    data_t &operator[](size_t i) { return this->data[0]; }
    data_t &operator()(size_t i) { return this->data[0]; }
    data_t *get_ptr(size_t i) { return 0; }
    void set(size_t i, data_t val) { return; }
    void set_all(double val) { return; }
    uvector_cx_view_tlate<data_t,complex_t> &operator+=
      (const uvector_cx_view_tlate<data_t,complex_t> &x) {
      return *this;
    }
    uvector_cx_view_tlate<data_t,complex_t> &operator-=
      (const uvector_cx_view_tlate<data_t,complex_t> &x) {
      return *this;
    }
    uvector_cx_view_tlate<data_t,complex_t> &operator*=(const data_t &y) {
      return *this;
    }
    //@}

#endif

  };

  /** \brief Create a vector from a subvector of another
   */
  template<class data_t, class complex_t> 
    class uvector_cx_const_subvector_tlate :
  public uvector_cx_view_tlate<data_t,complex_t> {
  public:
    /** \brief Create a vector from \c orig 
     */
    uvector_cx_const_subvector_tlate
      (const uvector_cx_view_tlate<data_t,complex_t> &orig, 
       size_t offset, size_t n) {
      if (offset+n-1<orig.sz) {
	this->data=orig.data+2*offset;
	this->sz=n;
	this->owner=0;
      } else {
	this->sz=0;
	O2SCL_ERR("Subvector failed in uvector_cx_subvector().",1);
      }
    }

#ifndef DOXYGENP

  protected:

    /** \name Ensure \c const by hiding non-const members
     */
    //@{
    data_t &operator[](size_t i) { return this->data[0]; }
    data_t &operator()(size_t i) { return this->data[0]; }
    data_t *get_ptr(size_t i) { return 0; }
    void set(size_t i, data_t val) { return; }
    void set_all(double val) { return; }
    uvector_cx_view_tlate<data_t,complex_t> &operator+=
      (const uvector_cx_view_tlate<data_t,complex_t> &x) {
      return *this;
    }
    uvector_cx_view_tlate<data_t,complex_t> &operator-=
      (const uvector_cx_view_tlate<data_t,complex_t> &x) {
      return *this;
    }
    uvector_cx_view_tlate<data_t,complex_t> &operator*=(const data_t &y) {
      return *this;
    }
    //@}

#endif

  };

  /// uvector_cx typedef
  typedef uvector_cx_tlate<double,gsl_complex> uvector_cx;
  /// uvector_cx_view typedef
  typedef uvector_cx_view_tlate<double,gsl_complex> uvector_cx_view;
  /// uvector_cx_array typedef
  typedef uvector_cx_array_tlate<double,gsl_complex> uvector_cx_array;
  /// uvector_cx_subvector typedef
  typedef uvector_cx_subvector_tlate<double,gsl_complex> uvector_cx_subvector;
  /// uvector_cx_const_array typedef
  typedef uvector_cx_const_array_tlate<double,gsl_complex> 
    uvector_cx_const_array;
  /// uvector_cx_const_subvector typedef
  typedef uvector_cx_const_subvector_tlate<double,gsl_complex> 
    uvector_cx_const_subvector;
  
  /** \brief A operator for naive vector output

      This outputs all of the vector elements. All of these are
      separated by one space character, though no trailing space or \c
      endl is sent to the output.
  */
  template<class data_t, class complex_t>
    std::ostream &operator<<
    (std::ostream &os, 
     const uvector_cx_view_tlate<data_t,complex_t> &v) {
    if (v.size()>0) {
      for(size_t i=0;i<v.size()-1;i++) {
	os << v[i] << ' ';
      }
      os << v[v.size()-1];
    } else {
      os << "<empty>";
    }
    return os;
  }

#ifndef DOXYGENP
}
#endif

#endif

