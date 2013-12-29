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
#ifndef O2SCL_OVECTOR_CX_TLATE_H
#define O2SCL_OVECTOR_CX_TLATE_H

/** \file ovector_cx_tlate.h
    \brief File for definitions of complex vectors
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <complex>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>

#include <o2scl/err_hnd.h>
#include <o2scl/ovector_tlate.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief A vector view of double-precision numbers
   */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    class ovector_cx_view_tlate : public vparent_t {
    
    public:
    
    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same vector
    ovector_cx_view_tlate(const ovector_cx_view_tlate &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same vector
    ovector_cx_view_tlate& operator=(const ovector_cx_view_tlate &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      

      return *this;
    }
    //@}

    ~ovector_cx_view_tlate() {};
    
    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    complex_t &operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::operator[]. Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(complex_t *)vparent_t::data;
      }
#endif
      return *(complex_t *)(vparent_t::data+i*vparent_t::stride*2);
    }
    
    /** \brief Array-like indexing 
    */
    const complex_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::operator[] const. Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(complex_t *)vparent_t::data;
      }
#endif
      return *(complex_t *)(vparent_t::data+i*vparent_t::stride*2);
    }
    
    /** \brief Array-like indexing 
    */
    complex_t &operator()(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::operator(). Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(complex_t *)vparent_t::data;
      }
#endif
      return *(complex_t *)(vparent_t::data+i*vparent_t::stride*2);
    }
    
    /** \brief Array-like indexing 
    */
    const complex_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::operator() const. Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(complex_t *)vparent_t::data;
      }
#endif
      return *(complex_t *)(vparent_t::data+i*vparent_t::stride*2);
    }
       
    /** \brief Get (with optional range-checking) */
    complex_t get(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::get(). Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(complex_t *)vparent_t::data;
      }
#endif
      return *(complex_t *)(vparent_t::data+i*vparent_t::stride*2);
    }
    
    /** \brief Get real part (with optional range-checking) */
    data_t real(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::real(). Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(data_t *)vparent_t::data;
      }
#endif
      return *(data_t *)(vparent_t::data+i*vparent_t::stride*2);
    }
    
    /** \brief Get imaginary part (with optional range-checking) */
    data_t imag(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::imag(). Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return *(data_t *)vparent_t::data;
      }
#endif
      return *(data_t *)(vparent_t::data+i*vparent_t::stride*2+1);
    }
    
    /** \brief Get STL-like complex number (with optional range-checking) */
    std::complex<data_t> get_stl(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::get_stl(). Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	std::complex<data_t> zero(0,0);
	return zero;
      }
#endif
      return *(std::complex<data_t> *)(vparent_t::data+i*vparent_t::stride*2);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    complex_t *get_ptr(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::get_ptr(). Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return (complex_t *)vparent_t::data;
      }
#endif
      return (complex_t *)(vparent_t::data+i*vparent_t::stride*2);
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const complex_t *get_const_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		 +" in ovector_cx_view_tlate::get_const_ptr(). Size: "+
		 itos(vparent_t::size)+
		 " (index should be less than size).").c_str(),gsl_eindex);
	return (complex_t *)vparent_t::data;
      }
#endif
      return (const complex_t *)(vparent_t::data+i*vparent_t::stride*2);
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, const complex_t &val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_cx_view_tlate::set(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
      }
#endif
      vparent_t::data[i*vparent_t::stride*2]=GSL_REAL(val);
      vparent_t::data[i*vparent_t::stride*2+1]=GSL_IMAG(val);
      return;
    }

    /** \brief Set (with optional range-checking) */
    void set_stl(size_t i, const std::complex<data_t> &d) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_cx_view_tlate::set_stl(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
      }
#endif
      vparent_t::data[i*vparent_t::stride*2]=d.real();
      vparent_t::data[i*vparent_t::stride*2+1]=d.imag();
      return;
    }

    /** \brief Set (with optional range-checking) */
    void set(size_t i, data_t vr, data_t vi) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_cx_view_tlate::set(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
      }
#endif
      vparent_t::data[i*vparent_t::stride*2]=vr;
      vparent_t::data[i*vparent_t::stride*2+1]=vi;
      return;
    }

    /** \brief Set all of the value to be the value \c val */
    void set_all(const complex_t &g) {
      for(size_t i=0;i<vparent_t::size;i++) {
	data_t gr=GSL_REAL(g), gi=GSL_IMAG(g);
	vparent_t::data[i*vparent_t::stride]=gr;
	vparent_t::data[i*vparent_t::stride+1]=gi;
      }
      return;
    }

    /** \brief Method to return vector size 
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t size() const {
      return vparent_t::size;
    }

    /** \brief Method to return vector stride

	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t stride() const {
      return vparent_t::stride;
    }
    //@}

    /// \name Other methods
    //@{
    /// Return true if this object owns the data it refers to
    bool is_owner() const {
      if (vparent_t::owner==1) return true;
      return false;
    }
    //@}

    /// \name Arithmetic
    //@{
    /** \brief operator+= */
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator+=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &x) {
      size_t lsize=x.size();
      if (lsize>vparent_t::size) lsize=vparent_t::size;
      for(size_t i=0;i<lsize;i++) {
	vparent_t::data[i*vparent_t::stride*2]+=x.data[i*x.stride()*2];
	vparent_t::data[i*vparent_t::stride*2+1]+=x.data[i*x.stride()*2+1];
      }
      
      return *this;
    }
    
    /** \brief operator-= */
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator-=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &x) {
      size_t lsize=x.size();
      if (lsize>vparent_t::size) lsize=vparent_t::size;
      for(size_t i=0;i<lsize;i++) {
	vparent_t::data[i*vparent_t::stride*2]-=x.data[i*x.stride()*2];
	vparent_t::data[i*vparent_t::stride*2+1]-=x.data[i*x.stride()*2+1];
      }
      
      return *this;
    }
    
    /** \brief operator+= */
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator+=
      (const complex_t &x) {
      data_t re=*(data_t *)x;
      data_t im=*(((data_t *)x)+1);
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride*2]+=re;
	vparent_t::data[i*vparent_t::stride*2+1]+=im;
      }
      
      return *this;
    }

    /** \brief operator-= */
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator-=
      (const complex_t &x) {
      data_t re=*(data_t *)x;
      data_t im=*(((data_t *)x)+1);
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride*2]-=re;
	vparent_t::data[i*vparent_t::stride*2+1]-=im;
      }
      
      return *this;
    }

    /** \brief operator*= */
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator*=
      (const complex_t &x) {
      data_t re=*(data_t *)x;
      data_t im=*(((data_t *)x)+1);
      for(size_t i=0;i<vparent_t::size;i++) {
	data_t rold=vparent_t::data[i*vparent_t::stride*2];
	data_t iold=vparent_t::data[i*vparent_t::stride*2+1];
	vparent_t::data[i*vparent_t::stride*2]=re*rold-im*iold;
	vparent_t::data[i*vparent_t::stride*2+1]=re*iold+im*rold;
      }
      
      return *this;
    }

    /** \brief operator+= */
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator+=
      (const data_t &x) {
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride*2]+=x;
      }
      
      return *this;
    }

    /** \brief operator-= */
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator-=
      (const data_t &x) {
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride*2]-=x;
      }
      
      return *this;
    }

    /** \brief operator*= */
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator*=
      (const data_t &x) {
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride*2]*=x;
	vparent_t::data[i*vparent_t::stride*2+1]*=x;
      }
      
      return *this;
    }
    //@}
    
    /// Conjugate the vector
    void conjugate() {
      for(size_t i=0;i<vparent_t::size;i++) {
	data_t tmp=(*this)[i].dat[1];
	(*this)[i].dat[1]=-tmp;
	//vparent_t::data[i*vparent_t::stride*2+1]*=-1.0;
      }
      return;
    }

    /** \brief Complex norm \f$ v^{\dagger} v \f$ */
    data_t norm() const {
      data_t result=0;
      for(size_t i=0;i<vparent_t::size;i++) {
	data_t re=vparent_t::data[i*vparent_t::stride*2];
	data_t im=vparent_t::data[i*vparent_t::stride*2+1];
	result+=re*re+im*im;
      }
      return sqrt(result);
    }

#ifndef DOXYGEN_INTERNAL

    protected:

    /** \brief Empty constructor provided for use by 
	ovector_cx_tlate(const ovector_cx_tlate &v) [protected] 
    */
    ovector_cx_view_tlate() {};

#endif
    
  };

  /** \brief A vector of double-precision numbers
  
      If the memory allocation fails, either in the constructor or in
      allocate(), then the error handler will be called, partially
      allocated memory will be freed, and the size will be reset to
      zero. You can test to see if the allocation succeeded using
      something like
      \code
      const size_t n=10;
      ovector_cx x(10);
      if (x.size()==0) cout << "Failed." << endl;
      \endcode

      \todo Add subvector_stride, const_subvector_stride

  */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    class ovector_cx_tlate : 
    public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {

    public:

    /// \name Standard constructor 
    //@{
    /** \brief Create an ovector_cx of size \c n with owner as 'true' */
    ovector_cx_tlate(size_t n=0) {

      vparent_t::data=0;
      vparent_t::size=0;
      vparent_t::stride=0;

      // This must be set to 1 even if n=0 so that future
      // calls to operator= work properly
      vparent_t::owner=1;

      if (n>0) {
	vparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (vparent_t::block) {
	  vparent_t::block->data=(data_t *)malloc(2*n*sizeof(data_t));
	  if (vparent_t::block->data) {
	    vparent_t::block->size=n;
	    vparent_t::data=vparent_t::block->data;
	    vparent_t::size=n;
	    vparent_t::stride=1;
	  } else {
	    std::free(vparent_t::block);
	    O2SCL_ERR("No memory for data in ovector_cx_tlate constructor",
		    gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in ovector_cx_tlate contructor",
		  gsl_enomem);
	}
      }
    }
    //@}
    
    /// \name Copy constructors
    //@{
    /// Deep copy constructor, allocate new space and make a copy
    ovector_cx_tlate(const ovector_cx_tlate &v) : 
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t>() {
      
      vparent_t::data=0;
      vparent_t::size=0;
      vparent_t::stride=0;
      vparent_t::block=0;
      vparent_t::owner=0;
      
      size_t n=v.size();
      if (n>0) {
	vparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (vparent_t::block) {
	  vparent_t::block->data=(data_t *)malloc(2*n*sizeof(data_t));
	  if (vparent_t::block->data) {
	    vparent_t::block->size=n;
	    vparent_t::data=vparent_t::block->data;
	    vparent_t::size=n;
	    vparent_t::stride=1;
	    vparent_t::owner=1;
	    for(size_t i=0;i<n;i++) {
	      vparent_t::data[2*i]=v.data[2*i*v.stride()];
	      vparent_t::data[2*i+1]=v.data[2*i*v.stride()+1];
	    }
	  } else {
	    std::free(vparent_t::block);
	    O2SCL_ERR("No memory for data in ovector_cx_tlate constructor",
		    gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in ovector_cx_tlate contructor",
		  gsl_enomem);
	}
      } else {
	vparent_t::size=0;
	vparent_t::stride=0;
      }
    }
    
    /// Deep copy constructor, allocate new space and make a copy
    ovector_cx_tlate
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &v) : 
      ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t>() {
      size_t n=v.size();
      if (n>0) {
	vparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (vparent_t::block) {
	  vparent_t::block->data=(data_t *)malloc(2*n*sizeof(data_t));
	  if (vparent_t::block->data) {
	    vparent_t::block->size=n;
	    vparent_t::data=vparent_t::block->data;
	    vparent_t::size=n;
	    vparent_t::stride=1;
	    vparent_t::owner=1;
	    for(size_t i=0;i<n;i++) {
	      vparent_t::data[2*i]=v.data[2*i*v.stride()];
	      vparent_t::data[2*i+1]=v.data[2*i*v.stride()+1];
	    }
	  } else {
	    std::free(vparent_t::block);
	    O2SCL_ERR("No memory for data in ovector_cx_tlate constructor",
		    gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in ovector_cx_tlate contructor",
		  gsl_enomem);
	}
      } else {
	vparent_t::size=0;
	vparent_t::stride=0;
      }
    }
    
    /** \brief Deep copy constructor, if owner is true, allocate space and
	make a new copy, otherwise, just copy into the view
    */
    ovector_cx_tlate& operator=(const ovector_cx_tlate &v) {

      // Check for self-assignment
      if (this==&v) return *this;

      size_t size=v.size();
      if (vparent_t::owner) {
	allocate(size);
      } else {
	if (vparent_t::size!=size) {
	  O2SCL_ERR("Sizes don't match in ovector_cx_tlate::operator=()",
		  gsl_ebadlen);
	  return *this;
	}
      }
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[2*i*vparent_t::stride]=v.data[2*i*v.stride()];
	vparent_t::data[2*i*vparent_t::stride+1]=v.data[2*i*v.stride()+1];
      }
      return *this;
    }

    /** \brief Deep copy constructor, if owner is true, allocate space and
	make a new copy, otherwise, just copy into the view
    */
    ovector_cx_tlate& operator=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &v) {

      // Check for self-assignment
      if (this==&v) return *this;

      size_t size=v.size();
      if (vparent_t::owner) {
	allocate(size);
      } else {
	if (vparent_t::size!=size) {
	  O2SCL_ERR("Sizes don't match in ovector_cx_tlate::operator=()",
		  gsl_ebadlen);
	  return *this;
	}
      }
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[2*i*vparent_t::stride]=v.data[2*i*v.stride()];
	vparent_t::data[2*i*vparent_t::stride+1]=v.data[2*i*v.stride()+1];
      }
      return *this;
    }
    //@}

    ~ovector_cx_tlate() {
      if (vparent_t::size>0) {
	if (vparent_t::owner==1) {
	  if (vparent_t::block->size>0) {
	    std::free(vparent_t::block->data);
	  }
	  std::free(vparent_t::block);
	  vparent_t::size=0;
	}
      }
    }
    
    /// \name Memory allocation
    //@{
    /** \brief Allocate memory for size \c n after freeing any memory
	presently in use
    */
    void allocate(size_t nsize) {
      if (vparent_t::size>0) free();

      if (nsize>0) {
	vparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (vparent_t::block) {
	  vparent_t::block->data=(data_t *)malloc(2*nsize*sizeof(data_t));
	  if (vparent_t::block->data) {
	    vparent_t::block->size=nsize;
	    vparent_t::data=vparent_t::block->data;
	    vparent_t::size=nsize;
	    vparent_t::stride=1;
	    vparent_t::owner=1;
	  } else {
	    std::free(vparent_t::block);
	    O2SCL_ERR("No memory for data in ovector_cx_tlate::allocate()",
		      gsl_enomem);
	  }
	} else {
	  O2SCL_ERR("No memory for block in ovector_cx_tlate::allocate()",
		    gsl_enomem);
	}
      } else {
	O2SCL_ERR("Zero size in ovector_cx::allocate()",gsl_einval);
      }
      return;
    }

    /** \brief Free the memory 
    
	This function will safely do nothing if used without first
	allocating memory or if called multiple times in succession.
    */
    void free() {
      if (vparent_t::size>0) {
	if (vparent_t::owner==1) {
	  if (vparent_t::block->size>0) {
	    std::free(vparent_t::block->data);
	  }
	  std::free(vparent_t::block);
	}
	vparent_t::size=0;
	vparent_t::stride=0;
      }
      return;
    }
    //@}
    
    /// \name Other methods
    //@{
    /** \brief Return a gsl vector_cx */
    vparent_t *get_gsl_vector_complex() { return this; }
    
    /** \brief Return a gsl vector_cx */
    const vparent_t *get_gsl_vector_complex_const() const
      { return this; }
    //@}
    
  };

  /** \brief Create a vector from an array 
  */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    class ovector_cx_array_tlate : 
    public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {

    public:

    /** \brief Create a vector from \c dat with size \c n */
    ovector_cx_array_tlate(size_t n, complex_t *dat) {
      if (n>0) {
	vparent_t::block=0;
	vparent_t::data=dat;
	vparent_t::size=n;
	vparent_t::stride=1;
	vparent_t::owner=0;
      }
    }
  };

  /** \brief Create a vector from an array with a stride 
  */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    class ovector_cx_array_stride_tlate : 
    public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {

    public:

    /** \brief Create a vector from \c dat with size \c n and stride \c s 
     */
    ovector_cx_array_stride_tlate(size_t n, complex_t *dat, size_t s) {
      if (n>0) {
	vparent_t::block=0;
	vparent_t::data=dat;
	vparent_t::size=n;
	vparent_t::stride=s;
	vparent_t::owner=0;
      }
    }
  };

  /** \brief Create a vector from a subvector of another
  */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    class ovector_cx_subvector_tlate : 
    public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {

    public:

    /** \brief Create a vector from \c orig */
    ovector_cx_subvector_tlate
      (ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &orig, 
       size_t offset, size_t n) {
      
      vparent_t::size=0;
      vparent_t::stride=0;
      vparent_t::data=0;

      if (offset+n-1<orig.size()) {
	vparent_t::block=0;
	vparent_t::data=orig.data+offset*2*orig.stride();
	vparent_t::size=n;
	vparent_t::stride=orig.stride();
	vparent_t::owner=0;
      } else {
	O2SCL_ERR("Subvector failed in ovector_cx_sub_view().",1);
      }
    }
  };

  /** \brief Create a vector from an array 
  */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    class ovector_cx_const_array_tlate :
    public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {
    public:
    /** \brief Create a vector from \c dat with size \c n */
    ovector_cx_const_array_tlate(size_t n, const complex_t *dat) {
      if (n>0) {
	vparent_t::block=0;
	// We have to do an explicit cast here, but we prevent the
	// user from changing the data.
	vparent_t::data=(data_t *)dat;
	vparent_t::size=n;
	vparent_t::stride=1;
	vparent_t::owner=0;
      }
    }
    
    ~ovector_cx_const_array_tlate() {};
    
#ifndef DOXYGEN_INTERNAL

    protected:
    
    /** \name These are inaccessible to ensure the vector is \c const.
    */
    //@{
    data_t &operator[](size_t i) { return vparent_t::data[0]; }
    data_t &operator()(size_t i) { return vparent_t::data[0]; }
    data_t *get_ptr(size_t i) { return 0; }
    void set(size_t i, data_t val) { return 0; }
    void set_all(double val) { return 0; }
    vparent_t *get_gsl_vector() { return 0; }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator+=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &x) {
      return *this;
    }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator-=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &x) {
      return *this;
    }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator*=
      (const data_t &y) {
      return *this;
    }
    //@}
      
#endif

  };

  /** \brief Create a vector from an array_stride 
  */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    class ovector_cx_const_array_stride_tlate : 
    public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {

    public:

    /** \brief Create a vector from \c dat with size \c n */
    ovector_cx_const_array_stride_tlate(size_t n, const complex_t *dat, 
					size_t s) {
      if (n>0) {
	vparent_t::block=0;
	// We have to do an explicit cast here, but we prevent the
	// user from changing the data.
	vparent_t::data=(data_t *)dat;
	vparent_t::size=n;
	vparent_t::stride=s;
	vparent_t::owner=0;
      }
    }

#ifndef DOXYGEN_INTERNAL

    protected:
    
    /** \name These are inaccessible to ensure the vector is \c const.
    */
    //@{
    data_t &operator[](size_t i) { return vparent_t::data[0]; }
    data_t &operator()(size_t i) { return vparent_t::data[0]; }
    data_t *get_ptr(size_t i) { return 0; }
    void set(size_t i, data_t val) { return 0; }
    void set_all(double val) { return 0; }
    vparent_t *get_gsl_vector() { return 0; }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator+=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &x) {
      return *this;
    }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator-=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &x) {
      return *this;
    }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator*=
      (const data_t &y) {
      return *this;
    }
    //@}

#endif

  };

  /** \brief Create a vector from a subvector of another
  */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    class ovector_cx_const_subvector_tlate :
    public ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> {

    public:

    /** \brief Create a vector from \c orig 
     */
    ovector_cx_const_subvector_tlate
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &orig, 
       size_t offset, size_t n) {
      if (offset+n-1<orig.size()) {
	vparent_t::block=0;
	vparent_t::data=orig.data+offset*2*orig.stride();
	vparent_t::size=n;
	vparent_t::stride=orig.stride();
	vparent_t::owner=0;
      } else {
	vparent_t::size=0;
	O2SCL_ERR("Subvector failed in ovector_cx_subvector().",1);
      }
    }

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \name These are inaccessible to ensure the vector is \c const.
    */
    //@{
    data_t &operator[](size_t i) { return vparent_t::data[0]; }
    data_t &operator()(size_t i) { return vparent_t::data[0]; }
    data_t *get_ptr(size_t i) { return 0; }
    void set(size_t i, data_t val) { return 0; }
    void set_all(double val) { return 0; }
    vparent_t *get_gsl_vector() { return 0; }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator+=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &x) {
      return *this;
    }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator-=
      (const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &x) {
      return *this;
    }
    ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &operator*=
      (const data_t &y) {
      return *this;
    }
    //@}

#endif

  };

  /** \brief Create a real vector from the real parts of a complex vector
   */
  template<class data_t, class vparent_t, class block_t, 
    class cvparent_t, class cblock_t, class complex_t> 
    class ovector_cx_real_tlate :
    public ovector_view_tlate<data_t,vparent_t,block_t> {

    public:

    /** \brief Create a real vector from the real parts of a complex vector
     */
    ovector_cx_real_tlate
      (ovector_cx_view_tlate<data_t,cvparent_t,cblock_t,complex_t> &x) {
      vparent_t::block=0;
      vparent_t::data=x.data;
      vparent_t::size=x.size();
      vparent_t::stride=2*x.stride();
      vparent_t::owner=false;
    }
  };

  /** \brief Create a imaginary vector from the imaginary parts of a 
      complex vector
   */
  template<class data_t, class vparent_t, class block_t, 
    class cvparent_t, class cblock_t, class complex_t> 
    class ovector_cx_imag_tlate :
    public ovector_view_tlate<data_t,vparent_t,block_t> {

    public:

    /** \brief Create a imaginary vector from the imaginary parts of a 
	complex vector
     */
    ovector_cx_imag_tlate
      (ovector_cx_view_tlate<data_t,cvparent_t,cblock_t,complex_t> &x) {
      vparent_t::block=0;
      vparent_t::data=x.data+1;
      vparent_t::size=x.size();
      vparent_t::stride=2*x.stride();
      vparent_t::owner=false;
    }
  };
  
  /** \brief A complex vector where the memory allocation is performed in 
      the constructor
  */
#ifdef DOXYGENP
  template<size_t N=0> class ofvector_cx : 
    public ovector_cx_tlate<data_t,vparent_t,block_t,complex_t>
#else
    template<size_t N=0> class ofvector_cx : 
    public ovector_cx_tlate<double,gsl_vector_complex,gsl_block_complex,
    gsl_complex>
#endif
    {
      public:
      ofvector_cx() : ovector_cx_tlate<double,gsl_vector_complex,
      gsl_block_complex,gsl_complex>(N) {
      }
    };

  /// ovector_cx typedef
  typedef ovector_cx_tlate<double,gsl_vector_complex,gsl_block_complex,
    gsl_complex> ovector_cx;
  /// ovector_cx_view typedef
  typedef ovector_cx_view_tlate<double,gsl_vector_complex,gsl_block_complex,
    gsl_complex> ovector_cx_view;
  /// ovector_cx_array typedef
  typedef ovector_cx_array_tlate<double,gsl_vector_complex,gsl_block_complex,
    gsl_complex> ovector_cx_array;
  /// ovector_cx_array_stride typedef
  typedef ovector_cx_array_stride_tlate<double,gsl_vector_complex,
    gsl_block_complex,gsl_complex> 
    ovector_cx_array_stride;
  /// ovector_cx_subvector typedef
  typedef ovector_cx_subvector_tlate<double,gsl_vector_complex,
    gsl_block_complex,gsl_complex> 
    ovector_cx_subvector;
  /// ovector_cx_const_array typedef
  typedef ovector_cx_const_array_tlate<double,gsl_vector_complex,
    gsl_block_complex,gsl_complex> 
    ovector_cx_const_array;
  /// ovector_cx_const_array_stride typedef
  typedef ovector_cx_const_array_stride_tlate<double,gsl_vector_complex,
    gsl_block_complex,gsl_complex> 
    ovector_cx_const_array_stride;
  /// ovector_cx_const_subvector typedef
  typedef ovector_cx_const_subvector_tlate<double,gsl_vector_complex,
    gsl_block_complex,gsl_complex> 
    ovector_cx_const_subvector;
  /// ovector_cx_real typedef
  typedef ovector_cx_real_tlate<double,gsl_vector,gsl_block,
    gsl_vector_complex,gsl_block_complex,gsl_complex> ovector_cx_real;
  /// ovector_cx_imag typedef
  typedef ovector_cx_imag_tlate<double,gsl_vector,gsl_block,
    gsl_vector_complex,gsl_block_complex,gsl_complex> ovector_cx_imag;

  /** \brief Conjugate a vector */
  template<class data_t,class vparent_t,class block_t, class complex_t>
    ovector_cx_tlate<data_t,vparent_t,block_t,complex_t>
    conjugate(ovector_cx_tlate<data_t,vparent_t,block_t,complex_t> &v) {
    ovector_tlate<data_t,vparent_t,block_t> result(v.size());
    for(size_t i=0;i<v.size();i++) {
      double *p=vparent_t::data+i*vparent_t::stride*2+1;
      (*p)=-(*p);
    }
    return result;
  }

  /** \brief A operator for naive vector output

      This outputs all of the vector elements in the form \c (r,i).
      All of these are separated by one space character, though no
      trailing space or \c endl is sent to the output.
  */
  template<class data_t, class vparent_t, class block_t, class complex_t> 
    std::ostream &operator<<
    (std::ostream &os, 
     const ovector_cx_view_tlate<data_t,vparent_t,block_t,complex_t> &v) {
    for(size_t i=0;i<v.size()-1;i++) {
      os << '(' << v.real(i) << ',' << v.imag(i) << ") ";
    }
    os << '(' << v.real(v.size()-1) << ',' 
       << v.imag(v.size()-1) << ") ";
    return os;
  }

#ifndef DOXYGENP
}
#endif

#endif

