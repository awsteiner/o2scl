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
#ifndef O2SCL_OVECTOR_TLATE_H
#define O2SCL_OVECTOR_TLATE_H

/** \file ovector_tlate.h
    \brief Vector classes

    See the \ref vecmat_section section of the User's guide
    for general information about vector classes in \o2 .

    \future Clean up maybe by moving, for example, ovector reverse
    classes to a different header file
    \future Define ovector_uint classes?
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sys.h>

#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>
#include <o2scl/uvector_tlate.h>
#include <o2scl/array.h>
#include <o2scl/vector.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief Norm object for gsl vectors
      
      This object is principally for use inside the ovector templates
      e.g.. \ref ovector_tlate . For most applications, this object
      need not be instantiated directly by the end-user. 

     For an empty vector, this function returns zero and does not
      call the error handler.
  */
  class gsl_vector_norm : public gsl_vector {
    
  public:
    
    /// Compute the norm of a gsl_vector
    double norm() const {
      
      double scale = 0.0;				 
      double ssq = 1.0;				 
      
      if (this->size <= 0) {				 
	return 0.0;				 
      } else if (this->size == 1) {			 
	return fabs((this->data)[0]);
      }						 
      
      for (size_t i = 0; i < this->size; i++) {		 
	const double x = (this->data)[i*this->stride];
	
	if (x != 0.0) {				 
	  const double ax = fabs(x);		 
	  
	  if (scale < ax) {				    
	    ssq = 1.0 + ssq*(scale / ax)*(scale / ax);    
	    scale = ax;					    
	  } else {					    
	    ssq += (ax / scale)*(ax / scale);		    
	  }						    
	}							    
	
      }							    
      
      return scale*std::sqrt(ssq);				    
    }
  };

  /* Forward declaration of ovector_base_tlate (This is used to allow
     the iterator to friend ovector_base; undocumented)
   */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_base_tlate;
  
  /** \brief A const vector view with finite stride

      See the \ref vecmat_section section of the User's guide
      for general information about vector classes in \o2 .
      
      See the explicit instantations of 
      this class, \ref ovector_const_view and
      \ref ovector_int_const_view .
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_const_view_tlate : 
  public vparent_t {
    
  public:

    // Forward declaration of iterator class (undocumented)
    class iterator;

    /** \brief A const iterator for ovectors
	
	\note Experimental.
	
	\todo Default constructor and iterator typedefs
    */
    class const_iterator {

    protected:

      /// Pointer to current
      data_t *dp;

      /// Stride
      size_t stride;

      /// Internally create an iterator directly from a pointer
      const_iterator(data_t *p, size_t s) {
	dp=p;
	stride=s;
      }

      /// Grant access to ovector classes for begin() and end() functions
      friend class ovector_const_view_tlate<data_t,vparent_t,block_t>;

      /// Grant access to ovector classes for begin() and end() functions
      friend class ovector_base_tlate<data_t,vparent_t,block_t>;

    public:
      /// Copy-constructor for a const iterator
      const_iterator(const iterator &it) {
	dp=it.dp;
	stride=it.stride;
      }

      /// Equality
      bool operator==(const const_iterator &it) const {
	return (dp==it.dp);
      }

      /// Inequality
      bool operator!=(const const_iterator &it) const {
	return (dp!=it.dp);
      }

      /// Inequality
      bool operator!=(const iterator &it) const {
	return (dp!=it.dp);
      }

      /// Less than
      bool operator<(const const_iterator &it) const {
	return (dp<it.dp);
      }

      /// Greater than
      bool operator>(const const_iterator &it) const {
	return (dp>it.dp);
      }

      /// Prefix increment
      const_iterator operator++() {
	dp+=stride;
	return *this;
      }

      /// Prefix decrement
      const_iterator operator--() {
	dp-=stride;
	return *this;
      }

      /// Postfix increment
      const_iterator operator++(int) {
	dp+=stride;
	return *this;
      }

      /// Postfix decrement
      const_iterator operator--(int) {
	dp-=stride;
	return *this;
      }

      /// Dereference - return the corresponding vector element
      const data_t operator*() const {
	return *dp;
      }

      /// Move forward
      const_iterator operator+=(size_t n) {
	dp+=n*stride;
	return *this;
      }

      /// Move backwards
      const_iterator operator-=(size_t n) {
	dp-=n*stride;
	return *this;
      }
    };

    /** \brief An iterator for ovectors
	
	\note Experimental.
     */
    class iterator : public const_iterator {

     protected:

      /// Internally create an iterator directly from a pointer
    iterator(data_t *p, size_t s) : const_iterator(p,s) {
      }

      /// Grant access to ovector classes for begin() and end() functions
      friend class ovector_const_view_tlate<data_t,vparent_t,block_t>;

      /// Grant access to ovector classes for begin() and end() functions
      friend class ovector_base_tlate<data_t,vparent_t,block_t>;

    public:

      /// Copy constructor
      iterator(const const_iterator &cit) {
	this->dp=cit.dp;
	this->stride=cit.stride;
      }

      /// Dereference - return the corresponding vector element
      const data_t operator*() const {
	return *(this->dp);
      }

      /// Dereference - return the corresponding vector element
      data_t operator*() {
	return *(this->dp);
      }
    };

    /// An iterator for the beginning of the vector
    const_iterator begin() const {
      return const_iterator(vparent_t::data,vparent_t::stride);
    }

    /// An iterator for the end of the vector
    const_iterator end() const {
      return const_iterator(vparent_t::data+vparent_t::size*vparent_t::stride,
			    vparent_t::stride);
    }

    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same vector
    ovector_const_view_tlate(const ovector_const_view_tlate &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      
    }

    /// Shallow copy constructor - create a new view of the same vector
    ovector_const_view_tlate& operator=(const ovector_const_view_tlate &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      

      return *this;
    }

    /** \brief Shallow copy constructor - view a unit-stride vector
     */
    ovector_const_view_tlate(const uvector_const_view_tlate<data_t> &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=1;
      vparent_t::owner=0;      
    }
    
    /// Shallow copy constructor - view a unit-stride vector
    ovector_const_view_tlate& operator=
      (const uvector_const_view_tlate<data_t> &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=1;
      vparent_t::owner=0;      

      return *this;
    }
    //@}

    ~ovector_const_view_tlate() {};

    /// \name Get methods
    //@{
    /** \brief Array-like indexing (with optional range-checking)
     */
    const data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
        O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
                   +" in ovector_const_view_tlate::operator[]. Size: "+
                   szttos(vparent_t::size)+
                   " (index should be less than size).").c_str(),gsl_eindex);
        return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }
    
    /** \brief Array-like indexing (with optional range-checking)
     */
    const data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
        O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
                   +" in ovector_const_view_tlate::operator(). Size: "+
                   szttos(vparent_t::size)+
                   " (index should be less than size).").c_str(),gsl_eindex);
        return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }

    /** \brief Get (with optional range-checking) */
    data_t get(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_const_view_tlate::get(). Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const data_t *get_const_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
        O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
                   +" in ovector_const_view_tlate::get_const_ptr(). Size: "+
                   szttos(vparent_t::size)+
                   " (index should be less than size).").c_str(),gsl_eindex);
        return vparent_t::data;
      }
#endif
      return (const data_t *)(vparent_t::data+i*vparent_t::stride);
    }
    
    /** \brief Method to return vector size 

        If no memory has been allocated, this will quietly 
        return zero.
    */
    size_t size() const {
      return vparent_t::size;
    }

    /** \brief Method to return capacity

        Analogous to <tt>std::vector<>.capacity()</tt>.
    */
    size_t capacity() const {
      if (vparent_t::block) return vparent_t::block->size;
      return 0;
    }

    /** \brief Method to return vector stride

        If no memory has been allocated, this will quietly return
        zero.
    */
    size_t stride() const {
      return vparent_t::stride;
    }
    //@}

    /// \name Other methods
    //@{
    /** \brief Return true if this object owns the data it refers to

        This can be used to determine if an object is a "vector_view",
        or a "vector". If is_owner() is true, then it is an \ref
        ovector_tlate object.

        If any \o2 class creates a \ref ovector_tlate object in which
        \ref is_owner() returns false, then it is a bug and should be
        reported.
    */
    bool is_owner() const {
      if (vparent_t::owner==1) return true;
      return false;
    }

    /** \brief Exhaustively look through the vector for a
        particular value and return the closest match

        This can only fail if the vector is empty or if \em all of the
        entries in the vector are not finite. In these cases the
        function calls the error handler and returns 0.

        If more than one entry is the same distance from <tt>x0</tt>,
        this function returns the entry with smallest index.
    */
    size_t lookup(const data_t x0) const {
      if (vparent_t::size==0) {
        O2SCL_ERR("Empty vector ovector_const_view_tlate::lookup().",
                  gsl_einval);
        return 0;
      }
      const ovector_const_view_tlate<data_t,vparent_t,block_t> *a=this;
      size_t row=0, i=0, nvar=size();
      while(!gsl_finite((*a)[i]) && i<nvar-1) i++;
      if (i==nvar-1) {
        O2SCL_ERR2("Entire array not finite in ",
                   "ovector_const_view_tlate::lookup()",gsl_einval);
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

    /** \brief Find the maximum element 

        This can only fail if \em all of the entries in the array are
        not finite or if the vector is empty, in which case it calls
        the error handler and returns 0.
    */
    data_t max() const {
      data_t maxval=0;
      if (vparent_t::size>0) {
        bool setb=false;
        for(size_t i=0;i<vparent_t::size;i++) {
          double v=vparent_t::data[i];
          if (gsl_finite(v)) {
            if (setb==false) {
              maxval=v;
              setb=true;
            } else if (v>maxval) {
              maxval=v;
            }
          }
        }
        if (setb==false) {
          O2SCL_ERR("No finite values in ovector_const_view_tlate::max().",
                    gsl_einval);
          return 0.0;
        }
      } else {
        O2SCL_ERR("Empty vector in ovector_const_view_tlate::max().",
                  gsl_einval);
        return 0.0;
      }
      return maxval;
    }

    /** \brief Find the location of the maximum element 

        This can only fail if \em all of the entries in the array are
        not finite or if the vector is empty, in which case it calls
        the error handler and returns 0.
    */
    size_t max_index() const {

      if (vparent_t::size>0) {

	// These next two variables don't need initialization, but
	// we do so here to avoid warnings. These initial values
	// are essentially ignored by the code below.
	data_t maxval=vparent_t::data[0];
	size_t loc=0;
	
        bool setb=false;
        for(size_t i=0;i<vparent_t::size;i++) {
          double v=vparent_t::data[i];
          if (gsl_finite(v)) {
            if (setb==false) {
              maxval=v;
              loc=i;
              setb=true;
            } else if (v>maxval) {
              maxval=v;
              loc=i;
            }
          }
        }

        if (setb==false) {
          O2SCL_ERR2("No finite values in ",
		     "ovector_const_view_tlate::max_index().",gsl_einval);
          return 0;
        }

	return loc;

      }
       
      O2SCL_ERR("Empty vector ovector_const_view_tlate::max_index().",
		gsl_einval);
      return 0;
    }

    /** \brief Find the minimum element 

        This can only fail if \em all of the entries in the array are
        not finite or if the vector is empty, in which case it calls
        the error handler and returns 0.
    */
    data_t min() const {
      data_t minval=0;
      if (vparent_t::size>0) {
        bool setb=false;
        for(size_t i=0;i<vparent_t::size;i++) {
          double v=vparent_t::data[i];
          if (gsl_finite(v)) {
            if (setb==false) {
              minval=v;
              setb=true;
            } else if (v<minval) {
              minval=v;
            }
          }
        }
        if (setb==false) {
          O2SCL_ERR("No finite values in ovector_const_view_tlate::min().",
                    gsl_einval);
          return 0.0;
        }
      } else {
        O2SCL_ERR("Empty vector ovector_const_view_tlate::min().",
                  gsl_einval);
        return 0.0;
      }
      return minval;
    }

    /** \brief Find the location of the minimum element 

        This can only fail if \em all of the entries in the array are
        not finite or if the vector is empty, in which case it calls
        the error handler and returns 0.
    */
    size_t min_index() const {
      
      if (vparent_t::size>0) {

	// These next two variables don't need initialization, but
	// we do so here to avoid warnings. These initial values
	// are essentially ignored by the code below.
	data_t maxval=vparent_t::data[0];
	size_t loc=0;
	
        bool setb=false;
        for(size_t i=0;i<vparent_t::size;i++) {
          double v=vparent_t::data[i];
          if (gsl_finite(v)) {
            if (setb==false) {
              maxval=v;
              loc=i;
              setb=true;
            } else if (v<maxval) {
              maxval=v;
              loc=i;
            }
          }
        }

        if (setb==false) {
          O2SCL_ERR2("No finite values in ",
		     "ovector_const_view_tlate::min_index().",gsl_einval);
          return 0;
        }
	
	return loc;
      } 
      
      O2SCL_ERR("Empty vector ovector_const_view_tlate::min_index().",
		gsl_einval);
      return 0;
    }
    //@}



#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor provided for use by 
	ovector_view_tlate(const ovector_view_tlate &v)
    */
    ovector_const_view_tlate() {};

#endif
    
  };

  /** \brief A base class for ovector and ovector_view

      This class provides a base class for ovector and ovector_view,
      mostly useful for creating function arguments which accept
      either ovector or ovector_view. 

      See the \ref vecmat_section section of the User's guide
      for general information about vector classes in \o2 .
      
      See the explicit instantations of 
      this class, \ref ovector_base and
      \ref ovector_int_base .
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_base_tlate : 
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {
    
  public:

    /// Desc
    typedef typename ovector_const_view_tlate<data_t,
      vparent_t,block_t>::const_iterator const_iterator;
    /// Desc
    typedef typename ovector_const_view_tlate<data_t,
      vparent_t,block_t>::iterator iterator;
    
    /// Desc
    const_iterator begin() const {
      return const_iterator(vparent_t::data,vparent_t::stride);
    }
    /// Desc
    const_iterator end() const {
      return const_iterator(vparent_t::data+vparent_t::size*vparent_t::stride,
			    vparent_t::stride);
    }
    /// Desc
    iterator begin() {
      return iterator(vparent_t::data,vparent_t::stride);
    }
    /// Desc
    iterator end() {
      return iterator(vparent_t::data+vparent_t::size*vparent_t::stride,
		      vparent_t::stride);
    }

    /// \name Copy constructors
    //@{

    /// Shallow copy constructor - create a new view of the same vector
  ovector_base_tlate(ovector_base_tlate &v)
    : ovector_const_view_tlate<data_t,vparent_t,block_t>() {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      
    }
    
    /// Shallow copy constructor - create a new view of the same vector
    ovector_base_tlate& operator=(ovector_base_tlate &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      

      return *this;
    }

    //@}

    ~ovector_base_tlate() {};

    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    const data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_base_tlate::operator[]. Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }
    
    /** \brief Array-like indexing with operator()
    */
    const data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_base_tlate::operator(). Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }

    /** \brief Array-like indexing 
    */
    data_t &operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_base_tlate::operator[]. Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }
    
    /** \brief Array-like indexing with operator()
    */
    data_t &operator()(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_base_tlate::operator(). Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_base_tlate::get_ptr(). Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data;
      }
#endif
      return vparent_t::data+i*vparent_t::stride;
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, data_t val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_base_tlate::set(). Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),
		  gsl_eindex);
      }
#endif
      vparent_t::data[i*vparent_t::stride]=val;
      return;
    }

    /** \brief Set all of the value to be the value \c val 

	If the vector is empty, this function does not perform
	any assignment and does not call the error handler.
    */
    void set_all(data_t val) {
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride]=val;
      }
      return;
    }
    //@}

    /** \name Arithmetic */
    //@{
    /** \brief operator+= 
	  
	This operator only operates on elements which are present
	in both vectors, i.e. elements which are present in one vector
	but missing in the other will be ignored. If one of the
	two vectors is empty, this function does nothing and
	does not call the error handler.
    */
    ovector_base_tlate<data_t,vparent_t,block_t> &operator+=
      (const ovector_base_tlate<data_t,vparent_t,block_t> &x) {
      size_t lsize=x.size();
      if (lsize>this->size()) lsize=this->size();
      for(size_t i=0;i<lsize;i++) (*this)[i]+=x[i];
      
      return *this;
    }
    
    /** \brief operator-= 
	
	This operator only operates on elements which are present
	in both vectors, i.e. elements which are present in one vector
	but missing in the other will be ignored. If one of the
	two vectors is empty, this function does nothing and
	does not call the error handler.
    */
    ovector_base_tlate<data_t,vparent_t,block_t> &operator-=
      (const ovector_base_tlate<data_t,vparent_t,block_t> &x) {
      size_t lsize=x.size();
      if (lsize>this->size()) lsize=this->size();
      for(size_t i=0;i<lsize;i++) (*this)[i]-=x[i];
      
      return *this;
    }

    /** \brief operator+= 

	If the vector is empty, this function does not perform
	any modifications and does not call the error handler.
    */
    ovector_base_tlate<data_t,vparent_t,block_t> &operator+=(const data_t &x) {
      for(size_t i=0;i<this->size();i++) (*this)[i]+=x;
      return *this;
    }
    
    /** \brief operator-= 

	If the vector is empty, this function does not perform
	any modifications and does not call the error handler.
    */
    ovector_base_tlate<data_t,vparent_t,block_t> &operator-=(const data_t &x) {
      for(size_t i=0;i<this->size();i++) (*this)[i]-=x;
      return *this;
    }
    
    /** \brief operator*= 

	If the vector is empty, this function does not perform
	any modifications and does not call the error handler.
    */
    ovector_base_tlate<data_t,vparent_t,block_t> &operator*=
      (const data_t &y) {
      for(size_t i=0;i<this->size();i++) (*this)[i]*=y;
	
      return *this;
    }
    //@}
    
    /// \name Dynamic casting to the base type
    //@{
    /** \brief Return a gsl vector */
    vparent_t *get_gsl_vector() { 
      vparent_t *x=dynamic_cast<vparent_t *>(this);
      if (x==0) {
	O2SCL_ERR2("Dynamic cast failed in ovector_base_tlate::",
		   "get_gsl_vector().",gsl_esanity);
      }
      return x; 
    };

    /** \brief Return a \c const gsl vector 
     */
    const vparent_t *get_gsl_vector_const() const { 
      const vparent_t *x=dynamic_cast<const vparent_t *>(this);
      if (x==0) {
	O2SCL_ERR2("Dynamic cast failed in ovector_base_tlate::",
		   "get_gsl_vector_const().",gsl_esanity);
      }
      return x; 
    }
    //@}

  protected:

    /** \brief Empty constructor for use by children

	\comment
	This is private because I want to prevent the casual end-user
	from creating ovector_base objects. End-users should really be
	focused on making normal vectors and views.
	\endcomment
    */
    ovector_base_tlate() {};

  };
  
  /** \brief A vector view with finite stride

      See the \ref vecmat_section section of the User's guide
      for general information about vector classes in \o2 .
      
      See the explicit instantations of 
      this class, \ref ovector_view and
      \ref ovector_int_view .
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_view_tlate : 
  public ovector_base_tlate<data_t,vparent_t,block_t> {
    
  public:

    /// \name Copy constructors
    //@{
    /// Shallow copy constructor - create a new view of the same vector
  ovector_view_tlate(const ovector_view_tlate &v) :
    ovector_base_tlate<data_t,vparent_t,block_t>() {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      
    }
      
    /// Shallow copy constructor - create a new view of the same vector
    ovector_view_tlate& operator=
      (const ovector_view_tlate &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      
      
      return *this;
    }
    
    /// Shallow copy constructor for non-views
  ovector_view_tlate(ovector_base_tlate<data_t,vparent_t,block_t> &v) :
    ovector_base_tlate<data_t,vparent_t,block_t>() {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      
    }
    
    /// Shallow copy constructor for non-views
    ovector_view_tlate& operator=
      (ovector_base_tlate<data_t,vparent_t,block_t> &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;      
      
      return *this;
    }
    
    /// Shallow copy constructor for unit-stride vectors
    ovector_view_tlate
      (uvector_base_tlate<data_t> &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=1;
      vparent_t::owner=0;      
    }
    
    /// Shallow copy constructor for unit-stride vectors
    ovector_view_tlate& operator=
      (uvector_base_tlate<data_t> &v) {
      vparent_t::block=0;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=1;
      vparent_t::owner=0;      
      
      return *this;
    }
    //@}

    /** \name Get and set methods

	\comment
	We have to redefine these here because they are different than
	the forms given in ovector_base. 
	\endcomment
    */
    //@{
    /** \brief Array-like indexing 
    */
    data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_view_tlate::operator[]. Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }
    
    /** \brief Array-like indexing with operator()
    */
    data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_view_tlate::operator(). Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[i*vparent_t::stride];
    }

    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_view_tlate::get_ptr(). Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data;
      }
#endif
      return vparent_t::data+i*vparent_t::stride;
    }
    
    /** \brief Set (with optional range-checking) */
    void set(size_t i, data_t val) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in ovector_view_tlate::set(). Size: "+
		   szttos(vparent_t::size)+
		   " (index should be less than size).").c_str(),
		  gsl_eindex);
      }
#endif
      vparent_t::data[i*vparent_t::stride]=val;
      return;
    }

    /** \brief Set all of the value to be the value \c val 

	If the vector is empty, this function does not perform
	any assignment and does not call the error handler.
    */
    void set_all(double val) const {
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride]=val;
      }
      return;
    }
    //@}
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Empty constructor provided for use by 
	ovector_view_tlate(const ovector_view_tlate &v)
    */
    ovector_view_tlate() {};

#endif

  };

  /** \brief A vector with finite stride

      See the \ref vecmat_section section of the User's guide
      for general information about vector classes in \o2 .
      
      See the explicit instantations of 
      this class, \ref ovector and
      \ref ovector_int .

      There are several useful methods which are defined in the parent
      class, \ref ovector_const_view_tlate . There is also an <<
      operator for this class documented "Functions" section of \ref
      uvector_tlate.h.

      <b>Design notes</b>

      At present, the owner variable can be used to distinguish
      between ovectors and ovector_views: for views, owner is always
      zero, but for normal ovectors, owner is 1 (even for empty
      vectors).

      The reserve(), pop_back(), and push_back() methods require the
      ability to have empty vectors which still have memory allocated
      for them, so the proper test for whether or not memory is
      allocated is to test whether or not 'block' is zero. This means
      we should never have a case where the block is non-zero (i.e.
      there is memory allocated there) but the size of the block is
      zero. This is the test used in the destructor and the free()
      method. The free method also sets block to zero accordingly.

  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_tlate : 
  public ovector_base_tlate<data_t,vparent_t,block_t> {

#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// An internal sanity check to ensure correctness
    void intl_sanity_check(size_t ix) {
      if (vparent_t::owner==0 || (vparent_t::block==0 && vparent_t::size>0) ||
	  (vparent_t::block!=0 && (vparent_t::block->size==0 || 
				   vparent_t::size>vparent_t::block->size))) {
	std::cout << "owner: " << vparent_t::owner 
		  << " block: " << vparent_t::block 
		  << " size: " << vparent_t::size 
		  << " bsize: " << vparent_t::block->size << " "
		  << ix << std::endl;
	O2SCL_ERR("Internal error in ovector_tlate.",gsl_esanity);
      }
      return;
    }

    /** \brief An internal allocate function

	\note This function does nothing if \c nsize is zero. Also,
	this does not free any already allocated memory first (unlike
	the user-interface version allocate()).
    */
    void intl_allocate(size_t nsize) {
      
#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(0);
#endif
      
      if (nsize>0) {
	vparent_t::block=(block_t *)malloc(sizeof(block_t));
	if (vparent_t::block) {
	  vparent_t::block->data=(data_t *)malloc(nsize*sizeof(data_t));
	  if (vparent_t::block->data) {
	    vparent_t::block->size=nsize;
	    vparent_t::data=vparent_t::block->data;
	    vparent_t::size=nsize;
	    vparent_t::stride=1;
	    vparent_t::owner=1;
	  } else {
	    std::free(vparent_t::block);
	    // If the block allocation failed, make sure to 
	    // return block to zero so that destructors know
	    // that memory has not been allocated
	    vparent_t::block=0;
	    O2SCL_ERR(((std::string)"No memory for data (nsize="+
		       szttos(nsize)+") in ovector::allocate()").c_str(),
		      gsl_enomem);
	  }
	} else {
	  O2SCL_ERR(((std::string)"No memory for block (nsize="+
		     szttos(nsize)+") in ovector::allocate()").c_str(),
		    gsl_enomem);
	}
      }

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(1);
#endif

      return;
    }

    /// The internal free function
    void intl_free() {

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(2);
#endif

      if (vparent_t::block!=0 && vparent_t::block->size>0) {
	std::free(vparent_t::block->data);
	vparent_t::block->data=0;
	std::free(vparent_t::block);
	vparent_t::block=0;
      }
      vparent_t::size=0;
      vparent_t::stride=0;
      return;
    }


    /** \brief An internal init() function
	
	This function calls intl_allocate() first, and then copies the
	data from the vector \c v.
    */
    template<class alt_vec_t> void intl_init(size_t n, alt_vec_t &v) {
      
      intl_allocate(n);
      for(size_t i=0;i<n;i++) {
	vparent_t::data[i*vparent_t::stride]=v[i];
      }
      
      return;
    }

    /// An internal assignment function for <tt>operator=()</tt>
    template<class alt_vec_t> ovector_tlate &intl_assign
      (size_t n, alt_vec_t &v) {

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(3);
#endif
      
      // Check for self-assignment
      if (this==&v) return *this;
      
      allocate(n);
      for(size_t i=0;i<n;i++) {
	vparent_t::data[i*vparent_t::stride]=v[i];
      }
      
      return *this;
    }

#endif    

  public:
    
    /// \name Standard constructor
    //@{
    /** \brief Create an ovector of size \c n with owner as 'true' */
    ovector_tlate(size_t n=0) {
      
      // Start with an empty vector (Avoid uninitialized variable errors.)
      vparent_t::data=0;
      vparent_t::size=0;
      vparent_t::stride=0;
      vparent_t::block=0;
      vparent_t::owner=1;

      if (n>0) intl_allocate(n);

    }
    //@}
    
    /// \name Copy constructors
    //@{
    /** \brief Deep copy constructor

	\comment
	Note I think this has to be separate from the constructor
	for ovector_base_tlate objects otherwise the compiler will
	synthesize its own.
	\endcomment
    */
  ovector_tlate(const ovector_tlate &v) :
    ovector_base_tlate<data_t,vparent_t,block_t>() {
      
      // Start with an empty vector (Avoid uninitialized variable errors.)
      vparent_t::data=0;
      vparent_t::size=0;
      vparent_t::stride=0;
      vparent_t::block=0;
      vparent_t::owner=1;

      intl_init(v.size(),v);

    }

    /// Deep copy constructor for generic vectors
    template<class alt_vec_t>
      ovector_tlate(size_t nv, alt_vec_t &v) : 
    ovector_base_tlate<data_t,vparent_t,block_t>() {
      
      // Start with an empty vector (Avoid uninitialized variable errors.)
      vparent_t::data=0;
      vparent_t::size=0;
      vparent_t::stride=0;
      vparent_t::block=0;
      vparent_t::owner=1;

      intl_init(nv,v);

    }

    /// Deep copy constructor for other related vectors
  ovector_tlate(const ovector_const_view_tlate<data_t,vparent_t,block_t> &v) : 
    ovector_base_tlate<data_t,vparent_t,block_t>() {
      
      // Start with an empty vector (Avoid uninitialized variable errors.)
      vparent_t::data=0;
      vparent_t::size=0;
      vparent_t::stride=0;
      vparent_t::block=0;
      vparent_t::owner=1;

      intl_init(v.size(),v);

    }

    /** \brief Deep copy constructor, if owner is true, allocate
	space and make a new copy, otherwise, just copy into the
	view
    */
    ovector_tlate& operator=(const ovector_tlate &v) {

      return intl_assign(v.size(),v);
    }

    /** \brief Deep copy constructor, if owner is true, allocate
	space and make a new copy, otherwise, just copy into the
	view
    */
    ovector_tlate& operator=
      (const ovector_const_view_tlate<data_t,vparent_t,block_t> &v) {
      
      return intl_assign(v.size(),v);
    }

    /** \brief Deep copy constructor, if owner is true, allocate
	space and make a new copy, otherwise, just copy into the
	view
    */
    ovector_tlate& operator=
      (const uvector_const_view_tlate<data_t> &v) {

      return intl_assign(v.size(),v);
    }
    //@}

    ~ovector_tlate() {

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(4);
#endif

      if (vparent_t::owner==1) {
	if (vparent_t::block!=0 && vparent_t::block->size>0) {
	  std::free(vparent_t::block->data);
	}
	if (vparent_t::block!=0) {
	  std::free(vparent_t::block);
	}
      }
      
    }

    /** \brief Reallocate memory for size \c nsize

	If \c nsize is zero, then all memory is deallocated. If the
	new size is less than or equal to the old, then the size is
	updated and the capacity is unchanged. If the new size is
	larger than the old, new memory is allocated and the old
	contents are copied into the new memory (up to the lesser of
	the old and new size).
    */
    void resize(size_t new_size) {

      // New size is zero
      if (new_size==0) {
	intl_free();
      }

      // Old size was zero
      if (vparent_t::size==0) {
	intl_allocate(new_size);
	return;
      }

      // We already have the capacity, so just update the size.
      // This also handles the case where the old and new sizes
      // are the same.
      if (new_size<=vparent_t::size) {
	vparent_t::size=new_size;
	return;
      }

      // Allocate new memory in vparent_t::block
      vparent_t::block->data=(data_t *)malloc
	(new_size*sizeof(data_t));
      vparent_t::block->size=new_size;
      
      // Copy the data to the new memory
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::block->data[i]=vparent_t::data[i];
      }
      
      // Delete the old memory
      std::free(vparent_t::data);
      
      // Reset data to the new memory and update size
      vparent_t::data=vparent_t::block->data;
      vparent_t::size=new_size;
      
      return;
    }
    
    /// \name Memory allocation
    //@{
    /** \brief Allocate memory for size \c n after freeing any memory
	currently in use

	Note that this automatically deallocates any previously
	allocated memory before attempting to allocate more. If this
	allocation fails (i.e. because we ran out of memory) then the
	original vector will still have been deallocated.
    */
    void allocate(size_t nsize) {

      intl_free();
      intl_allocate(nsize);

      return;
    }

    /** \brief Free the memory 
	
	This function will safely do nothing if used without first
	allocating memory or if called multiple times in succession.
    */
    void free() {
      intl_free();
      return; 
    }
    //@}

    /// \name Stack-like operations
    //@{

    /// Add a value to the end of the vector
    void push_back(data_t val) {

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(6);
#endif

      if (vparent_t::block==0) {
	
	// Empty, so make a 1-element vector
	allocate(1);
	vparent_t::data[0]=val;
	  
      } else if (vparent_t::block->size==vparent_t::size) {
	  
	// Allocate new memory in vparent_t::block
	vparent_t::block->data=(data_t *)malloc
	  (2*vparent_t::block->size*sizeof(data_t));
	vparent_t::block->size*=2;
		
	// Copy the data to the new memory
	for(size_t i=0;i<vparent_t::size;i++) {
	  vparent_t::block->data[i]=vparent_t::data[i];
	}
	  
	// Delete the old memory
	std::free(vparent_t::data);
	  
	// Reset data to the new memory
	vparent_t::data=vparent_t::block->data;

	// Add the new value to the end of the array and 
	// increment size
	vparent_t::block->data[vparent_t::size]=val;
	vparent_t::size++;

      } else {

	// Add the new value to the end of the array and 
	// increment size
	vparent_t::block->data[vparent_t::size]=val;
	vparent_t::size++;

      }

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(7);
#endif

      return;
    }

    /** \brief Reserve memory by increasing capacity
    
	Increase the maximum capacity of the vector so that calls
	to push_back() do not need to automatically increase the
	capacity.

	If the argument \c cap is smaller than the present vector size
	given by \ref size(), then this function does nothing and the
	error handler is not called.
    */
    void reserve(size_t cap) {
	    
#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(8);
#endif

      // Do nothing if we are reserving memory for less than the
      // current size
	      
      if (cap>vparent_t::size) {

	bool initially_empty=false;

	// If it's initially_empty, we need to allocate a block
	if (!vparent_t::block) {
	  initially_empty=true;
	  vparent_t::block=(block_t *)malloc(sizeof(block_t));
	}

	// If it's empty, or if the block size is too small,
	// reserve space. On the other hand, if the block size
	// is already large enough, do nothing.
	if (initially_empty || vparent_t::block->size<cap) {

	  // Allocate new memory in vparent_t::block
	  vparent_t::block->data=(data_t *)malloc(cap*sizeof(data_t));
	  vparent_t::block->size=cap;
	  
	  // Copy the data to the new memory
	  for(size_t i=0;i<vparent_t::size;i++) {
	    vparent_t::block->data[i]=vparent_t::data[i];
	  }
	  
	  // Delete the old memory
	  if (!initially_empty) {
	    std::free(vparent_t::data);
	  }
	  
	  // Reset data to the new memory
	  vparent_t::data=vparent_t::block->data;
	}
      }

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(9);
#endif

      return;
    }

    /** \brief Return the last value and shrink the vector size by one
    */
    data_t pop_back() {

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(10);
#endif

      if (vparent_t::size==0) {
	O2SCL_ERR("Attempted to pop an empty vector in ovector::pop_back().",
		  gsl_einval);
	return 0.0;
      }

      // Otherwise, just decrease the size by one and return the last
      // element
      vparent_t::size--;

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(11);
#endif

      return vparent_t::data[vparent_t::size];
    }
    //@}

    /// \name Other methods
    //@{
    /// Remove element with index \c ix and decrease the vector size by one
    void erase(size_t ix) {

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(12);
#endif

      if (ix<vparent_t::size) {
	for(size_t i=ix+1;i<vparent_t::size;i++) {
	  vparent_t::data[i-1]=vparent_t::data[i];
	}
	vparent_t::size--;
      } else {
	O2SCL_ERR("Tried to erase() past end in ovector_tlate::erase().",
		  gsl_einval);
      }

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(13);
#endif

      return;
    }
    
    /** \brief Sort the vector and ensure all elements are unique by 
	removing duplicates
    */
    void sort_unique() {

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(14);
#endif
      
      o2scl::vector_sort<ovector_tlate<data_t,vparent_t,block_t>,data_t >
	(vparent_t::size,*this);
      for(size_t i=1;i<vparent_t::size;i++) {
	if ((*this)[i]==(*this)[i-1]) {
	  erase(i-1);
	  i--;
	}
      }

#if O2SCL_NO_RANGE_CHECK
#else
      intl_sanity_check(15);
#endif

      return;
    }
    //@}
    
  };

  /** \brief Create a vector from an array 
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_array_tlate : 
  public ovector_view_tlate<data_t,vparent_t,block_t> {

  public:

    /** \brief Create a vector from \c dat with size \c n */
    ovector_array_tlate(size_t n, data_t *dat) {
      vparent_t::block=0;
      vparent_t::owner=0;
      if (n>0) {
	vparent_t::data=dat;
	vparent_t::size=n;
	vparent_t::stride=1;
      } else {
	vparent_t::data=0;
	vparent_t::size=0;
	vparent_t::stride=0;
      }
    }

  };

  /** \brief Create a vector from an array with a stride
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_array_stride_tlate : 
  public ovector_view_tlate<data_t,vparent_t,block_t> {

  public:

    /** \brief Create a vector from \c dat with size \c n and stride \c s 
     */
    ovector_array_stride_tlate(size_t n, data_t *dat, size_t s) {
      vparent_t::block=0;
      vparent_t::owner=0;
      if (n>0) {
	vparent_t::data=dat;
	vparent_t::size=n;
	vparent_t::stride=s;
      } else {
	vparent_t::data=0;
	vparent_t::size=0;
	vparent_t::stride=0;
      }
    }

  };

  /** \brief Create a vector from a subvector of another

      The constructor will fail if the original vector is empty, or if
      the user requests a subvector which includes elements beyond the
      end of the original vector.
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_subvector_tlate : 
  public ovector_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector from \c orig 
     */
    ovector_subvector_tlate
      (ovector_base_tlate<data_t,vparent_t,block_t> &orig, 
       size_t offset, size_t n) {
      if (orig.size()>0 && offset+n<orig.size()+1) {
	vparent_t::block=0;
	vparent_t::data=orig.data+offset*orig.stride();
	vparent_t::size=n;
	vparent_t::stride=orig.stride();
	vparent_t::owner=0;
      } else {
	std::string err="Failed because vector is of size "+
	  szttos(orig.size())+" but requested size "+szttos(n)+
	  " with offset "+szttos(offset)+" in ovector_subvector "+
	  "constructor.";
	O2SCL_ERR(err.c_str(),gsl_efailed);
	vparent_t::block=0;
	vparent_t::data=0;
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::owner=0;
      }
    }
  };

  /** \brief Create a const vector from an array 

      The constructor will fail if the size argument 
      \c n is zero or if the pointer \c dat is 0.

      \future Use const_cast here?
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_const_array_tlate :
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {

  public:

    /** \brief Create a vector from \c dat with size \c n */
    ovector_const_array_tlate(size_t n, const data_t *dat) {
      if (n>0 && dat!=0) {
	vparent_t::block=0;
	// We have to do an explicit cast here, but we prevent the
	// user from changing the data.
	vparent_t::data=(data_t *)dat;
	vparent_t::size=n;
	vparent_t::stride=1;
	vparent_t::owner=0;
      } else {
	O2SCL_ERR2("Empty array in ovector_const_array_tlate() ",
		   "constructor.",gsl_efailed);
	vparent_t::block=0;
	vparent_t::data=0;
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::owner=0;
      }
    }
    
    ~ovector_const_array_tlate() {};
    
  };

  /** \brief Create a const vector from an array with a stride 

      The constructor will fail if the size argument \c n is zero, the
      stride \c s is zero, or if the pointer \c dat is 0.
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_const_array_stride_tlate : 
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {

  public:

    /** \brief Create a vector from \c dat with size \c n */
    ovector_const_array_stride_tlate(size_t n, const data_t *dat, size_t s) {
      if (n>0 && s>0 && dat!=0) {
	vparent_t::block=0;
	// We have to do an explicit cast here, but we prevent the
	// user from changing the data.
	vparent_t::data=(data_t *)dat;
	vparent_t::size=n;
	vparent_t::stride=s;
	vparent_t::owner=0;
      } else {
	O2SCL_ERR2("Empty array in ovector_const_array_stride_tlate() ",
		   "constructor.",gsl_efailed);
	vparent_t::block=0;
	vparent_t::data=0;
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::owner=0;
      }
    }

  };

  /** \brief Create a const vector from a subvector of another vector
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_const_subvector_tlate :
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector from \c orig 
     */
    ovector_const_subvector_tlate
      (const ovector_const_view_tlate<data_t,vparent_t,block_t> &orig, 
       size_t offset, size_t n) {
      if (offset+n-1<orig.size()) {
	vparent_t::block=0;
	vparent_t::data=orig.data+offset*orig.stride();
	vparent_t::size=n;
	vparent_t::stride=orig.stride();
	vparent_t::owner=0;
      } else {
	std::string err="Failed because vector is of size "+
	  szttos(orig.size())+" but requested size "+szttos(n)+
	  " with offset "+szttos(offset)+" in ovector_const_subvector "+
	  "constructor.";
	O2SCL_ERR(err.c_str(),gsl_efailed);
	vparent_t::block=0;
	vparent_t::data=0;
	vparent_t::size=0;
	vparent_t::stride=0;
	vparent_t::owner=0;
      }
    }

  };

  /// ovector typedef
  typedef ovector_tlate<double,gsl_vector_norm,gsl_block> ovector;
  /// ovector_base typedef
  typedef ovector_base_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_base;
  /// ovector_view typedef
  typedef ovector_view_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_view;
  /// ovector_const_view typedef
  typedef ovector_const_view_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_const_view;
  /// ovector_array typedef
  typedef ovector_array_tlate<double,gsl_vector_norm,gsl_block> ovector_array;
  /// ovector_array_stride typedef
  typedef ovector_array_stride_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_array_stride;
  /// ovector_subvector typedef
  typedef ovector_subvector_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_subvector;
  /// ovector_const_array typedef
  typedef ovector_const_array_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_const_array;
  /// ovector_const_array_stride typedef
  typedef ovector_const_array_stride_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_const_array_stride;
  /// ovector_const_subvector typedef
  typedef ovector_const_subvector_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_const_subvector;

  /// ovector_int typedef
  typedef ovector_tlate<int,gsl_vector_int,gsl_block_int> ovector_int;
  /// ovector_int_base typedef
  typedef ovector_base_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_base;
  /// ovector_int_view typedef
  typedef ovector_view_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_view;
  /// ovector_int_const_base typedef
  typedef ovector_const_view_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_const_view;
  /// ovector_int_array typedef
  typedef ovector_array_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_array;
  /// ovector_int_array_stride typedef
  typedef ovector_array_stride_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_array_stride;
  /// ovector_int_subvector typedef
  typedef ovector_subvector_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_subvector;
  /// ovector_int_const_array typedef
  typedef ovector_const_array_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_const_array;
  /// ovector_int_const_array_stride typedef
  typedef ovector_const_array_stride_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_const_array_stride;
  /// ovector_int_const_subvector typedef
  typedef ovector_const_subvector_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_const_subvector;

  /** \brief A operator for naive vector output

      This outputs all of the vector elements. All of these are
      separated by one space character, though no trailing space or \c
      endl is sent to the output. If the vector is empty, nothing
      is done.
  */
  template<class data_t, class vparent_t, class block_t> 
    std::ostream &operator<<
    (std::ostream &os, 
     const ovector_const_view_tlate<data_t,vparent_t,block_t> &v) {
    if (v.size()>0) {
      for(size_t i=0;i<v.size()-1;i++) {
	os << v[i] << ' ';
      }
      os << v[v.size()-1];
    }
    return os;
  }
  
  /** \brief A simple class to provide an \c allocate() function
      for \ref ovector

      See the \ref vec_for_tlate_subsect section of the \o2
      User's guide for more information.

      \future Could (or should?) the allocate() functions be adapted to
      return an integer error value?
  */
  class ovector_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(ovector &o, size_t i) { o.allocate(i); }
    /// Free memory
    void free(ovector &o) { o.free(); }
  };

  /** \brief A simple class to provide an \c allocate() function
      for \ref ovector_int

      See the \ref vec_for_tlate_subsect section of the \o2
      User's guide for more information.
  */
  class ovector_int_alloc {
  public:
    /// Allocate \c v for \c i elements
    void allocate(ovector_int &o, size_t i) { o.allocate(i); }
    /// Free memory
    void free(ovector_int &o) { o.free(); }
  };

  /** \brief A vector where the memory allocation is performed in 
      the constructor

      This can be useful, for example to easily make C-style arrays of
      \ref ovector objects with a fixed size. For example,
      \code
      ofvector<8> x[10];
      \endcode
      would mean <tt>x</tt> is a 10-dimensional array of 
      \ref ovector object with initial length 8.

      \future Consider making allocate() and free() functions
      private for this class?
  */
#ifdef DOXYGENP
  template<size_t N=0> class ofvector : 
  public ovector_tlate
#else
    template<size_t N=0> class ofvector : 
  public ovector_tlate<double,gsl_vector,gsl_block>
#endif
    {
    public:
    ofvector() : ovector_tlate<double,gsl_vector,gsl_block>(N) {
      }
    };
  
#ifndef DOXYGENP
}
#endif

#endif

