/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
#ifndef O2SCL_SEARCH_VEC_H
#define O2SCL_SEARCH_VEC_H

/** \file search_vec.h
    \brief File defining \ref o2scl::search_vec and 
    \ref o2scl::search_vec_ext 
*/

#include <iostream>
#include <string>
#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Searching class for monotonic data with caching
      
      A searching class for monotonic vectors. A caching system
      similar to \c gsl_interp_accel is used.
      
      To find the interval containing a value, use find(). If you
      happen to know in advance that the vector is increasing or
      decreasing, then you can use find_inc() or find_dec() instead.
      To ignore the caching and just use straight binary search, you
      can use the functions in \ref vector.h .

      Alternatively, if you just want to find the index with the
      element closest to a specified value, use ordered_lookup(). 
      
      The functions find_inc(), find_dec() and find() are designed to
      return the lower index of an interval containing the desired
      point, and as a result will never return the last index of a
      vector, e.g. for a vector of size <tt>n</tt> they always return
      a number between <tt>0</tt> and <tt>n-2</tt> inclusive. See \ref
      o2scl::search_vec_ext for an alternative.

      \note The behavior of these functions is undefined if some of
      the user-specified data is not finite or not strictly monotonic.
      Two adjacent data points should not be equal. This class does
      not verify that the user-specified data has these properties.

      \note Because results are cached, this class is not thread-safe
      and cannot be used simultaneously by different threads. (This
      holds even for member functions marked const, because the cache 
      data member is marked as mutable.)

      \note This class does not store a copy of the data, but only a
      pointer to it. This means that one can safely modify the data
      after the constructor is called, so long as one does not make
      the vector smaller (as the cache might then point to a value
      outside the new vector) and so long as the new vector is still
      monotonic. Copy constructors are also private to prevent 
      confusing situations which arise when bit-copying pointers. 
  */
  template<class vec_t> class search_vec {

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Storage for the most recent index
	
        \note This is marked mutable to ensure const-correctness is 
	straightforward.
    */
    mutable size_t cache;

    /// The vector to be searched
    const vec_t *v;
    
    /// The vector size
    size_t n;

#endif

  public:

    /** \brief Create a blank searching object
     */
  search_vec() : v(0), n(0) {
    }

    /** \brief Create a searching object with vector \c x of size \c nn
     */
  search_vec(size_t nn, const vec_t &x) : v(&x), n(nn) {
      if (nn<2) {
	std::string str=((std::string)"Vector too small (size=")+
	  o2scl::szttos(nn)+") in search_vec::search_vec().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
      cache=0;
    }

    /** \brief Set the vector to be searched 
     */
    void set_vec(size_t nn, const vec_t &x) {
      if (nn<2) {
	std::string str=((std::string)"Vector too small (size=")+
	  o2scl::szttos(nn)+") in search_vec::set_vec().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
      cache=0;
      v=&x;
      n=nn;
    }
    
    /** \brief Search an increasing or decreasing vector for the
	interval containing <tt>x0</tt>
	
	This function is identical to find_inc() if the data is
	increasing, and find_dec() if the data is decreasing. 
    */
    size_t find(const double x0) const {
#if !O2SCL_NO_RANGE_CHECK
      if (cache>=n) {
	O2SCL_ERR("Cache mis-alignment in search_vec::find().",
		  exc_esanity);
      }
#endif
      if ((*v)[0]<(*v)[n-1]) return find_inc(x0);
      return find_dec(x0);
    }

    /** \brief Search an increasing vector for the interval
	containing <tt>x0</tt>

	This function is a cached version of \ref vector_bsearch_inc()
	, analogous to <tt>gsl_interp_accel_find()</tt>, except
	that it does not internally record cache hits and 
	misses. 

    */
    size_t find_inc(const double x0) const {
      if (x0<(*v)[cache]) {
	cache=vector_bsearch_inc<vec_t,double>(x0,*v,0,cache);
      } else if (x0>=(*v)[cache+1]) {
	cache=vector_bsearch_inc<vec_t,double>(x0,*v,cache,n-1);
      }
#if !O2SCL_NO_RANGE_CHECK
      if (cache>=n) {
	O2SCL_ERR("Cache mis-alignment in search_vec::find_inc().",
		  exc_esanity);
      }
#endif
      return cache;
    }
    
    /** \brief Search a decreasing vector for the interval
	containing <tt>x0</tt>

	This function is a cached version of \ref vector_bsearch_dec()
	.  The operation of this function is undefined if the data is
	not strictly monotonic, i.e. if some of the data elements are
	equal. 
    */
    size_t find_dec(const double x0) const {
      if (x0>(*v)[cache]) {
	cache=vector_bsearch_dec<vec_t,double>(x0,*v,0,cache);
      } else if (x0<=(*v)[cache+1]) {
	cache=vector_bsearch_dec<vec_t,double>(x0,*v,cache,n-1);
      }
#if !O2SCL_NO_RANGE_CHECK
      if (cache>=n) {
	O2SCL_ERR("Cache mis-alignment in search_vec::find_dec().",
		  exc_esanity);
      }
#endif
      return cache;
    }

    /** \brief Find the index of x0 in the ordered array \c x 

	This returns the index i for which x[i] is as close as
	possible to x0 if x[i] is either increasing or decreasing.

	If you have a non-monotonic vector, you can use \ref
	vector_lookup() instead.

	Generally, if there are two adjacent entries with the same
	value, this function will return the entry with the smaller
	index. 

	\future This function just uses the <tt>find</tt> functions
	and then adjusts the answer at the end if necessary. It might
	be possible to improve the speed by rewriting this from
	scratch.
    */
    size_t ordered_lookup(const double x0) const {
      if (n<1) {
	std::string str=((std::string)"Not enough data (n=")+
	  o2scl::szttos(n)+") in search_vec::ordered_lookup().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }

      size_t row;

      if ((*v)[0]<=(*v)[n-1]) {

	// Increasing case

	if (x0>=(*v)[n-1]) {
	  row=n-1;
	} else { 
	  row=find_inc(x0);
	  if (row<n-1 && fabs((*v)[row+1]-x0)<fabs((*v)[row]-x0)) row++;
	}
    
      } else {

	// Decreasing case
    
	if (x0<=(*v)[n-1]) {
	  row=n-1;
	} else {
	  row=find_dec(x0);
	  if (row<n-1 && fabs((*v)[row+1]-x0)<fabs((*v)[row]-x0)) row++;
	}
      }

      return row;
    }

#ifndef DOXYGEN_INTERNAL

  private:

    search_vec<vec_t>(const search_vec<vec_t> &);
    search_vec<vec_t>& operator=(const search_vec<vec_t>&);

#endif

  };

  /** \brief An extended search_vec which is allowed to return 
      the last element
  */
  template<class vec_t> class search_vec_ext {
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Storage for the most recent index
	
        \note This is marked mutable to ensure const-correctness is 
	straightforward.
    */
    mutable size_t cache;

    /// The vector to be searched
    const vec_t *v;
    
    /// The vector size
    size_t n;

#endif

  public:

    /** \brief Create a blank searching object
     */
  search_vec_ext() : v(0), n(0) {
    }

    /** \brief Create a searching object for vector \c x of size 
	\c nn

	\comment
	Note that this constructor does not call the parent
	constructor because that requires nn<2 while this
	class really only requires nn<1.
	\endcomment

	\future Ensure this is fully tested for vectors with
	only one element.
    */
  search_vec_ext(size_t nn, const vec_t &x) : v(&x), n(nn) {
      if (nn<1) {
	std::string str=((std::string)"Vector too small (n=")+
	  o2scl::szttos(nn)+") in search_vec_ext::search_vec_ext().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
      cache=0;
    }
    
    /** \brief Search an increasing or decreasing vector for the interval
	containing <tt>x0</tt>
    */
    size_t find(const double x0) const {
#if !O2SCL_NO_RANGE_CHECK
      if (this->cache>=this->n) {
	O2SCL_ERR("Cache mis-alignment in search_vec_ext::find().",
		  exc_esanity);
      }
#endif
      if ((*this->v)[0]<(*this->v)[this->n-1]) return find_inc(x0);
      return find_dec(x0);
    }

    /** \brief Search an increasing vector for the interval
	containing <tt>x0</tt>
    */
    size_t find_inc(const double x0) const {
      if (x0<(*this->v)[this->cache]) {
	this->cache=vector_bsearch_inc<vec_t,double>
	  (x0,*this->v,0,this->cache);
      } else if (this->cache<this->n-1 && x0>=(*this->v)[this->cache+1]) {
	this->cache=vector_bsearch_inc<vec_t,double>
	  (x0,*this->v,this->cache,this->n);
      }
#if !O2SCL_NO_RANGE_CHECK
      if (this->cache>=this->n) {
	O2SCL_ERR("Cache mis-alignment in search_vec_ext::find_inc().",
		  exc_esanity);
      }
#endif
      return this->cache;
    }
    
    /** \brief Search a decreasing vector for the interval
	containing <tt>x0</tt>
    */
    size_t find_dec(const double x0) const {
      if (x0>(*this->v)[this->cache]) {
	this->cache=vector_bsearch_dec<vec_t,double>
	  (x0,*this->v,0,this->cache);
      } else if (this->cache<this->n-1 && x0<=(*this->v)[this->cache+1]) {
	this->cache=vector_bsearch_dec<vec_t,double>
	  (x0,*this->v,this->cache,this->n);
      }
#if !O2SCL_NO_RANGE_CHECK
      if (this->cache>=this->n) {
	O2SCL_ERR("Cache mis-alignment in search_vec_ext::find_dec().",
		  exc_esanity);
      }
#endif
      return this->cache;
    }

#ifndef DOXYGEN_INTERNAL

  private:

    search_vec_ext<vec_t>(const search_vec_ext<vec_t> &);
    search_vec_ext<vec_t>& operator=(const search_vec_ext<vec_t>&);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
