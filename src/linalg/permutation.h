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
/* permutation/permute_source.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifndef O2SCL_PERMUTATION_H
#define O2SCL_PERMUTATION_H

/** \file permutation.h
    \brief Permutation class and output operator
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A class for representing permutations 

      To apply a permutation to a user-specified
      vector, see \ref o2scl::permutation::apply().

  */
  class permutation {
    
  protected:

    size_t size_;
    boost::numeric::ublas::unbounded_array<size_t> data;

  public:

    /// Create a permutation of size \c dim
    permutation(size_t dim=0) {
      size_=dim;
      if (dim>0) data.resize(dim);
    }

    /// \name Copy constructors
    //@{
    permutation(const permutation &v) {
      size_=v.size();
      if (size_>0) data.resize(size_);
      for(size_t i=0;i<size_;i++) {
	data[i]=v.data[i];
      }
    }
    
    permutation& operator=(const permutation &v) {
      
      // Check for self-assignment
      if (this==&v) return *this;

      size_=v.size();
      if (size_>0) data.resize(size_);
      for(size_t i=0;i<size_;i++) {
	data[i]=v.data[i];
      }
      return *this;
    }
    //@}

    ~permutation() {
    };

    /** \brief Array-like indexing 
    */
    size_t &operator[](size_t i) {
#if !O2SCL_NO_RANGE_CHECK
      if (i>=size_) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		   +" in permutation::operator[]. Size: "+
		   szttos(size_)+
		   " (index should be less than size).").c_str(),exc_eindex);
	return data[0];
      }
#endif
      return data[i];
    }
    
    /** \brief Array-like indexing 
    */
    const size_t &operator[](size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
      if (i>=size_) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		 +" in permutation::operator[]. Size: "+
		 szttos(size_)+
		 " (index should be less than size).").c_str(),exc_eindex);
	return data[0];
      }
#endif
      return data[i];
    }
    
    /** \brief Array-like indexing 
    */
    size_t &operator()(size_t i) {
#if !O2SCL_NO_RANGE_CHECK
      if (i>=size_) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		 +" in permutation::operator(). Size: "+
		 szttos(size_)+
		 " (index should be less than size).").c_str(),exc_eindex);
	return data[0];
      }
#endif
      return data[i];
    }
    
    /** \brief Array-like indexing 
    */
    const size_t &operator()(size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
      if (i>=size_) {
	O2SCL_ERR((((std::string)"Array index ")+szttos(i)+" out of bounds"
		 +" in permutation::operator(). Size: "+
		 szttos(size_)+
		 " (index should be less than size).").c_str(),exc_eindex);
	return data[0];
      }
#endif
      return data[i];
    }
    
    /** \brief Get (with optional range-checking) */
    size_t get(size_t i) const {
#if !O2SCL_NO_RANGE_CHECK
      if (i>=size_) {
	O2SCL_ERR((((std::string)"Permutation index ")+
		   szttos(i)+" out of bounds"+
		   " in permutation::get(). Size: "+
		   szttos(size_)+
		   " (index should be less than size).").c_str(),exc_eindex);
	return 0;
      }
#endif
      return data[i];
    }
    
    /** \brief Set (with optional range-checking) */
    int set(size_t i, size_t val) {
#if !O2SCL_NO_RANGE_CHECK
      if (i>=size_) {
	O2SCL_ERR((((std::string)"Permutation index ")+szttos(i)+
		     " out of bounds"+" in permutation::set(). Size: "+
		     szttos(size_)+
		     " (index should be less than size).").c_str(),
		      exc_eindex);
      }
#endif
      data[i]=val;
      return 0;
    }

    /// Initialize permutation to the identity
    int init() {
      for(size_t i=0;i<size_;i++) data[i]=i;
      return 0;
    }

    /** \brief Return permutation size 
	
	If no memory has been allocated, this will quietly 
	return zero.
    */
    size_t size() const {
      return size_;
    }

    /** \brief Allocate memory for a permutation of size \c dim
    */
    int allocate(size_t dim) {
      if (size_!=dim && size_>0) free();
      size_=dim;
      if (dim>0) {
	data.resize(dim);
      }
      return 0;
    }

    /** \brief Free the memory 
	
	This function will safely do nothing if used without first
	allocating memory or if called multiple times in succession.
    */
    int free() {
      if (size_>0) {
	data.resize(0);
	size_=0;
      }
      return 0;
    }

    /// Resize
    void resize(size_t dim) {
      free();
      allocate(dim);
    }
    //@}

    /// Swap two elements of a permutation
    int swap(const size_t i, const size_t j) {
      size_t tmp=data[i];
      data[i]=data[j];
      data[j]=tmp;
      return 0;
    }

    /// Check to see that a permutation is valid
    bool valid() const {
      for(size_t i=0;i<size_;i++) {
	if (data[i]>size_) return false;
	for(size_t j=0;j<i;j++) {
	  if (data[i]==data[j]) return false;
	}
      }
      return true;
    }

    /// Reverse the permutation
    int reverse() {
      size_t i;
      for (i = 0; i < (size_ / 2); i++){
	size_t j = size_ - i - 1;
	
	size_t tmp = this->data[i] ;
	this->data[i] = this->data[j] ;
	this->data[j] = tmp ;
      }
      return 0;
    }

    /// Compute the inverse of a permutation
    permutation inverse() const {
      permutation p(size_);
      for(size_t i=0;i<size_;i++) {
	p.data[data[i]]=i;
      }
      return p;
    }
    
    /// Apply the permutation to a vector
    template<class vec_t> int apply(vec_t &v) const {
      size_t i, k, pk;
      for(i=0;i<size_;i++) {
	k=data[i];
	while (k>i) k=data[k];
	if (k<i) continue;
	// Now have k==i, i.e. the least in its cycle
	pk=data[k];
	if (pk==i) continue;
	// Shuffle the elements of the cycle
	{
	  double t=v[i];
	  while (pk!=i) {
	    double r1=v[pk];
	    v[k]=r1;
	    k=pk;
	    pk=data[k];
	  }
	  v[k]=t;
	}
      }
      return 0;
    }

    /// Apply the inverse permutation to a vector
    template<class vec_t> int apply_inverse(vec_t &v) const {
      size_t i, k, pk;
      for(i=0;i<size_;i++) {
	k=data[i];
	while (k>i) k=data[k];
	if (k<i) continue;
	// Now have k==i, i.e. the least in its cycle
	pk=data[k];
	if (pk==i) continue;
	// Shuffle the elements of the cycle
	{
	  double t=v[k];
	  while (pk!=i) {
	    double r1=v[pk];
	    v[pk]=t;
	    t=r1;
	    k=pk;
	    pk=data[k];
	  }
	  v[pk]=t;
	}
      }
      return 0;
    }

    // End of permutation class
  };

  /** \brief Output operator for permutations

      A space is output between the permutation elements but no
      space or endline character is output after the last element.

      If the size is zero, this function outputs nothing and does
      not call the error handler.
  */
  std::ostream &operator<<(std::ostream &os, const permutation &p);

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

