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
#ifndef O2SCL_OVECTOR_REV_TLATE_H
#define O2SCL_OVECTOR_REV_TLATE_H

/** \file ovector_rev_tlate.h
    \brief Reversed vector classes
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
#include <o2scl/ovector_tlate.h>
#include <o2scl/array.h>
#include <o2scl/vector.h>

#ifndef DOXYGENP
namespace o2scl {
#endif

  /** \brief Reversed view of a vector

      \note Note that you can't reverse a reversed vector, and this
      is why this class does not have a constructor of the form
      <tt>ovector_reverse(ovector_reverse &v)</tt>.
      
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_reverse_tlate : 
  public ovector_base_tlate<data_t,vparent_t,block_t> {

  public:

    /** \brief Create a vector from \c dat with size \c n and stride \c s 
     */
    ovector_reverse_tlate(ovector_tlate<data_t,vparent_t,block_t> &v) {
      vparent_t::block=v.block;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;
    }

    /** \brief Create a vector from \c dat with size \c n and stride \c s 
     */
    ovector_reverse_tlate(ovector_view_tlate<data_t,vparent_t,block_t> &v) {
      vparent_t::block=v.block;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;
    }

    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    data_t &operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_reverse_tlate::operator[]. Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Array-like indexing 
    */
    const data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_reverse_tlate::operator[]. Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Array-like indexing 
    */
    data_t &operator()(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_reverse_tlate::operator(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Array-like indexing 
    */
    const data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_reverse_tlate::operator(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Get (with optional range-checking) */
    data_t get(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_reverse_tlate::get(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_reverse_tlate::get_ptr(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data;
      }
#endif
      return vparent_t::data+(vparent_t::size-1-i)*vparent_t::stride;
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const data_t *get_const_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_reverse_tlate::get_const_ptr(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (const data_t *)vparent_t::data;
      }
#endif
      return (const data_t *)(vparent_t::data+
			      (vparent_t::size-1-i)*vparent_t::stride);
    }
    
    /** \brief Set (with optional range-checking) */
    int set(size_t i, data_t val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR_RET((((std::string)"Array index ")+itos(i)+" out of bounds"
		       +" in ovector_reverse_tlate::set(). Size: "+
		       itos(vparent_t::size)+
		       " (index should be less than size).").c_str(),
		      gsl_eindex);
      }
#endif
      vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride]=val;
      return 0;
    }

    /** \brief Set all of the value to be the value \c val */
    int set_all(double val) {
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride]=val;
      }
      return 0;
    }

  };

  /** \brief Reversed view of a vector

      \warning At present, reversing a reversed vector does not give
      the original ordering.

      \future I think that maybe in order to ensure that this isn't
      created from an already reversed vector, this class has to be
      separated from the hierarchy altogether. However, I think this
      might break the smart interpolation stuff. 
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_const_reverse_tlate : 
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {
  public:

    /** \brief Create a vector from \c dat with size \c n and stride \c s 
     */
    ovector_const_reverse_tlate
      (const ovector_const_view_tlate<data_t,vparent_t,block_t> &v) {
      vparent_t::block=v.block;
      vparent_t::data=v.data;
      vparent_t::size=v.size();
      vparent_t::stride=v.stride();
      vparent_t::owner=0;
    }

    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing
    */
    const data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_const_reverse_tlate::operator[]. Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }

    /** \brief Array-like indexing
    */
    const data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_const_reverse_tlate::operator(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }

    /** \brief Get (with optional range-checking) */
    data_t get(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_const_reverse_tlate::get(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }

    /** \brief Get pointer (with optional range-checking) */
    const data_t *get_const_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_const_reverse_tlate::get_const_ptr(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (const data_t *)vparent_t::data;
      }
#endif
      return (const data_t *)(vparent_t::data+
			      (vparent_t::size-1-i)*vparent_t::stride);
    }

  };

  /** \brief Reversed view of a subvector

      \warning At present, reversing a reversed vector does not give
      the original ordering.
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_subvector_reverse_tlate : 
  public ovector_base_tlate<data_t,vparent_t,block_t> {
  public:
    /** \brief Create a vector from \c dat with size \c n and stride \c s 
     */
    ovector_subvector_reverse_tlate
      (ovector_base_tlate<data_t,vparent_t,block_t> &v, size_t offset, 
       size_t n) {
      vparent_t::data=0;
      vparent_t::size=0;
      vparent_t::stride=0;
      if (offset+n-1<v.size()) {
	vparent_t::block=0;
	vparent_t::data=v.data+offset*v.stride();
	vparent_t::size=n;
	vparent_t::stride=v.stride();
	vparent_t::owner=0;
      } else {
	O2SCL_ERR("Failed in ovector_sub_reverse() constructor.",gsl_einval);
      }
    }

    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    data_t &operator[](size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_subvector_reverse_tlate::operator[]. Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Array-like indexing 
    */
    const data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_subvector_reverse_tlate::operator[]. Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Array-like indexing 
    */
    data_t &operator()(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_subvector_reverse_tlate::operator(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Array-like indexing 
    */
    const data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_subvector_reverse_tlate::operator(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Get (with optional range-checking) */
    data_t get(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_subvector_reverse_tlate::get(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Get pointer (with optional range-checking) */
    data_t *get_ptr(size_t i) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"
		   +" in ovector_subvector_reverse_tlate::get_ptr(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data;
      }
#endif
      return vparent_t::data+(vparent_t::size-1-i)*vparent_t::stride;
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const data_t *get_const_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds in "
		   +"ovector_subvector_reverse_tlate::get_const_ptr(). Size: "+
		   itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (const data_t *)vparent_t::data;
      }
#endif
      return (const data_t *)(vparent_t::data+
			      (vparent_t::size-1-i)*vparent_t::stride);
    }
    
    /** \brief Set (with optional range-checking) */
    int set(size_t i, data_t val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR_RET((((std::string)"Array index ")+itos(i)+" out of bounds "
		       +" in ovector_subvector_reverse_tlate::"+
		       "get_const_ptr(). Size: "+itos(vparent_t::size)+
		       " (index should be less than size).").c_str(),
		      gsl_eindex);
      }
#endif
      vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride]=val;
      return 0;
    }

    /** \brief Set all of the value to be the value \c val */
    int set_all(double val) {
      for(size_t i=0;i<vparent_t::size;i++) {
	vparent_t::data[i*vparent_t::stride]=val;
      }
      return 0;
    }

  };

  /** \brief Reversed view of a const subvector
      
      \warning At present, reversing a reversed vector does not give
      the original ordering.
  */
  template<class data_t, class vparent_t, class block_t> 
    class ovector_const_subvector_reverse_tlate : 
  public ovector_const_view_tlate<data_t,vparent_t,block_t> {

  public:

    /** \brief Create a vector from \c dat with size \c n and stride \c s 
     */
    ovector_const_subvector_reverse_tlate
      (const ovector_const_view_tlate<data_t,vparent_t,block_t> &v, 
       size_t offset,
       size_t n) {
      if (offset+n-1<v.size()) {
	vparent_t::block=0;
	vparent_t::data=v.data+offset*v.stride();
	vparent_t::size=n;
	vparent_t::stride=v.stride();
	vparent_t::owner=0;
      } else {
	O2SCL_ERR("Subvector failed in ovector_sub_reverse().",gsl_einval);
      }
    }

    /// \name Get and set methods
    //@{
    /** \brief Array-like indexing 
    */
    const data_t &operator[](size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"+
		   " in ovector_const_subvector_reverse_tlate::operator[]."+
		   " Size: "+itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Array-like indexing 
    */
    const data_t &operator()(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"+
		   " in ovector_const_subvector_reverse_tlate::operator()."+
		   " Size: "+itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Get (with optional range-checking) */
    data_t get(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds"+
		   " in ovector_const_subvector_reverse_tlate::get()."+
		   " Size: "+itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return vparent_t::data[0];
      }
#endif
      return vparent_t::data[(vparent_t::size-1-i)*vparent_t::stride];
    }
    
    /** \brief Get pointer (with optional range-checking) */
    const data_t *get_const_ptr(size_t i) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (i>=vparent_t::size) {
	O2SCL_ERR((((std::string)"Array index ")+itos(i)+" out of bounds in "+
		   "ovector_const_subvector_reverse_tlate::get_const_ptr()."+
		   " Size: "+itos(vparent_t::size)+
		   " (index should be less than size).").c_str(),gsl_eindex);
	return (const data_t *)vparent_t::data;
      }
#endif
      return (const data_t *)(vparent_t::data+
			      (vparent_t::size-1-i)*vparent_t::stride);
    }
    
  };

  /// ovector_reverse typedef
  typedef ovector_reverse_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_reverse;
  /// ovector_const_reverse typedef
  typedef ovector_const_reverse_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_const_reverse;
  /// ovector_subvector_reverse typedef
  typedef ovector_subvector_reverse_tlate<double,gsl_vector_norm,gsl_block> 
    ovector_subvector_reverse;
  /// ovector_const_subvector_reverse typedef
  typedef ovector_const_subvector_reverse_tlate<double,gsl_vector_norm,
    gsl_block> ovector_const_subvector_reverse;
  

  /// ovector_int_reverse typedef
  typedef ovector_reverse_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_reverse;
  /// ovector_int_const_reverse typedef
  typedef ovector_const_reverse_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_const_reverse;
  /// ovector_int_subvector_reverse typedef
  typedef ovector_subvector_reverse_tlate<int,gsl_vector_int,gsl_block_int> 
    ovector_int_subvector_reverse;
  /// ovector_int_const_subvector_reverse typedef
  typedef ovector_const_subvector_reverse_tlate<int,gsl_vector_int,
    gsl_block_int> ovector_int_const_subvector_reverse;


#ifndef DOXYGENP
}
#endif

#endif
