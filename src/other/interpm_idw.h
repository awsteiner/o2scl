/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_INTERPM_IDW_H
#define O2SCL_INTERPM_IDW_H

/** \file interpm_idw.h
    \brief File defining \ref o2scl::interpm_idw
*/

#include <iostream>
#include <string>
#include <cmath>

#include <gsl/gsl_combination.h>

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multi-dimensional interpolation by inverse distance
      weighting

      This class performs interpolation on a multi-dimensional data
      set specified as a series of scattered points using the inverse
      distance-weighted average of the nearest three points. The
      function \ref set_data() takes as input: the number of
      dimensions, the number of points which specify the data, and a
      "vector of vectors", e.g. <tt>std::vector<std::vector<double>
      ></tt> which contains the data for all the points.
  */
  template<class vec_t> class interpm_idw {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    
    interpm_idw() {
      data_set=false;
      scales.resize(1);
      scales[0]=1.0;
    }

    /// Distance scales for each coordinate
    ubvector scales;

    /** \brief Initialize the data for the interpolation

	The object \c vecs should be a vector (of size <tt>dim+1</tt>)
	of vectors (all of size <tt>n_points</tt>). It may have any
	type for which the data can be accessed through 
	<tt>operator[][]</tt>. 
     */
    template<class vec_vec_t>
      void set_data(size_t dim, size_t n_points, vec_vec_t &vecs) {

      if (n_points<3) {
	O2SCL_ERR2("Must provide at least three points in ",
		       "interpm_idw::set_data()",exc_efailed);
      }
      np=n_points;
      nd=dim;
      ptrs.resize(dim+1);
      for(size_t i=0;i<dim+1;i++) {
	ptrs[i]=&vecs[i];
      }
      data_set=true;

      return;
    }

    /** \brief Perform the interpolation 
    */
    template<class vec2_t> double operator()(vec2_t &x) const {
      return eval(x);
    }

    /** \brief Perform the interpolation 
    */
    template<class vec2_t> double eval(vec2_t &x) const {
      
      if (data_set==false) {
	O2SCL_ERR("Data not set in interpm_idw::eval_points().",
		  exc_einval);
      }

      // Find the three closest points by
      // exhaustively searching the data
      
      // Put in initial points
      size_t i1=0; 
      size_t i2=1; 
      size_t i3=2;
      double c1=dist(i1,x);
      double c2=dist(i2,x);
      double c3=dist(i3,x);

      // Sort initial points
      if (c2<c1) {
	if (c3<c2) {
	  // 321
	  swap(i1,c1,i3,c3);
	} else if (c3<c1) {
	  // 231
	  swap(i1,c1,i2,c2);
	  swap(i2,c2,i3,c3);
	} else {
	  // 213
	  swap(i1,c1,i2,c2);
	}
      } else {
	if (c3<c1) {
	  // 312
	  swap(i1,c1,i3,c3);
	  swap(i2,c2,i3,c3);
	} else if (c3<c2) {
	  // 132
	  swap(i3,c3,i2,c2);
	}
	// 123
      }

      // Go through remaining points and sort accordingly
      for(size_t j=3;j<np;j++) {
	size_t i4=j;
	double c4=dist(i4,x);
	if (c4<c1) {
	  swap(i4,c4,i3,c3);
	  swap(i3,c3,i2,c2);
	  swap(i2,c2,i1,c1);
	} else if (c4<c2) {
	  swap(i4,c4,i3,c3);
	  swap(i3,c3,i2,c2);
	} else if (c4<c3) {
	  swap(i4,c4,i3,c3);
	}
      }

      // The function values of the three-closest points
      double f1=(*(ptrs[nd]))[i1];
      double f2=(*(ptrs[nd]))[i2];
      double f3=(*(ptrs[nd]))[i3];

      // Check if any of the distances is zero
      if (c1==0.0) {
	return f1;
      } else if (c2==0.0) {
	return f2;
      } else if (c3==0.0) {
	return f3;
      }

      // Return the inverse-distance weighed average
      double norm=1.0/c1+1.0/c2+1.0/c3;
      return (f1/c1+f2/c2+f3/c3)/norm;
    }
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The number of points
    size_t np;
    /// The number of dimensions
    size_t nd;
    /// Desc
    std::vector<vec_t *> ptrs;
    /// True if the data has been specified
    bool data_set;
    
    /// Compute the distance between \c x and the point at index \c index
    template<class vec2_t> double dist(size_t index, vec2_t &x) const {
      double ret=0.0;
      size_t nscales=scales.size();
      for(size_t i=0;i<nd;i++) {
	ret+=pow((x[i]-(*(ptrs[i]))[index])/scales[i%nscales],2.0);
      }
      return sqrt(ret);
    }

    /// Swap points 1 and 2.
    int swap(size_t &index_1, double &dist_1, size_t &index_2, 
	     double &dist_2) const {

      size_t index_temp;
      double dist_temp;
      
      index_temp=index_1; dist_temp=dist_1;
      index_1=index_2; dist_1=dist_2;
      index_2=index_temp; dist_2=dist_temp;
      
      return 0;
    }
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



