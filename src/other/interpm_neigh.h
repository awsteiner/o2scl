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
#ifndef O2SCL_INTERPM_NEIGH_H
#define O2SCL_INTERPM_NEIGH_H

#include <iostream>
#include <string>
#include <cmath>

#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Nearest-neighbor interpolation in arbitrary number of 
      dimensions
  */
  template<class vec_t> class interpm_neigh {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;

    interpm_neigh() {
      data_set=false;
      scales.resize(1);
      scales[0]=1.0;
    }

    /// Distance scales for each coordinate
    ubvector scales;

    /** \brief Initialize the data for the neigh interpolation
     */
    template<class vec_vec_t>
      void set_data(size_t dim, size_t n_points, vec_vec_t &vecs) {

      if (n_points<1) {
	O2SCL_ERR2("Must provide at least one point in ",
		   "interpm_neigh::set_data()",exc_efailed);
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
    template<class vec2_t> double eval(vec2_t &x) const {
      
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp_planar::eval_points().",
		  exc_einval);
      }
      
      // Check that some scales have been given
      size_t nscales=scales.size();
      if (nscales==0) {
	O2SCL_ERR("No scales in interpm_neigh::eval().",exc_einval);
      }

      // Exhaustively search the data
      size_t index_min=0;
      double dist_min=0.0;
      for(size_t i=0;i<nd;i++) {
	dist_min+=pow((x[i]-(*(ptrs[i]))[index_min])/scales[i%nscales],2.0);
      }
      for(size_t index=1;index<np;index++) {
	double dist=0.0;
	for(size_t i=0;i<nd;i++) {
	  dist+=pow((x[i]-(*(ptrs[i]))[index])/scales[i%nscales],2.0);
	}
	if (dist<dist_min) {
	  dist_min=dist;
	  index_min=index;
	}
      }
      
      // Return the function value
      return (*(ptrs[nd]))[index_min];
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
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



