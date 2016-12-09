/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#include <o2scl/vec_stats.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multi-dimensional interpolation by inverse distance
      weighting

      This class performs interpolation on a multi-dimensional data
      set specified as a series of scattered points using the inverse
      distance-weighted average of nearby points. The
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
      order=3;
    }

    /** \brief Set the number of closest points to use
	for each interpolation (default 3)
    */
    void set_order(size_t n) {
      if (n==0) {
	O2SCL_ERR("Order cannot be zero in interpm_idw.",
		  o2scl::exc_einval);
      }
      order=n;
      return;
    }

    /** \brief Set the scales for the distance metric
     */
    template<class vec2_t> void set_scales(size_t n, vec2_t &v) {
      if (n==0) {
	O2SCL_ERR("Scale vector size cannot be zero in interpm_idw.",
		  o2scl::exc_einval);
      }
      scales.resize(n);
      o2scl::vector_copy(n,v,scales);
      return;
    }
    
    /** \brief Initialize the data for the interpolation

	The object \c vecs should be a vector (of size <tt>n_in+n_out</tt>)
	of vectors (all of size <tt>n_points</tt>). It may have be
	any time which allows the use of <tt>std::swap</tt> for
	each vector in the list. 
    */
    template<class vec_vec_t>
      void set_data(size_t n_in, size_t n_out, size_t n_points,
		    vec_vec_t &vecs, bool auto_scale=true) {

      if (n_points<3) {
	O2SCL_ERR2("Must provide at least three points in ",
		   "interpm_idw::set_data()",exc_efailed);
      }
      if (n_in<1) {
	O2SCL_ERR2("Must provide at least one input column in ",
		   "interpm_idw::set_data()",exc_efailed);
      }
      if (n_out<1) {
	O2SCL_ERR2("Must provide at least one output column in ",
		   "interpm_idw::set_data()",exc_efailed);
      }
      np=n_points;
      nd_in=n_in;
      nd_out=n_out;
      ptrs.resize(n_in+n_out);
      for(size_t i=0;i<n_in+n_out;i++) {
	std::swap(ptrs[i],vecs[i]);
      }
      data_set=true;

      if (auto_scale) {
	scales.resize(n_in);
	for(size_t i=0;i<n_in;i++) {
	  scales[i]=fabs(o2scl::vector_max_value<vec_t,double>(ptrs[i])-
			 o2scl::vector_min_value<vec_t,double>(ptrs[i]));
	}
      }

      return;
    }

    /** \brief Initialize the data for the interpolation

	The object \c vecs should be a vector (of size <tt>n_in+1</tt>)
	of vectors (all of size <tt>n_points</tt>). It may have be
	any time which allows the use of <tt>std::swap</tt> for
	each vector in the list. 
    */
    template<class vec_vec_t>
      void set_data(size_t n_in, size_t n_points,
		    vec_vec_t &vecs, bool auto_scale=true) {
      set_data(n_in,1,n_points,vecs,auto_scale);
      return;
    }

    /** \brief Perform the interpolation over the first function
     */
    template<class vec2_t> double operator()(const vec2_t &x) const {
      return eval(x);
    }

    /** \brief Perform the interpolation over the first function
     */
    template<class vec2_t> double eval(const vec2_t &x) const {
    
      if (data_set==false) {
	O2SCL_ERR("Data not set in interpm_idw::eval().",
		  exc_einval);
      }
    
      // Compute distances
      std::vector<double> dists(np);
      for(size_t i=0;i<np;i++) {
	dists[i]=dist(i,x);
      }

      // Find closest points
      std::vector<size_t> index;
      o2scl::vector_smallest_index<std::vector<double>,double,
	std::vector<size_t> >(dists,order,index);

      // Check if the closest distance is zero
      if (dists[index[0]]<=0.0) {
	return ptrs[nd_in][index[0]];
      }

      // Compute normalization
      double norm=0.0;
      for(size_t i=0;i<order;i++) {
	norm+=1.0/dists[index[i]];
      }

      // Compute the inverse-distance weighted average
      double ret=0.0;
      for(size_t i=0;i<order;i++) {
	ret+=ptrs[nd_in][index[i]]/dists[index[i]];
      }
      ret/=norm;

      // Return the average
      return ret;
    }
    
    /** \brief Perform the interpolation over the first function
	with uncertainty
     */
    template<class vec2_t> void eval_err(const vec2_t &x, double &val,
					 double &err) const {
      
      if (data_set==false) {
	O2SCL_ERR("Data not set in interpm_idw::eval_err().",
		  exc_einval);
      }
      
      // Compute distances
      std::vector<double> dists(np);
      for(size_t i=0;i<np;i++) {
	dists[i]=dist(i,x);
      }

      // Find closest points
      std::vector<size_t> index;
      o2scl::vector_smallest_index<std::vector<double>,double,
	std::vector<size_t> >(dists,order+1,index);

      if (dists[index[0]]<=0.0) {

	// If the closest distance is zero, just set the value
	val=ptrs[nd_in][index[0]];
	err=0.0;
	return;

      } else {

	std::vector<double> vals(order+1);

	for(size_t j=0;j<order+1;j++) {

	  // Compute normalization
	  double norm=0.0;
	  for(size_t i=0;i<order+1;i++) {
	    if (i!=j) norm+=1.0/dists[index[i]];
	  }
	  
	  // Compute the inverse-distance weighted average
	  vals[j]=0.0;
	  for(size_t i=0;i<order+1;i++) {
	    if (i!=j) {
	      vals[j]+=ptrs[nd_in][index[i]]/dists[index[i]];
	    }
	  }
	  vals[j]/=norm;

	}

	val=vals[order];
	err=o2scl::vector_stddev(vals);

      }

      return;
    }
    
    /** \brief Perform the interpolation over all the functions,
	storing the result in \c y
    */
    template<class vec2_t, class vec3_t>
      void eval(vec2_t &x, vec3_t &y) const {
      
      if (data_set==false) {
	O2SCL_ERR("Data not set in interpm_idw::eval().",
		  exc_einval);
      }

      // Compute distances
      std::vector<double> dists(np);
      for(size_t i=0;i<np;i++) {
	dists[i]=dist(i,x);
      }

      // Find closest points
      std::vector<size_t> index;
      o2scl::vector_smallest_index<std::vector<double>,double,
	std::vector<size_t> >(dists,order,index);

      // Check if the closest distance is zero
      if (dists[index[0]]<=0.0) {
	for(size_t i=0;i<nd_out;i++) {
	  y[i]=ptrs[index[0]][nd_in+i];
	}
      }

      // Compute normalization
      double norm=0.0;
      for(size_t i=0;i<order;i++) {
	norm+=1.0/dists[index[i]];
      }

      // Compute the inverse-distance weighted averages
      for(size_t j=0;j<nd_out;j++) {
	y[j]=0.0;
	for(size_t i=0;i<order;i++) {
	  y[j]+=ptrs[nd_in][index[i]]/dists[index[i]];
	}
	y[j]/=norm;
      }

      return;
    }
    
    /** \brief Perform the interpolation over all the functions
	with uncertainties
     */
    template<class vec2_t, class vec3_t>
      void eval_err(const vec2_t &x, vec3_t &val, vec3_t &err) const {
      
      if (data_set==false) {
	O2SCL_ERR("Data not set in interpm_idw::eval_err().",
		  exc_einval);
      }
      
      // Compute distances
      std::vector<double> dists(np);
      for(size_t i=0;i<np;i++) {
	dists[i]=dist(i,x);
      }

      // Find closest points
      std::vector<size_t> index;
      o2scl::vector_smallest_index<std::vector<double>,double,
	std::vector<size_t> >(dists,order+1,index);

      if (dists[index[0]]<=0.0) {

	// If the closest distance is zero, just set the values and
	// errors
	for(size_t k=0;k<nd_out;k++) {
	  val[k]=ptrs[nd_in+k][index[0]];
	  err[k]=0.0;
	}
	return;

      } else {
	
	for(size_t k=0;k<nd_out;k++) {
	  
	  std::vector<double> vals(order+1);
	  
	  for(size_t j=0;j<order+1;j++) {
	    
	    // Compute normalization
	    double norm=0.0;
	    for(size_t i=0;i<order+1;i++) {
	      if (i!=j) norm+=1.0/dists[index[i]];
	    }
	    
	    // Compute the inverse-distance weighted average
	    vals[j]=0.0;
	    for(size_t i=0;i<order+1;i++) {
	      if (i!=j) {
		vals[j]+=ptrs[nd_in+k][index[i]]/dists[index[i]];
	      }
	    }
	    vals[j]/=norm;
	    
	  }
	  
	  val[k]=vals[order];
	  err[k]=o2scl::vector_stddev(vals);
	  
	}

      }

      return;
    }

#ifdef O2SCL_NEVER_DEFINED
    
    template<class vec2_t, class vec3_t>
      void derivs_err(const vec2_t &x, size_t ix,
		      vec3_t &derivs, vec3_t &errs) const {
      
      if (data_set==false) {
	O2SCL_ERR("Data not set in interpm_idw::eval_err().",
		  exc_einval);
      }

      o2scl_linalg::linear_solver_HH<> lshh;
      
      // Compute distances
      std::vector<double> dists(np);
      for(size_t i=0;i<np;i++) {
	dists[i]=dist(i,x);
      }

      // Find nd_in+1 closest points
      std::vector<size_t> index;
      o2scl::vector_smallest_index<std::vector<double>,double,
	std::vector<size_t> >(dists,nd_in+1,index);

      if (dist[0]<=0.0) {
	O2SCL_ERR("Derivative algorithm fails if a distance is zero.",
		  o2scl::exc_einval);
      }

      // Unit vector storage
      std::vector<ubvector> units(nd_in+1);
      // Difference vector norms
      std::vector<double> diff_norms(nd_in+1);

      for(size_t i=0;i<nd_in+1;i++) {

	// Assign unit vector elements
	units[i].resize(nd_in);
	for(size_t j=0;j<nd_in;j++) {
	  units[i][j]=ptrs[j][index[i]]-x[j];
	}

	// Normalize the unit vectors
	diff_norms[i]=o2scl::vector_norm(units[i]);
	for(size_t j=0;j<nd_in;j++) {
	  units[i][j]/=diff_norms[i];
	}

      }

      std::vector<ubvector> ders(nd_in+1);
      
      // Go through each set of points
      for(size_t i=0;i<nd_in+1;i++) {

	ders[i].resize(nd_in);

	// Construct the matrix and vector for the solver
	ubmatrix m(nd_in,nd_in);
	ubvector v(nd_in);
	size_t jj=0;
	for(size_t j=0;j<nd_in+1;j++) {
	  if (j!=i) {
	    for(size_t k=0;k<nd_in;k++) {
	      m(jj,k)=units[j][k];
	    }
	    v[jj]=(ptrs[ix+nd_in][index[i]]-f[i])/diff_norms[jj];
	  }
	  jj++;
	}

	// Solve to compute the derivatives
	lshh.solve(nd_in,m,v,ders[i]);
	
      }

      // Rearranged derivative object
      std::vector<ubvector> ders2(nd_in);
      
      for(size_t i=0;i<nd_in;i++) {

	// Rearrange derivatives
	ders2[i].resize(nd_in+1);
	for(size_t j=0;j<nd_in+1;j++) {
	  ders2[i][j]=ders[j][i];
	}

	// Compute mean and standard deviation
	derivs[i]=o2scl::vector_mean(ders[i]);
	errs[i]=o2scl::vector_stddev(ders[i]);
      }
      
      return;
    }
    
#endif
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// Distance scales for each coordinate
    ubvector scales;

    /// The number of points
    size_t np;
    /// The number of dimensions of the inputs
    size_t nd_in;
    /// The number of dimensions of the outputs
    size_t nd_out;
    /// A vector of pointers holding the data
    std::vector<vec_t> ptrs;
    /// True if the data has been specified
    bool data_set;
    /// Number of points to include in each interpolation (default 3)
    size_t order;
    
    /// Compute the distance between \c x and the point at index \c index
    template<class vec2_t> double dist(size_t index, const vec2_t &x) const {
      double ret=0.0;
      size_t nscales=scales.size();
      for(size_t i=0;i<nd_in;i++) {
	ret+=pow((x[i]-ptrs[i][index])/scales[i%nscales],2.0);
      }
      return sqrt(ret);
    }

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



