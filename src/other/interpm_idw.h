/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_combination.h>

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>
#include <o2scl/vec_stats.h>
#include <o2scl/linear_solver.h>
#include <o2scl/columnify.h>
#include <o2scl/table.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multi-dimensional interpolation by inverse distance
      weighting

      This class is experimental, particularly the evaluation of
      derivatives.

      This class performs interpolation on a multi-dimensional data
      set specified as a series of scattered points using the inverse
      distance-weighted average of nearby points. 

      The function \ref set_data() takes as input: the number of input
      dimensions, the number of output functions, the number of points
      which specify the data, and a "vector of vectors" which contains
      the data for all the points. The vector of vectors must be of a
      type which allows std::swap on individual elements (which are of
      type <tt>vec_t</tt>). The vector of vectors should be structured
      so that the first index selects the input or output function
      and the second index selects the data point. 

      The number of nearby points which are averaged defaults to 3 and
      can be changed by \ref set_points(). To obtain interpolation
      uncertainties, this class finds one additional nearby point and
      returns the standard deviation of the interpolated value for all
      subsets which are missing one nearby point. One-point
      interpolation corresponds to nearest-neighbor interpolation.

      This class requires a distance metric to weight the
      interpolation, and a Euclidean distance is used. By default, the
      length scales in each direction are automatically determined by
      extent of the data (absolute value of max minus min in each
      direction), but the user can specify the length scales manually
      in \ref set_scales() .

      First derivatives can be obtained using \ref derivs_err() , but
      these derivatives are not obtained from the same approximation
      used in the interpolation. That is, the derivatives returned are
      not equal to exact derivatives from the interpolated function
      (as is the case in, e.g., cubic spline interpolation in one
      dimension). This will typically only be particularly noticable
      near discontinuities.

      If the user specifies an array of pointers, the data can be
      changed between calls to the interpolation, but data points
      cannot be added (as set data separately stores the total number
      of data points) without a new call to \ref set_data(). Also, the
      automatically-determined length scales may need to be recomputed
      by calling \ref auto_scale().

      Increasing the value of \c n_extra away from zero allows the
      interpolation to ignore points in the data set which are
      degenerate because they are too close to each other. Points with
      a distance (normalized by the scales) less than \ref min_dist
      are automatically considered degenerate and only the single
      point closest to the requested coordinate is considered.
      Increasing the value of \c n_extra increases the computational
      time required to compute the nearest points which are
      nondegenerate.

      \todo Make verbose output consistent between the various
      eval() functions.

      \future Share code between the various functions
  */
  template<class mat_t=const_matrix_view_table<> >
    class interpm_idw {
    
  protected:

  public:

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    
  interpm_idw() {
    data_set=false;
    scales.resize(1);
    scales[0]=1.0;
    points=3;
    verbose=0;
    n_extra=0;
    min_dist=1.0e-6;
    dist_expo=2.0;
  }

  /** \brief Exponent in computing distance (default 2.0)
   */
  double dist_expo;
  
  /** \brief Verbosity parameter (default 0)
   */
  int verbose;

  /** \brief The number of extra nearest neighbors
      to include to avoid degeneracies (default 0)
  */
  size_t n_extra;

  /** \brief The minimum distance to consider points as
      non-degenerate (default \f$ 10^{-6} \f$ )
  */
  double min_dist;

  /// \name Get and set functions
  //@{
  /** \brief Set the number of closest points to use
      for each interpolation (default 3)
  */
  void set_points(size_t n) {
    if (n==0) {
      O2SCL_ERR("Points cannot be zero in interpm_idw.",
		o2scl::exc_einval);
    }
    points=n;
    return;
  }

  /** \brief Set the scales for the distance metric
	
      All the scales must be positive and non-zero. The size of the
      vector \c (specified in \c n) must be larger than zero.
  */
  template<class vec2_t> void set_scales(size_t n, vec2_t &v) {
    if (n==0) {
      O2SCL_ERR("Scale vector size cannot be zero in interpm_idw.",
		o2scl::exc_einval);
    }
    for(size_t i=0;i<n;i++) {
      if (v[i]<=0.0) {
	O2SCL_ERR("Scale must be positive and non-zero in interpm_idw.",
		  o2scl::exc_einval);
      }
    }
    scales.resize(n);
    o2scl::vector_copy(n,v,scales);
    return;
  }
    
  /** \brief Initialize the data for the interpolation

      The object \c vecs should be a matrix with a
      first index of size <tt>n_in+n_out</tt> and a second 
      index of size <tt>n_points</tt>. It may have be
      any type which allows the use of <tt>operator(,)</tt>
      and <tt>std::swap</tt>.
  */
  void set_data(size_t n_in, size_t n_out, size_t n_points,
		mat_t &dat, bool auto_scale_flag=true) {
    
    if (n_points<points) {
      O2SCL_ERR2("Not enough points provided in ",
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
    std::swap(data,dat);
    data_set=true;

    if (auto_scale_flag) {
      auto_scale();
    }

    return;
  }

  /** \brief Get the data used for interpolation
   */
  void get_data(size_t &n_in, size_t &n_out, size_t &n_points,
		mat_t &dat) {
    n_points=np;
    n_in=nd_in;
    n_out=nd_out;
    std::swap(data,dat);
    data_set=false;
    n_points=0;
    n_in=0;
    n_out=0;
    return;
  }
    
  /** \brief Automatically determine the length scales from the
      data
  */
  void auto_scale() {
    scales.resize(nd_in);
    for(size_t i=0;i<nd_in;i++) {
      double min=data(i,0), max=min;
      for(size_t j=1;j<np;j++) {
	double val=data(i,j);
	if (val>max) max=val;
	if (val<min) min=val;
      }
      scales[i]=max-min;
    }
    return;
  }
    
  /** \brief Initialize the data for the interpolation
      for only one output function

      The object \c vecs should be a vector (of size <tt>n_in+1</tt>)
      of vectors (all of size <tt>n_points</tt>). It may be
      any type which allows the use of <tt>std::swap</tt> for
      each vector in the list. 
  */
  void set_data(size_t n_in, size_t n_points,
		mat_t &dat, bool auto_scale=true) {
    set_data(n_in,1,n_points,dat,auto_scale);
    return;
  }
  //@}

  /// \name Evaluate interpolation
  //@{
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
    std::vector<size_t> >(dists,points+n_extra,index);

    if (n_extra>0) {
      // Remove degenerate points to ensure accurate interpolation
      bool found=true;
      while (found==true) {
	found=false;
	// Find degenerate points and remove them
	for(size_t j=0;j<points+n_extra;j++) {
	  for(size_t k=j;k<points+n_extra;k++) {
	    double dist_jk=dist(j,k);
	    if (index.size()>points && dist_jk<min_dist) {
	      found=true;
	      index.erase(index.begin()+j);
	    }
	  }
	}
      }
    }
      
    // Check if the closest distance is zero
    if (dists[index[0]]<=0.0) {
      return data(nd_in,index[0]);
    }

    // Compute normalization
    double norm=0.0;
    for(size_t i=0;i<points;i++) {
      norm+=1.0/dists[index[i]];
    }

    // Compute the inverse-distance weighted average
    double ret=0.0;
    for(size_t i=0;i<points;i++) {
      ret+=data(nd_in,index[i])/dists[index[i]];
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
    std::vector<size_t> >(dists,points+1+n_extra,index);

    if (n_extra>0) {
      // Remove degenerate points to ensure accurate interpolation
      bool found=true;
      while (found==true) {
	found=false;
	// Find degenerate points and remove them
	for(size_t j=0;j<points+1+n_extra;j++) {
	  for(size_t k=j;k<points+1+n_extra;k++) {
	    double dist_jk=dist(j,k);
	    if (index.size()>points+1 && dist_jk<min_dist) {
	      found=true;
	      index.erase(index.begin()+j);
	    }
	  }
	}
      }
    }
      
    if (dists[index[0]]<=0.0) {

      // If the closest distance is zero, just set the value
      val=data(nd_in,index[0]);
      err=0.0;
      return;

    } else {

      std::vector<double> vals(points+1);

      for(size_t j=0;j<points+1;j++) {

	// Compute normalization
	double norm=0.0;
	for(size_t i=0;i<points+1;i++) {
	  if (i!=j) norm+=1.0/dists[index[i]];
	}
	  
	// Compute the inverse-distance weighted average
	vals[j]=0.0;
	for(size_t i=0;i<points+1;i++) {
	  if (i!=j) {
	    vals[j]+=data(nd_in,index[i])/dists[index[i]];
	  }
	}
	vals[j]/=norm;

      }

      val=vals[points];
      err=o2scl::vector_stddev(vals);

    }

    return;
  }
    
  /** \brief Perform the interpolation over all the functions,
      at point \c x, storing the result in \c y
  */
  template<class vec2_t, class vec3_t>
  void eval(const vec2_t &x, vec3_t &y) const {
      
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_idw::eval().",
		exc_einval);
    }

    if (verbose>0) {
      std::cout << "interpm_idw: input: ";
      for(size_t k=0;k<nd_in;k++) {
	std::cout << x[k] << " ";
      }
      std::cout << std::endl;
    }
      
    // Compute distances
    std::vector<double> dists(np);
    for(size_t i=0;i<np;i++) {
      dists[i]=dist(i,x);
    }

    // Find closest points
    std::vector<size_t> index;
    o2scl::vector_smallest_index<std::vector<double>,double,
    std::vector<size_t> >(dists,points,index);
    if (verbose>0) {
      for(size_t i=0;i<points;i++) {
	std::cout << "interpm_idw: closest point: ";
	for(size_t k=0;k<nd_in;k++) {
	  std::cout << data(k,index[i]) << " ";
	}
	std::cout << std::endl;
      }
    }
      
    if (n_extra>0) {
      // Remove degenerate points to ensure accurate interpolation
      bool found=true;
      while (found==true) {
	found=false;
	// Find degenerate points and remove them
	for(size_t j=0;j<points+n_extra;j++) {
	  for(size_t k=j;k<points+n_extra;k++) {
	    double dist_jk=dist(j,k);
	    if (index.size()>points && dist_jk<min_dist) {
	      found=true;
	      index.erase(index.begin()+j);
	    }
	  }
	}
      }
    }
      
    // Check if the closest distance is zero, if so, just
    // return the value
    if (dists[index[0]]<=0.0) {
      for(size_t i=0;i<nd_out;i++) {
	y[i]=data(nd_in+i,index[0]);
      }
      if (verbose>0) {
	std::cout << "interpm_idw: distance zero. "
	<< "Returning values at index: " << index[0] << std::endl;
	std::cout << "\t";
	o2scl::vector_out(std::cout,nd_out,y,true);
      }
      return;
    }

    // Compute normalization
    double norm=0.0;
    for(size_t i=0;i<points;i++) {
      norm+=1.0/dists[index[i]];
    }
    if (verbose>0) {
      std::cout << "interpm_idw: norm is " << norm << std::endl;
    }

    // Compute the inverse-distance weighted averages
    for(size_t j=0;j<nd_out;j++) {
      y[j]=0.0;
      for(size_t i=0;i<points;i++) {
	if (j==0 && verbose>0) {
	  std::cout << "interpm_idw: Point: ";
	  for(size_t k=0;k<nd_in;k++) {
	    std::cout << data(k,index[i]) << " ";
	  }
	  std::cout << std::endl;
	}
	y[j]+=data(nd_in+j,index[i])/dists[index[i]];
	if (verbose>0) {
	  std::cout << "interpm_idw: j,points,value,1/dist: "
		    << j << " " << i << " "
		    << data(nd_in+j,index[i]) << " "
		    << 1.0/dists[index[i]] << std::endl;
	}
      }
      y[j]/=norm;
      if (verbose>0) {
	std::cout << "interpm_idw: y[" << j << "]: " << y[j]
		  << std::endl;
      }
    }

    return;
  }
    
  /** \brief Perform the interpolation over all the functions
      giving uncertainties and the sorted index vector 

      The vector \c index is automatically resized to
      a size equal to n_points+1+n_extra
      be larger than 
  */
  template<class vec2_t, class vec3_t, class vec4_t>
  void eval_err_index(const vec2_t &x, vec3_t &val, vec4_t &err,
		      std::vector<size_t> &index) const {
      
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_idw::eval_err().",
		exc_einval);
    }
      
    // Compute distances
    std::vector<double> dists(np);
    for(size_t i=0;i<np;i++) {
      dists[i]=dist(i,x);
    }

    // Find closest points, note that index is automatically resized
    // by the vector_smallest_index function
    o2scl::vector_smallest_index<std::vector<double>,double,
    std::vector<size_t> >(dists,points+1+n_extra,index);

    if (n_extra>0) {
      // Remove degenerate points to ensure accurate interpolation
      bool found=true;
      while (found==true) {
	found=false;
	// Find degenerate points and remove them
	for(size_t j=0;j<points+1+n_extra;j++) {
	  for(size_t k=j;k<points+1+n_extra;k++) {
	    double dist_jk=dist(j,k);
	    if (index.size()>points+1 && dist_jk<min_dist) {
	      found=true;
	      index.erase(index.begin()+j);
	    }
	  }
	}
      }
    }

    if (dists[index[0]]<=0.0) {

      // If the closest distance is zero, just set the values and
      // errors
      for(size_t k=0;k<nd_out;k++) {
	val[k]=data(nd_in+k,index[0]);
	err[k]=0.0;
      }
      return;

    } else {
	
      for(size_t k=0;k<nd_out;k++) {
	  
	std::vector<double> vals(points+1);
	  
	for(size_t j=0;j<points+1;j++) {
	    
	  // Compute normalization
	  double norm=0.0;
	  for(size_t i=0;i<points+1;i++) {
	    if (i!=j) norm+=1.0/dists[index[i]];
	  }
	    
	  // Compute the inverse-distance weighted average
	  vals[j]=0.0;
	  for(size_t i=0;i<points+1;i++) {
	    if (i!=j) {
	      vals[j]+=data(nd_in+k,index[i])/dists[index[i]];
	    }
	  }
	  vals[j]/=norm;
	    
	}

	// Instead of using the average, we report the value as the
	// last element in the array, which is the interpolated
	// value from the closest points
	val[k]=vals[points];
	  
	err[k]=o2scl::vector_stddev(vals);
	  
      }

    }
      
    return;
  }
    
  /** \brief Perform the interpolation over all the functions
      with uncertainties
  */
  template<class vec2_t, class vec3_t, class vec4_t>
  void eval_err(const vec2_t &x, vec3_t &val, vec4_t &err) const {
    std::vector<size_t> index;
    return eval_err_index(x,val,err,index);
  }
  //@}

  /// \name Evaluate derivatives
  //@{
  /** \brief For one of the functions, compute the partial
      derivatives (and uncertainties) with respect to all of the
      inputs at one data point

      \note This function ignores the points chosen by \ref
      set_points() and always chooses to average derivative
      calculations determined from \c n_in+1 combinations of \c n_in
      points .

      \todo Use the mechanism provided by <tt>n_extra</tt> above
      to remove degenerate points. 

      \future This function requires an extra copy from
      "ders" to "ders2" which could be removed.
  */
  template<class vec3_t>
  void derivs_err(size_t func_index, size_t point_index, 
		  vec3_t &derivs, vec3_t &errs) const {
      
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_idw::derivs_err().",
		exc_einval);
    }
      
    // Set x equal to the specified point
    ubvector x(nd_in);
    for(size_t i=0;i<nd_in;i++) {
      x[i]=data(i,point_index);
    }
    // Set f equal to the value of the function at the specified point
    double f=data(nd_in+func_index,point_index);

    // The linear solver
    o2scl_linalg::linear_solver_HH<> lshh;
    
    // Compute distances
    std::vector<double> dists(np);
    for(size_t i=0;i<np;i++) {
      dists[i]=dist(i,x);
    }
  
    // Find closest (but not identical) points

    std::vector<size_t> index;
    size_t max_smallest=(nd_in+2)*2;
    if (max_smallest>np) max_smallest=np;
    if (max_smallest<nd_in+1) {
      O2SCL_ERR("Couldn't find enough nearby points.",o2scl::exc_einval);
    }

    if (verbose>0) {
      std::cout << "max_smallest: " << max_smallest << std::endl;
    }
      
    o2scl::vector_smallest_index<std::vector<double>,double,
    std::vector<size_t> >(dists,max_smallest,index);

    if (verbose>0) {
      for(size_t i=0;i<index.size();i++) {
	std::cout << "index[" << i << "] = " << index[i] << " "
		  << dists[index[i]] << std::endl;
      }
    }
      
    std::vector<size_t> index2;
    for(size_t i=0;i<max_smallest;i++) {
      if (dists[index[i]]>0.0) {
	index2.push_back(index[i]);
	if (index2.size()==nd_in+1) i=max_smallest;
      }
    }
    if (index2.size()<nd_in+1) {
      O2SCL_ERR("Couldn't find enough nearby points (2).",
		o2scl::exc_einval);
    }

    if (verbose>0) {
      for(size_t i=0;i<index2.size();i++) {
	std::cout << "index2[" << i << "] = " << index2[i] << " "
	<< dists[index2[i]] << std::endl;
      }
    }
      
    // Unit vector storage
    std::vector<ubvector> units(nd_in+1);
    // Difference vector norms
    std::vector<double> diff_norms(nd_in+1);
    // Storage for the derivative estimates
    std::vector<ubvector> ders(nd_in+1);
    // Matrix of unit vectors
    ubmatrix m(nd_in,nd_in);
    // Vector of function value differences
    ubvector v(nd_in);
    // Rearranged derivative object
    std::vector<ubvector> ders2(nd_in);
      
    for(size_t i=0;i<nd_in+1;i++) {

      // Assign unit vector elements
      units[i].resize(nd_in);
      for(size_t j=0;j<nd_in;j++) {
	units[i][j]=data(j,index2[i])-x[j];
      }

      // Normalize the unit vectors
      diff_norms[i]=o2scl::vector_norm<ubvector,double>(units[i]);
      for(size_t j=0;j<nd_in;j++) {
	units[i][j]/=diff_norms[i];
      }

    }

    // Verbose output of the closest points and their norms
    if (verbose>0) {
      std::cout << "Point:     ";
      for(size_t i=0;i<nd_in;i++) {
	std::cout << x[i] << " ";
      }
      std::cout << f << std::endl;
      for(size_t j=0;j<nd_in+1;j++) {
	std::cout << "Closest: " << j << " " << index2[j] << " ";
	for(size_t i=0;i<nd_in;i++) {
	  std::cout << data(i,index2[j]) << " ";
	}
	std::cout << data(nd_in+func_index,index2[j]) << " "
		  << diff_norms[j] << std::endl;
      }
      for(size_t j=0;j<nd_in+1;j++) {
	std::cout << "diff_norm: " << j << " " << diff_norms[j]
	<< std::endl;
      }
      // End of verbose output
    }
    
    // Go through each set of points
    for(size_t i=0;i<nd_in+1;i++) {

      ders[i].resize(nd_in);

      // Construct the matrix and vector for the solver
      size_t jj=0;
      for(size_t j=0;j<nd_in+1;j++) {
	if (j!=i) {
	  for(size_t k=0;k<nd_in;k++) {
	    m(jj,k)=units[j][k];
	  }
	  v[jj]=(data(nd_in+func_index,index2[j])-f)/diff_norms[j];
	  jj++;
	}
      }

      // Solve to compute the derivatives
      if (verbose>0) {
	std::cout << "m:" << std::endl;
	o2scl::matrix_out(std::cout,nd_in,nd_in,m);
	std::cout << "v:" << std::endl;
	o2scl::vector_out(std::cout,nd_in,v,true);
      }
      lshh.solve(nd_in,m,v,ders[i]);
      if (verbose>0) {
	std::cout << "Derivs:  " << i << " ";
	std::cout.setf(std::ios::showpos);
	for(size_t j=0;j<nd_in;j++) {
	  std::cout << ders[i][j] << " ";
	}
	std::cout.unsetf(std::ios::showpos);
	std::cout << std::endl;
      }

      // Go to next derivative estimate
    }
      
    for(size_t i=0;i<nd_in;i++) {

      // Rearrange derivatives
      ders2[i].resize(nd_in+1);
      for(size_t j=0;j<nd_in+1;j++) {
	ders2[i][j]=ders[j][i];
      }

      // Compute mean and standard deviation
      derivs[i]=o2scl::vector_mean(ders2[i]);
      errs[i]=o2scl::vector_stddev(ders2[i]);
    }
      
    return;
  }
  //@}
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
  /// The number of points
  size_t np;
  /// The number of dimensions of the inputs
  size_t nd_in;
  /// The number of dimensions of the outputs
  size_t nd_out;
  /// The copy of the data
  mat_t data;
  /// True if the data has been specified
  bool data_set;
  /// Number of points to include in each interpolation (default 3)
  size_t points;

  /// \name Distance determination [protected]
  //@{
  /// Distance scales for each coordinate
  ubvector scales;

  /** \brief Compute the distance between \c x and the point at
      index \c index
  */
  template<class vec2_t> double dist(size_t index,
				     const vec2_t &x) const {
    double ret=0.0;
    size_t nscales=scales.size();
    for(size_t i=0;i<nd_in;i++) {
      ret+=pow((x[i]-data(i,index))/scales[i%nscales],dist_expo);
    }
    return sqrt(ret);
  }

  /** \brief Compute the distance between two points in the
      data set
  */
  double dist(size_t j, size_t k) const {
    double ret=0.0;
    size_t nscales=scales.size();
    for(size_t i=0;i<nd_in;i++) {
      ret+=pow((data(i,j)-data(i,k))/scales[i%nscales],dist_expo);
    }
    return sqrt(ret);
  }
  //@}
    
#endif
    
  };
    
#ifndef DOXYGEN_NO_O2NS
}
#endif
    
#endif

