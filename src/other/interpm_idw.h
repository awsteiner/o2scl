/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
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
#include <o2scl/tensor.h>
#include <o2scl/interpm_base.h>

namespace o2scl {

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

      \verbatim embed:rst

      .. todo:

         In class interpm_idw:

         - Document the algorithm used in derivs_err.
         - Handle error cases and non-zero return values better.
         - Make verbose output consistent between the various
           eval() functions.
         - Future: Share code between the various functions.
         
      \endverbatim
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_x_t=o2scl::const_matrix_view_table<>,
           class mat_y_t=o2scl::const_matrix_view_table<> >
  class interpm_idw :
    public interpm_base<vec_t,mat_x_t,mat_y_t> {
    
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
      n_extra=0;
      min_dist=1.0e-6;
      dist_expo=2.0;
      rescale=true;
    }

    /** \brief Exponent in computing distance (default 2.0)
     */
    double dist_expo;
  
    /** \brief The number of extra nearest neighbors
        to include to avoid degeneracies (default 0)
    */
    size_t n_extra;

    /** \brief The minimum distance to consider points as
        non-degenerate (default \f$ 10^{-6} \f$ )
    */
    double min_dist;

    /** \brief If true, then the data will be automatically rescaled
        (default true)
    */
    bool rescale;
    
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
    virtual int set_data(size_t n_in, size_t n_out, size_t n_pts,
                         mat_x_t &dat_x, mat_y_t &dat_y) {
    
      if (n_pts<points) {
        O2SCL_ERR2("Not enough pts provided in ",
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
      this->n_points=n_pts;
      this->n_params=n_in;
      this->n_outputs=n_out;
      std::swap(data_x,dat_x);
      std::swap(data_y,dat_y);
      data_set=true;

      if (rescale) {
        auto_scale();
      }

      return 0;
    }

    /** \brief Get the data used for interpolation
     */
    void get_data(size_t &n_in, size_t &n_out, size_t &n_pts,
                  mat_x_t &dat_x, mat_y_t &dat_y) {
      n_pts=this->n_points;
      n_in=this->n_params;
      n_out=this->n_outputs;
      std::swap(data_x,dat_x);
      std::swap(data_y,dat_y);
      data_set=false;
      n_pts=0;
      n_in=0;
      n_out=0;
      return;
    }
    
    /** \brief Automatically determine the length scales from the
        data
    */
    void auto_scale() {
      scales.resize(this->n_params);
      for(size_t i=0;i<this->n_params;i++) {
        double min=data_x(0,i), max=min;
        for(size_t j=1;j<this->n_points;j++) {
          double val=data_x(j,i);
          if (val>max) max=val;
          if (val<min) min=val;
        }
        scales[i]=max-min;
	if (scales[i]==0.0) scales[i]=1.0;
      }
      return;
    }
    
    /** \brief Initialize the data for the interpolation
        for only one output function

        The object \c vecs should be a vector (of size <tt>n_in+1</tt>)
        of vectors (all of size <tt>n_pts</tt>). It may be
        any type which allows the use of <tt>std::swap</tt> for
        each vector in the list. 
    */
    void set_data(size_t n_in, size_t n_pts,
                  mat_x_t &dat_x, mat_y_t &dat_y) {
      set_data(n_in,1,n_pts,dat_x,dat_y);
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
    template<class vec2_t> double eval_one_tl(const vec2_t &x) const {
    
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_idw::eval().",
                  exc_einval);
      }
    
      // Compute distances
      std::vector<double> dists(this->n_points);
      for(size_t i=0;i<this->n_points;i++) {
        dists[i]=dist(i,x);
      }

      // Find closest points
      std::vector<size_t> index;
      o2scl::vector_smallest_index<std::vector<double>,double,
                                   std::vector<size_t> >
        (dists,points+n_extra,index);
      
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
        return data_y(index[0],0);
      }

      // Compute normalization
      double norm=0.0;
      for(size_t i=0;i<points;i++) {
        norm+=1.0/dists[index[i]];
      }

      // Compute the inverse-distance weighted average
      double ret=0.0;
      for(size_t i=0;i<points;i++) {
        ret+=data_y(index[i],0)/dists[index[i]];
      }
      ret/=norm;

      // Return the average
      return ret;
    }
    
    /** \brief Desc
     */
    virtual int eval(const vec_t &x, vec_t &y) const {
      return eval_tl(x,y);
    }
    
    /** \brief Perform the interpolation over the first function
        with uncertainty
    */
    template<class vec2_t> void eval_one_unc_tl
    (const vec2_t &x, double &val, double &err, double &extrap) {
      
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_idw::eval_err().",
                  exc_einval);
      }

      // Compute distances
      std::vector<double> dists(this->n_points);
      for(size_t i=0;i<this->n_points;i++) {
        dists[i]=dist(i,x);
      }

      // Find closest points
      std::vector<size_t> index;
      o2scl::vector_smallest_index<std::vector<double>,double,
                                   std::vector<size_t> >
        (dists,points+1+n_extra,index);

      if (this->verbose>1) {
        std::cout << "interpm_idw::eval_err(): n_extra is " << n_extra
                  << std::endl;
      }
      
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
                if (this->verbose>1) {
                  std::cout << "  Found degenerate point." << std::endl;
                }
              }
            }
          }
        }
      }
      
      if (dists[index[0]]<=0.0) {

        // If the closest distance is zero, just set the value
        if (this->verbose>1) {
          std::cout << "  Closest distance is zero." << std::endl;
        }
        val=data_y(index[0],0);
        err=0.0;
        extrap=0.0;
        return;

      } else {

        std::vector<double> vals(points+1);
        /// Distances between the selected points and the first point
        std::vector<double> dists_bw;
        /// Distances from the selected point to the user point
        std::vector<double> dists_to;

        // We construct points+1 estimates of the result and take the
        // average. Use the standard deviation for the uncertainty.
        for(size_t j=0;j<points+1;j++) {

          // Compute normalization
          double norm=0.0;
          for(size_t i=0;i<points+1;i++) {
            if (i!=j) norm+=1.0/dists[index[i]];
            dists_to.push_back(dists[index[i]]);
            if (i!=0) {
              double d=dist(index[0],index[i]);
              if (d!=0.0) {
                dists_bw.push_back(d);
              }
            }
          }
	  
          // Compute the inverse-distance weighted average
          vals[j]=0.0;
          for(size_t i=0;i<points+1;i++) {
            if (i!=j) {
              vals[j]+=data_y(index[i],0)/dists[index[i]];
            }
          }

          vals[j]/=norm;

        }

        val=vals[points];
        err=o2scl::vector_stddev(vals);
        if (dists_bw.size()>0) {
          double max_bw=vector_max_value<std::vector<double>,double>
	    (dists_bw);
          double min_to=vector_min_value<std::vector<double>,double>
	    (dists_to);
          extrap=min_to/max_bw;
        } else {
          extrap=0.0;
        }
        
      }

      return;
    }
    
    /** \brief Perform the interpolation over all the functions,
        at point \c x, storing the result in \c y
    */
    template<class vec2_t, class vec3_t>
    int eval_tl(const vec2_t &x, vec3_t &y) const {
      
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_idw::eval().",
                  exc_einval);
      }

      if (this->verbose>0) {
        std::cout << "interpm_idw: input: ";
        for(size_t k=0;k<this->n_params;k++) {
          std::cout << x[k] << " ";
        }
        std::cout << std::endl;
      }
      
      // Compute distances
      std::vector<double> dists(this->n_points);
      for(size_t i=0;i<this->n_points;i++) {
        dists[i]=dist(i,x);
      }

      // Find closest points
      std::vector<size_t> index;
      o2scl::vector_smallest_index<std::vector<double>,double,
                                   std::vector<size_t> >
        (dists,points,index);
      if (this->verbose>0) {
        for(size_t i=0;i<points;i++) {
          std::cout << "interpm_idw: closest point: ";
          for(size_t k=0;k<this->n_params;k++) {
            std::cout << data_x(index[i],k) << " ";
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
        for(size_t i=0;i<this->n_outputs;i++) {
          y[i]=data_y(index[0],i);
        }
        if (this->verbose>0) {
          std::cout << "interpm_idw: distance zero. "
                    << "Returning values at index: " << index[0]
		    << std::endl;
          std::cout << "\t";
          o2scl::vector_out(std::cout,this->n_outputs,y,true);
        }
        return 0;
      }

      // Compute normalization
      double norm=0.0;
      for(size_t i=0;i<points;i++) {
        norm+=1.0/dists[index[i]];
      }
      if (this->verbose>0) {
        std::cout << "interpm_idw: norm is " << norm << std::endl;
      }

      // Compute the inverse-distance weighted averages
      for(size_t j=0;j<this->n_outputs;j++) {
        y[j]=0.0;
        for(size_t i=0;i<points;i++) {
          if (j==0 && this->verbose>0) {
            std::cout << "interpm_idw: Point: ";
            for(size_t k=0;k<this->n_params;k++) {
              std::cout << data_x(index[i],k) << " ";
            }
            std::cout << std::endl;
          }
          y[j]+=data_y(index[i],j)/dists[index[i]];
          if (this->verbose>0) {
            std::cout << "interpm_idw: j,points,value,1/dist: "
                      << j << " " << i << " "
                      << data_y(index[i],j) << " "
                      << 1.0/dists[index[i]] << std::endl;
          }
        }
        y[j]/=norm;
        if (this->verbose>0) {
          std::cout << "interpm_idw: y[" << j << "]: " << y[j]
                    << std::endl;
        }
      }

      return 0;
    }

    /** \brief Desc
     */
    virtual int eval_unc(const vec_t &x, vec_t &y, vec_t &y_unc) const {
      std::vector<size_t> ix;
      std::vector<double> extrap;
      return eval_unc_tl_index(x,y,y_unc,ix,extrap);
    }
    
    /** \brief Perform the interpolation over all the functions
        giving uncertainties and the sorted index vector 

        The vector \c index is automatically resized to a size equal to
        n_points+1+n_extra.
    */
    template<class vec2_t, class vec3_t, class vec4_t>
    int eval_unc_tl_index(const vec2_t &x, vec3_t &val, vec4_t &err,
                          std::vector<size_t> &index,
                          std::vector<double> &extrap) const {
      
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_idw::eval_err().",
                  exc_einval);
      }
      
      if (this->verbose>1) {
        std::cout << "interpm_idw::eval_err_index(): n_extra: " << n_extra
                  << std::endl;
      }

      if (this->verbose>2) {
	std::cout << "x: ";
	o2scl::vector_out(std::cout,x,true);
      }

      extrap.resize(this->n_outputs);
      
      // Compute distances
      std::vector<double> dists(this->n_points);
      for(size_t i=0;i<this->n_points;i++) {
        dists[i]=dist(i,x);
	if (!std::isfinite(dists[i])) {
	  std::cout << "i,dists[i]: " << i << " " << dists[i] << std::endl;
	  std::cout << "x: ";
	  vector_out(std::cout,x,true);
	  std::cout << "data: ";
	  for(size_t jj=0;jj<this->n_params;jj++) {
	    std::cout << data_x(i,jj) << " ";
	  }
	  std::cout << std::endl;
	  std::cout << "scales: ";
	  vector_out(std::cout,scales,true);
	  O2SCL_ERR("Distance not finite.",o2scl::exc_einval);
	}
      }

      if (this->verbose>2) {
        std::cout << "interpm_idw::eval_err_index(): "
                  << "n_points,n_params,n_outputs: "
                  << this->n_points << " " << this->n_params << " "
                  << this->n_outputs
                  << std::endl;
      }

      // Find closest points, note that index is automatically resized
      // by the vector_smallest_index function
      o2scl::vector_smallest_index<std::vector<double>,double,
                                   std::vector<size_t> >
        (dists,points+1+n_extra,index);

      if (this->verbose>2) {
        std::cout << "  Indexes of closest points: ";
        o2scl::vector_out(std::cout,points+1+n_extra,index,true);
	std::cout << "  Distances:" << std::endl;
        for(size_t kk=0;kk<points+1+n_extra;kk++) {
          std::cout << "  " << kk << " " << dists[index[kk]] << std::endl;
        }
      }
    
      if (n_extra>0) {
        // Try to remove degenerate points to ensure accurate
        // interpolation
        bool found=true;
        while (found==true) {
          found=false;
          // Find degenerate points and remove them. Note that this
          // algorithm stops removing degenerate points if index.size()
          // has a size equal to points+1. This means that some
          // degeneracies may still remain if n_extra is not
          // sufficiently large.
          for(size_t j=0;j<points+1+n_extra;j++) {
            for(size_t k=j;k<points+1+n_extra;k++) {
              double dist_jk=dist(j,k);
              if (index.size()>points+1 && dist_jk<min_dist) {
                found=true;
                if (this->verbose>2) {
                  std::cout << "Erasing: " << j << std::endl;
                }
                index.erase(index.begin()+j);
              }
            }
          }
        }
      }

      if (this->verbose>2) {
	std::cout << "Closest: ";
	for(size_t j=0;j<this->n_params;j++) {
	  std::cout << data_x(index[0],j) << " ";
	}
	std::cout << std::endl;
	for(size_t j=0;j<this->n_outputs;j++) {
	  std::cout << data_y(index[0],j) << " ";
	}
	std::cout << std::endl;
      }
      
      if (dists[index[0]]<=0.0) {

        // If the closest distance is zero, just set the values and
        // errors
        for(size_t k=0;k<this->n_outputs;k++) {
          val[k]=data_y(index[0],k);
	  if (!std::isfinite(val[k])) {
	    std::cout << "n_extra,min_dist,k,n_outputs,val[k]: "
		      << n_extra << " " << min_dist << " " << k << " "
		      << this->n_outputs << " " << val[k] << std::endl;
	    O2SCL_ERR("Infinite value in interpm_idw() 1.",
		      o2scl::exc_efailed);
	  }
          err[k]=0.0;
          extrap[k]=0.0;
        }
        return 0;

      } else {
      
        for(size_t k=0;k<this->n_outputs;k++) {

          if (this->verbose>2) {
            std::cout << "  Output quantity " << k << " of " << this->n_outputs
                      << std::endl;
          }
        
          std::vector<double> vals(points+1);

          /// Distances between the selected points and the first point
          std::vector<double> dists_bw;
          /// Distances from the selected point to the user point
          std::vector<double> dists_to;
        
          // We construct points+1 estimates of the result and take the
          // average. Use the standard deviation for the uncertainty.
          for(size_t j=0;j<points+1;j++) {

            if (this->verbose>2) {
              std::cout << "  Point " << j << " of " << points+1
                        << std::endl;
            }
	    
            // Compute normalization
            double norm=0.0;
            for(size_t i=0;i<points+1;i++) {
              if (i!=j) norm+=1.0/dists[index[i]];
              dists_to.push_back(dists[index[i]]);
              if (i!=0) {
                double d=dist(index[0],index[i]);
                if (d!=0.0) {
                  dists_bw.push_back(d);
                }
              }
            }

            if (this->verbose>2) {
              std::cout << "  Norm: " << norm << std::endl;
            }
	    
            // Compute the inverse-distance weighted average
            vals[j]=0.0;
            for(size_t i=0;i<points+1;i++) {
              if (i!=j) {
                vals[j]+=data_y(index[i],k)/dists[index[i]];
                if (this->verbose>2) {
                  std::cout << "value, 1.0/dist: "
                            << data_y(index[i],k) << " "
                            << 1.0/dists[index[i]]
                            << std::endl;
                }
              }
            }

            vals[j]/=norm;
	    if (!std::isfinite(vals[j])) {
	      std::cout << "j,n_extra,min_dist,n_outputs,norm,vals[j]: "
			<< j << " " << n_extra << " " << min_dist
			<< " " << this->n_outputs << " " << norm << " "
			<< vals[j] << std::endl;
	      o2scl::vector_out(std::cout,dists,true);
	      O2SCL_ERR("Infinite value in interpm_idw 2.",
			o2scl::exc_efailed);
	    }

          }

          // Instead of using the average, we report the value as the
          // last element in the array, which is the interpolated
          // value from the closest points
          val[k]=vals[points];
	  
          err[k]=o2scl::vector_stddev(vals);

	  if (this->verbose>2) {
	    std::cout << "  Value, error: " << val[k] << " "
		      << err[k] << std::endl;
	  }
          
          if (dists_bw.size()>0) {
            double max_bw=vector_max_value<std::vector<double>,
                                           double>(dists_bw);
            double min_to=vector_min_value<std::vector<double>,
                                           double>(dists_to);
            extrap[k]=min_to/max_bw;
          } else {
            extrap[k]=0.0;
          }

          if (this->verbose>2) {
            std::cout << "Final value, err: " << val[k] << " "
                      << err[k] << std::endl;
          }

          // Proceed to the next dimension (end of loop over k)
        }

      }

      return 0;
    }
    
    /** \brief Perform the interpolation over all the functions
        with uncertainties
    */
    template<class vec2_t, class vec3_t, class vec4_t>
    int eval_unc_tl(const vec2_t &x, vec3_t &val, vec4_t &err) {
      std::vector<size_t> index;
      return eval_unc_tl_index(x,val,err,index);
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

        \verbatim embed:rst
        .. todo:: 

           - Use the mechanism provided by <tt>n_extra</tt> above
             to remove degenerate points. 

           - Future: This function requires an extra copy from
             "ders" to "ders2" which could be removed.

        \endverbatim
    */
    template<class vec3_t>
    void derivs_err(size_t func_index, size_t point_index, 
                    vec3_t &derivs, vec3_t &errs) const {
      
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_idw::derivs_err().",
                  exc_einval);
      }
      
      // Set x equal to the specified point
      ubvector x(this->n_params);
      for(size_t i=0;i<this->n_params;i++) {
        x[i]=data_x(point_index,i);
      }
      // Set f equal to the value of the function at the specified point
      double f=data_y(point_index,func_index);

      // The linear solver
      o2scl_linalg::linear_solver_HH<> lshh;
    
      // Compute distances
      std::vector<double> dists(this->n_points);
      for(size_t i=0;i<this->n_points;i++) {
        dists[i]=dist(i,x);
      }
  
      // Find closest (but not identical) points

      std::vector<size_t> index;
      size_t max_smallest=(this->n_params+2)*2;
      if (max_smallest>this->n_points) max_smallest=this->n_points;
      if (max_smallest<this->n_params+1) {
        O2SCL_ERR("Couldn't find enough nearby points.",
                  o2scl::exc_einval);
      }

      if (this->verbose>0) {
        std::cout << "max_smallest: " << max_smallest << std::endl;
      }
      
      o2scl::vector_smallest_index<std::vector<double>,double,
                                   std::vector<size_t> >
        (dists,max_smallest,index);

      if (this->verbose>0) {
        for(size_t i=0;i<index.size();i++) {
          std::cout << "index[" << i << "] = " << index[i] << " "
                    << dists[index[i]] << std::endl;
        }
      }
      
      std::vector<size_t> index2;
      for(size_t i=0;i<max_smallest;i++) {
        if (dists[index[i]]>0.0) {
          index2.push_back(index[i]);
          if (index2.size()==this->n_params+1) i=max_smallest;
        }
      }
      if (index2.size()<this->n_params+1) {
        O2SCL_ERR("Couldn't find enough nearby points (2).",
                  o2scl::exc_einval);
      }

      if (this->verbose>0) {
        for(size_t i=0;i<index2.size();i++) {
          std::cout << "index2[" << i << "] = " << index2[i] << " "
                    << dists[index2[i]] << std::endl;
        }
      }
      
      // Unit vector storage
      std::vector<ubvector> units(this->n_params+1);
      // Difference vector norms
      std::vector<double> diff_norms(this->n_params+1);
      // Storage for the derivative estimates
      std::vector<ubvector> ders(this->n_params+1);
      // Matrix of unit vectors
      ubmatrix m(this->n_params,this->n_params);
      // Vector of function value differences
      ubvector v(this->n_params);
      // Rearranged derivative object
      std::vector<ubvector> ders2(this->n_params);
      
      for(size_t i=0;i<this->n_params+1;i++) {

        // Assign unit vector elements
        units[i].resize(this->n_params);
        for(size_t j=0;j<this->n_params;j++) {
          units[i][j]=data_x(index2[i],j)-x[j];
        }

        // Normalize the unit vectors
        diff_norms[i]=o2scl::vector_norm<ubvector,double>(units[i]);
        for(size_t j=0;j<this->n_params;j++) {
          units[i][j]/=diff_norms[i];
        }

      }

      // Verbose output of the closest points and their norms
      if (this->verbose>0) {
        std::cout << "Point:     ";
        for(size_t i=0;i<this->n_params;i++) {
          std::cout << x[i] << " ";
        }
        std::cout << f << std::endl;
        for(size_t j=0;j<this->n_params+1;j++) {
          std::cout << "Closest: " << j << " " << index2[j] << " ";
          for(size_t i=0;i<this->n_params;i++) {
            std::cout << data_x(index2[j],i) << " ";
          }
          std::cout << data_y(index2[j],func_index) << " "
                    << diff_norms[j] << std::endl;
        }
        for(size_t j=0;j<this->n_params+1;j++) {
          std::cout << "diff_norm: " << j << " " << diff_norms[j]
                    << std::endl;
        }
        // End of verbose output
      }
    
      // Go through each set of points
      for(size_t i=0;i<this->n_params+1;i++) {

        ders[i].resize(this->n_params);

        // Construct the matrix and vector for the solver
        size_t jj=0;
        for(size_t j=0;j<this->n_params+1;j++) {
          if (j!=i) {
            for(size_t k=0;k<this->n_params;k++) {
              m(jj,k)=units[j][k];
            }
            v[jj]=(data_y(index2[j],func_index)-f)/diff_norms[j];
            jj++;
          }
        }

        // Solve to compute the derivatives
        if (this->verbose>0) {
          std::cout << "m:" << std::endl;
          o2scl::matrix_out(std::cout,this->n_params,this->n_params,m);
          std::cout << "v:" << std::endl;
          o2scl::vector_out(std::cout,this->n_params,v,true);
        }
        lshh.solve(this->n_params,m,v,ders[i]);
        if (this->verbose>0) {
          std::cout << "Derivs:  " << i << " ";
          std::cout.setf(std::ios::showpos);
          for(size_t j=0;j<this->n_params;j++) {
            std::cout << ders[i][j] << " ";
          }
          std::cout.unsetf(std::ios::showpos);
          std::cout << std::endl;
        }

        // Go to next derivative estimate
      }
      
      for(size_t i=0;i<this->n_params;i++) {

        // Rearrange derivatives
        ders2[i].resize(this->n_params+1);
        for(size_t j=0;j<this->n_params+1;j++) {
          ders2[i][j]=ders[j][i];
        }

        // Compute mean and standard deviation
        derivs[i]=o2scl::vector_mean(ders2[i]);
        errs[i]=o2scl::vector_stddev(ders2[i]);
      }
      
      return;
    }
    //@}
    
  protected:
    
    /// The copy of the x data 
    mat_x_t data_x;
    /// The copy of the y data
    mat_y_t data_y;
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
      for(size_t i=0;i<this->n_params;i++) {
        ret+=pow((x[i]-data_x(index,i))/scales[i%nscales],dist_expo);
      }
      return sqrt(ret);
    }

    /** \brief Compute the distance between two points in the
        data set
    */
    double dist(size_t j, size_t k) const {
      double ret=0.0;
      size_t nscales=scales.size();
      for(size_t i=0;i<this->n_params;i++) {
        ret+=pow((data_x(j,i)-data_x(k,i))/scales[i%nscales],dist_expo);
      }
      return sqrt(ret);
    }
    //@}
    
  };
  
}
    
#endif

