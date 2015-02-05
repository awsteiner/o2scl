/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#ifndef O2SCL_INTERP2_NEIGH_H
#define O2SCL_INTERP2_NEIGH_H

/** \file interp2_neigh.h
    \brief File defining \ref o2scl::interp2_neigh
*/

#include <iostream>
#include <string>
#include <cmath>

#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Nearest-neighbor interpolation in two dimensions

      This class performs nearest-neighbor interpolation when the data
      points are not arranged in a specified order (i.e. not on a
      grid). For a set of data \f$ {x_i,y_i,f_i} \f$, the value of \f$
      f \f$ is predicted given a new value of x and y. Distances are
      determined with
      \f[
      d_{ij} = \sqrt{\left(\frac{x_i-x_j}{\Delta x}\right)^2 +
      \left(\frac{y_i-y_j}{\Delta y}\right)^2}
      \f]
      The values \f$ \Delta_x \f$ and \f$ \Delta_y \f$ are specified
      in \ref x_scale and \ref y_scale, respectively. If these values
      are negative (the default) then they are computed with \f$
      \Delta x = x_{\mathrm{max}}-x_{\mathrm{min}} \f$ and \f$ \Delta
      y = y_{\mathrm{max}}-y_{\mathrm{min}} \f$ .

      This class stores pointers to the data, not a copy. The data can
      be changed between interpolations without an additional call to
      \ref set_data(), but the scales may need to be recomputed with
      \ref compute_scale().

      The vector type can be any type with a suitably defined \c
      operator[].
      
      \note This class operates by performing a \f$ {\cal O}(N) \f$
      brute-force search to find the closest points. 

      \future Make a parent class for this and \ref o2scl::interp2_planar.
  */
  template<class vec_t> class interp2_neigh {

  protected:

    /// The scale in the x direction
    double dx;

    /// The scale in the y direction
    double dy;

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    
    interp2_neigh() {
      data_set=false;
      x_scale=-1.0;
      y_scale=-1.0;
      dx=0.0;
      dy=0.0;
    }

    /// The user-specified x scale (default -1)
    double x_scale;

    /// The user-specified y scale (default -1)
    double y_scale;

    /// Find scaling
    void compute_scale() {
      if (x_scale<0.0) {
	double minx=(*ux)[0], maxx=(*ux)[0];
	for(size_t i=1;i<np;i++) {
	  if ((*ux)[i]<minx) minx=(*ux)[i];
	  if ((*ux)[i]>maxx) maxx=(*ux)[i];
	}
	dx=maxx-minx;
      } else {
	dx=x_scale;
      }
      if (y_scale<0.0) {
	double miny=(*uy)[0], maxy=(*uy)[0];
	for(size_t i=1;i<np;i++) {
	  if ((*uy)[i]<miny) miny=(*uy)[i];
	  if ((*uy)[i]>maxy) maxy=(*uy)[i];
	}
	dy=maxy-miny;
      } else {
	dy=y_scale;
      }

      if (dx<=0.0 || dy<=0.0) {
	O2SCL_ERR("No scale in interp2_planar::set_data().",exc_einval);
      }

      return;
    }

    /** \brief Initialize the data for the neigh interpolation

	This function will call the error handler if \c n_points
	is zero.
     */
    void set_data(size_t n_points, vec_t &x, vec_t &y, vec_t &f) {
      if (n_points<1) {
	O2SCL_ERR2("Must provide at least one point in ",
		   "interp2_neigh::set_data()",exc_efailed);
      }
      np=n_points;
      ux=&x;
      uy=&y;
      uf=&f;
      data_set=true;

      compute_scale();

      return;
    }
    
    /** \brief Perform the interpolation 
     */
    double eval(double x, double y) const {
      double x1, y1, f;
      size_t i1;
      eval_point(x,y,f,i1,x1,y1);
      return f;
    }

    /** \brief Perform the interpolation 
    */
    double operator()(double x, double y) const {
      return eval(x,y);
    }

    /** \brief Perform the planar interpolation using the first two
	elements of \c v as input
    */
    template<class vec2_t> double operator()(vec2_t &v) const {
      return eval(v[0],v[1]);
    }

    /** \brief Interpolation returning the closest point 

	This function interpolates \c x and \c y into the data
	returning \c f. It also returns the closest x- and y-values
	found.
    */
    void eval_point(double x, double y, double &f,
		    size_t &i1, double &x1, double &y1) const {
      
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp_planar::eval_points().",
		  exc_einval);
      }

      // Exhaustively search the data
      i1=0;
      double dist_min=pow((x-(*ux)[i1])/dx,2.0)+pow((y-(*uy)[i1])/dy,2.0);
      for(size_t index=1;index<np;index++) {
	double dist=pow((x-(*ux)[index])/dx,2.0)+pow((y-(*uy)[index])/dy,2.0);
	if (dist<dist_min) {
	  i1=index;
	  dist_min=dist;
	}
      }

      // Return the function value

      f=(*uf)[i1];

      return;
    }
    
#ifndef DOXYGEN_INTERNAL

  protected:
    
    /// The number of points
    size_t np;
    /// The x-values
    vec_t *ux;
    /// The y-values
    vec_t *uy;
    /// The f-values
    vec_t *uf;
    /// True if the data has been specified
    bool data_set;
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



