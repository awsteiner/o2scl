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
#ifndef O2SCL_INTERP2_PLANAR_H
#define O2SCL_INTERP2_PLANAR_H

#include <iostream>
#include <string>
#include <cmath>

#include <gsl/gsl_combination.h>

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Interpolate among two independent variables with planes

      This class is experimental.

      This is an analog of 1-d linear interpolation for two
      dimensions, which is particularly useful with the data points
      are not arranged in a specified order (i.e. on a grid). For a
      set of data \f$ {x_i,y_i,f_i} \f$, the value of \f$ f \f$ is
      predicted given a new value of x and y. In contrast to \ref
      o2scl::interp2_seq, the data need not be presented in a grid.
      This interpolation is performed by finding the plane that goes
      through three closest points in the data set. Distances are
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

      If the x- and y-values of the entire data set lie on a line,
      then the interpolation will fail and the error handler will be
      called. Colinearity is defined by a threshold, \ref thresh
      which defaults to \f$ 10^{-12} \f$.

      \comment
      (The following isn't true any more because of the scale
      computation)

      There is no caching so the numeric values of the data may be
      freely changed between calls to interp()
      \endcomment

      The vector type can be any type with a
      suitably defined function \c operator[]. 

      \note This class operates by performing a \f$ {\cal O}(N) \f$
      brute-force search to find the three closest points to the
      user-specified location which are not colinear.

      \future The generalization to more dimensions might be 
      straightforward.
  */
  template<class vec_t> class interp2_planar {

  protected:

    /// The scale in the x direction
    double dx;

    /// The scale in the y direction
    double dy;

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    
    interp2_planar() {
      data_set=false;
      thresh=1.0e-12;
      x_scale=-1.0;
      y_scale=-1.0;
      dx=0.0;
      dy=0.0;
    }

    /// Threshold for colinearity (default \f$ 10^{-12} \f$)
    double thresh;

    /// The user-specified x scale (default -1)
    double x_scale;

    /// The user-specified y scale (default -1)
    double y_scale;

    /** \brief Initialize the data for the planar interpolation
     */
    void set_data(size_t n_points, vec_t &x, vec_t &y, vec_t &f) {
      if (n_points<3) {
	O2SCL_ERR2("Must provide at least three points in ",
		       "interp2_planar::set_data()",exc_efailed);
      }
      np=n_points;
      ux=&x;
      uy=&y;
      uf=&f;
      data_set=true;

      // Find scaling
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
    
    /** \brief Perform the planar interpolation 
    */
    double eval(double x, double y) const {
      double x1, x2, x3, y1, y2, y3, f;
      size_t i1, i2, i3;
      eval_points(x,y,f,i1,x1,y1,i2,x2,y2,i3,x3,y3);
      return f;
    }

    /** \brief Perform the planar interpolation 
    */
    double operator()(double x, double y) const {
      return eval(x,y);
    }

    /** \brief Planar interpolation returning the closest points 

	This function interpolates \c x and \c y into the data
	returning \c f. It also returns the three closest x- and
	y-values used for computing the plane.
    */
    void eval_points(double x, double y, double &f,
		     size_t &i1, double &x1, double &y1, 
		     size_t &i2, double &x2, double &y2, 
		     size_t &i3, double &x3, double &y3) const {
      
      if (data_set==false) {
	O2SCL_ERR("Data not set in interp2_planar::eval_points().",
		  exc_einval);
      }

      // First, we just find the three closest points by
      // exhaustively searching the data
      
      // Put in initial points
      i1=0; i2=1; i3=2;
      double c1=sqrt(pow((x-(*ux)[0])/dx,2.0)+pow((y-(*uy)[0])/dy,2.0));
      double c2=sqrt(pow((x-(*ux)[1])/dx,2.0)+pow((y-(*uy)[1])/dy,2.0));
      double c3=sqrt(pow((x-(*ux)[2])/dx,2.0)+pow((y-(*uy)[2])/dy,2.0));

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
	double c4=sqrt(pow((x-(*ux)[i4])/dx,2.0)+pow((y-(*uy)[i4])/dy,2.0));
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

      // Solve for denominator:
      double a, b, c, den, z1, z2, z3;
      x1=(*ux)[i1];
      x2=(*ux)[i2];
      x3=(*ux)[i3];
      y1=(*uy)[i1];
      y2=(*uy)[i2];
      y3=(*uy)[i3];
      den=(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);

      // If the points are colinear (if den is too small), then
      // we do a complete sorting of the distance between the 
      // user-specified point and the data, and then go through
      // all combinations to find the closest triplet that is
      // not colinear.

      if (fabs(den)<thresh) {
	
	ubvector dist(np);
	ubvector_size_t order(np);
	for(size_t i=0;i<np;i++) {
	  dist[i]=sqrt(pow((x-(*ux)[i])/dx,2.0)+pow((y-(*uy)[i])/dy,2.0));
	}
	o2scl::vector_sort_index(np,dist,order);

	{
	  gsl_combination *cc=gsl_combination_alloc(np,3);

	  gsl_combination_init_first(cc);

	  bool done=false;
	  while (done==false) {
	    int status=gsl_combination_next(cc);

	    if (status!=gsl_failure) {
	      i1=order[gsl_combination_get(cc,0)];
	      i2=order[gsl_combination_get(cc,1)];
	      i3=order[gsl_combination_get(cc,2)];
	      x1=(*ux)[i1];
	      x2=(*ux)[i2];
	      x3=(*ux)[i3];
	      y1=(*uy)[i1];
	      y2=(*uy)[i2];
	      y3=(*uy)[i3];
	      den=(-x2*y1+x3*y1+x1*y2-x3*y2-x1*y3+x2*y3);

	      if (fabs(den)>thresh) {
		done=true;
	      }

	    } else {
	      /*
		If we've gone through the entire data set and the
		whole thing is colinear, then throw an exception.
	      */
	      O2SCL_ERR2("Interpolation failed in.",
			 "interp2_planar::eval_points().",exc_efailed);
	    }
	  }
	  
	  gsl_combination_free(cc);
	}

      }
      
      // Calculate the function value with the three closest
      // non-colinear points

      z1=(*uf)[i1];
      z2=(*uf)[i2];
      z3=(*uf)[i3];
      a=-(-y2*z1+y3*z1+y1*z2-y3*z2-y1*z3+y2*z3)/den;
      b=-(x2*z1-x3*z1-x1*z2+x3*z2+x1*z3-x2*z3)/den;
      c=-(x3*y2*z1-x2*y3*z1-x3*y1*z2+x1*y3*z2+x2*y1*z3-x1*y2*z3)/den;
      
      f=a*x+b*y+c;

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
    
    /// Swap points 1 and 2.
    int swap(size_t &i1, double &c1, size_t &i2, double &c2) const {
      int t;
      double tc;
      
      t=i1; tc=c1;
      i1=i2; c1=c2;
      i2=t; c2=tc;
      
      return 0;
    }
    
#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



