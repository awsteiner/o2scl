/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_PINSIDE_H
#define O2SCL_PINSIDE_H

/** \file pinside.h
    \brief File defining \ref o2scl::pinside 
*/

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/vector.h>

namespace o2scl {

  /** \brief Test line intersection and point inside polygon

      This is a fast and dirty implementation of the point inside
      polygon test from Jerome L. Lewis, SIGSCE Bulletin, 34 (2002)
      81.

      Note that an error in that article ("count-" should have been
      "count--") has been corrected here. 

      \future The inside() functions actually copy the points twice.
      This can be made more efficient.

      \comment
      See also
      http://www.ecse.rpi.edu/Homepages/wrf/Research/
      Short_Notes/pnpoly.html#Polyhedron
      which suggests the following code
      int pnpoly(int nvert, float *vertx, float *verty, float testx, 
      float testy) {
      int i, j, c = 0;
      for (i = 0, j = nvert-1; i < nvert; j = i++) {
      if ( ((verty[i]>testy) != (verty[j]>testy)) &&
      (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / 
      (verty[j]-verty[i]) + vertx[i]) )
      c = !c;
      }
      return c;
      }
      \endcomment

  */
  class pinside {

  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    
  protected:
  
    /// Internal point definition for \ref pinside
    struct point {
      double x;
      double y;
    };
  
    /// Internal line definition for \ref pinside
    struct line {
      point p1;
      point p2;
    };
  
    /// Test if line segments \c P and \c Q intersect
    int intersect(line P, line Q);
  
    /// Test if point \c t is inside polygon \c p of size \c N
    int inside(point t, point p[], int N);
  
  public:

    /** \brief Determine if two line segments intersect

	The function returns 1 if the line segment
	determined by the endpoints \f$ (x_1,y_1) \f$ 
	and \f$ (x_2,y_2) \f$ and the line segment 
	determined by the endpoints \f$ (x_3,y_3) \f$ 
	and \f$ (x_4,y_4) \f$ intersect, and 0 otherwise.
    */
    int intersect(double x1, double y1, double x2, double y2, 
		  double x3, double y3, double x4, double y4) {
      line P={{x1,y1},{x2,y2}};
      line Q={{x3,y3},{x4,y4}};
      return intersect(P,Q);
    }

    /** \brief Determine if point (x,y) is inside a polygon

	This returns 1 if the point (x,y) is inside the polygon
	defined by \c xa and \c ya, and 0 otherwise. 

	Note that if the point (x,y) is exactly on the polygon,
	then the result of this function is not well-defined
	and it will return either 0 or 1. 
    */
    int inside(double x, double y, const ubvector &xa,
	       const ubvector &ya);
  
    /** \brief Determine if point (x,y) is inside a polygon

	This returns 1 if the point (x,y) is inside the polygon
	defined by \c xa and \c ya, and 0 otherwise. 

	The parameter \c n should be the number of polygon points
	specified in vectors \c xa and \c ya.

	Note that if the point (x,y) is exactly on the polygon,
	then the result of this function is not well-defined
	and it will return either 0 or 1.
    */
    template<class vec_t>
      int inside(double x, double y, size_t n, const vec_t &xa, 
		 const vec_t &ya) {
      
      size_t ix;
      point t, *p=new point[n+1];
      t.x=x;
      t.y=y;

      // We have to copy the vectors so we can rearrange them because
      // they are const
      ubvector xb(n), yb(n);
      vector_copy(n,xa,xb);
      vector_copy(n,ya,yb);

      // Ensure that (yb[0],ya[0]) is the point with the smallest x
      // coordinate among all the points with the smallest y coordinate
      double xmin=xb[0];
      double ymin=yb[0];
      ix=0;
      for(size_t i=0;i<n;i++) {
	if (yb[i]<ymin) {
	  ymin=yb[i];
	  ix=i;
	}
      }
      for(size_t i=0;i<n;i++) {
	if (yb[i]==ymin && xb[i]<xmin) {
	  xmin=xb[i];
	  ix=i;
	}
      }
      vector_rotate<vec_t,double>(n,xb,ix);
      vector_rotate<vec_t,double>(n,yb,ix);

      // Copy to p[]
      for(size_t i=0;i<n;i++) {
	p[i+1].x=xb[i];
	p[i+1].y=yb[i];
      }
    
      int ret=inside(t,p,n);
      delete[] p;

      return ret;
    }

    /// Perform some simple testing
    int test(test_mgr &t);

  };

}

#endif
