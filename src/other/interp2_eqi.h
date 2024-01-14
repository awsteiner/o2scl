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
#ifndef O2SCL_INTERP2_EQI_H
#define O2SCL_INTERP2_EQI_H

/** \file interp2_eqi.h
    \brief File defining \ref o2scl::interp2_eqi
*/

#include <iostream>
#include <string>

namespace o2scl {

  /** \brief Two-dimensional interpolation for equally-spaced intervals

      \note This class is unfinished.

      This implements the relations from Abramowitz and Stegun:
      \f[
      f(x_0+p h,y_0+q k)=
      \f]
      3-point
      \f[
      (1-p-q) f_{0,0}+p f_{1,0}+q f_{0,1}
      \f]
      4-point
      \f[
      (1-p)(1-q) f_{0,0}+p(1-q)f_{1,0}+q(1-p)f_{0,1}+pqf_{1,1}
      \f]
      6-point
      \f[
      \frac{q(q-1)}{2}f_{0,-1}+\frac{p(p-1)}{2}f_{-1,0}+
      (1+pq-p^2-q^2)f_{0,0}+\frac{p(p-2q+1)}{2}f_{1,0}+
      \frac{q(q-2p+1)}{2}f_{0,1}+pqf_{1,1}
      \f]
  */
  class interp2_eqi {

  public:
    
    interp2_eqi();

    /** \brief Perform the 2-d interpolation 
     */
    double interp(double x, double y);

    /** \brief Offset in x-direction 
     */
    double xoff;
    
    /** \brief Offset in y-direction 
     */
    double yoff;

    /** \brief Set the interpolation type

        - 3: 3-point
        - 4: 4-point
        - 6: 6-point (default)
     */
    int set_type(int type) {
      if (type<=0) itype=6;
      else if (type<=3) itype=3;
      else if (type>=5) itype=6;
      else type=4;
      return 0;
    }

  protected:

    int itype;

  };
  
}

#endif



