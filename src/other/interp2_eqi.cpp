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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/interp2_eqi.h>

using namespace std;
using namespace o2scl;

interp2_eqi::interp2_eqi() {
  itype=6;
}

double interp2_eqi::interp(double x, double y) {
  //  double p=xoff/h, q=yoff/k;
  /*
      \f[
      f(x_0+p h,y_0+q k)=
      \f]
      3-point
      \f[
      =(1-p-q) f_{0,0}+p f_{1,0}+q f_{0,1}
      \f]
      4-point
      \f[
      =(1-p)(1-q) f_{0,0)+p(1-q)f_{1,0}+q(1-p)f_{0,1}+pqf_{1,1}
      \f]
      6-point
      \f[
      =\frac{q(q-1)}{2}f_{0,-1}+\frac{p(p-1)}{2}f_{-1,0}+
      (1+pq-p^2-q^2)f_{0,0}+\frac{p(p-2q+1)}{2}f_{1,0}+
      \frac{q(q-2p+1)}{2}f_{0,1}+pqf_{1,1}
      \f]
  */

  //  return result;
  return 0.0;
}

