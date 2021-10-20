/*
  -------------------------------------------------------------------
  
  Copyright (C) 2010-2021, Andrew W. Steiner
  
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

#include <o2scl/fract.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int o2scl::nrf_z4m1(size_t nv, const ubvector &x, ubvector &f,
                    ubmatrix &J) {
  double x2=x[0]*x[0];
  double y2=x[1]*x[1];
  f[0]=-1+x2*x2-6*x2*y2+y2*y2;
  f[1]=4*x[0]*x[1]*(x2-y2);
  J(0,0)=4*x2*x[0]-12*x[0]*y2;
  J(0,1)=-12*x2*x[1]+4*y2*x[1];
  J(1,0)=8*x2*x[1]+4*x[1]*(x2-y2);
  J(1,1)=-8*x[0]*y2+4*x[0]*(x2-y2);
  return 0;
}
  
int o2scl::itf_mandel(ubvector &z, ubvector c) {
  double zr=c[0]+z[0]*z[0]-z[1]*z[1];
  double zi=c[1]+2.0*z[0]*z[1];
  z[0]=zr;
  z[1]=zi;
  return 0;
}
  
