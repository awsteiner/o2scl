/*
  -------------------------------------------------------------------
  
  Copyright (C) 2009-2017, Marco Cammarata and Andrew W. Steiner
  
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
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>

#include <o2scl/smooth_gsl.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  size_t N=100;
  gsl_vector *x, *y, *e, *ys;
  x=gsl_vector_alloc(N);
  y=gsl_vector_alloc(N);
  e=gsl_vector_alloc(N);
  ys=gsl_vector_alloc(N);

  for(size_t i=0;i<N;i++) {
    gsl_vector_set(x,i,((double)i)*0.05);
    gsl_vector_set(y,i,sin(gsl_vector_get(x,i))+
		   0.05*sin(1.0e8*gsl_vector_get(x,i)));
    gsl_vector_set(e,i,0.05);
  }
  
  smooth_gsl gs(x);
  
  gs.smooth_data(y,e,ys);
  
  for(size_t i=0;i<N;i++) {
    if (i>3 && i!=63) {
      t.test_rel(gsl_vector_get(ys,i),
		 sin(gsl_vector_get(x,i)),0.02,"smooth_gsl");
    }
  }

  t.report();
  return 0;
}
