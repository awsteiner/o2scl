/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <iostream>
#include <hal/table.h>
#include <hal/minimize.h>
#include <hal/base_ioc.h>

using namespace std;
using namespace hal;

int main(void) {
  table t;
  t.line_of_names("x lb clb co cco");
  cout.setf(ios::scientific);
  for(double x=3.0;x<=7.0;x+=0.02) {
    double line[5]={x,lower_bound(x,5.0,1.0,10.0),
		    cont_lower_bound(x,5.0,1.0,10.0),
		    constraint(x,5.0,1.0,10.0),
		    cont_constraint(x,5.0,1.0,10.0)};
    t.line_of_data(5,line);
  }
  
  base_ioc bio;
  collection co;
  text_out_file *tof=new text_out_file("ex_cnstr.out");
  co.out_one(tof,"table","t",&t);
  delete tof;
  return 0;
}
  
