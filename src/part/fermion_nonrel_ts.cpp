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
#include <o2scl/test_mgr.h>

#include <o2scl/fermion_nonrel.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
  fermion_nonrel fnr;
  fermion f(1.0,2.0);

  part_calibrate_class pcc;

  // Third argument is false because we don't want to test pairs
  double v1=pcc.part_calibrate<fermion,fermion_nonrel>
    (f,fnr,false,true,"../../data/o2scl/fermion_nr_cal.o2",true,1,true);
  t.test_abs(v1,0.0,1.0e-12,"calibrate");

  t.set_output_level(2);
  t.report();

  return 0;
}
