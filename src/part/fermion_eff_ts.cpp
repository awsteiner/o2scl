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
#include <o2scl/fermion_eff.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/classical.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  fermion_eff fe;
  fermion f;
  part_calibrate_class pcc;

  /*
  double q=pcc.part_calibrate<fermion,fermion_eff>
    (f,fe,true,"../../data/o2scl/fermion_deriv_cal.o2",false,true,1,true);
  
  t.test_rel(q,0.0,0.5,"calibrate");
  */
  
  t.report();
  return 0;
}
  

