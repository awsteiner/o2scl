/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/exception.h>

using namespace std;
using namespace o2scl;

int main(void) {
  cout.setf(ios::scientific);
  test_mgr t;
  t.set_output_level(1);
  err_hnd_cpp ee;
  
  err_hnd=&ee;

  cout << err_hnd->get_str() << endl;
  try {
    O2SCL_ERR("Exception test",1);
  } catch (std::exception) {
    cout << err_hnd->get_str() << endl;
    err_hnd->reset();
  }
  cout << err_hnd->get_str() << endl;

  t.report();
  return 0;
}
