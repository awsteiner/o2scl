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
#include <o2scl/err_hnd.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int function1() {
  O2SCL_ERR("This is the third error.",2);
  return 0;
}

int function2() {
  O2SCL_ERR("This is the fourth error.",3);
  return 0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  if (false) {

    O2SCL_ERR("This is the first error.",1);
    cout << err_hnd->get_str() << endl;
    
    O2SCL_ERR("Adding a second error.",2);
    cout << err_hnd->get_str() << endl;
    
    cout << function1() << endl;
    cout << err_hnd->get_str() << endl;
    
    cout << function2() << endl;
    cout << err_hnd->get_str() << endl;
    cout << endl;

  }

  t.report();
  return 0;
}
