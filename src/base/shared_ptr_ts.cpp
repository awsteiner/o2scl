/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2014, Andrew W. Steiner

  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------
*/
#include <iostream>

#include <o2scl/test_mgr.h>
#include <o2scl/shared_ptr.h>

using namespace std;
using namespace o2scl;

int x;

class obj {
public:
  obj() {
    cout << "Constructor." << endl;
    x*=2;
  }
  ~obj() {
    cout << "Destructor." << endl;
    x*=3;
  }
};

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  x=1;

  // Test the aliased version
  {
    cout << "Point 1" << endl;
    o2_shared_ptr<obj>::type ptr(new obj);
    t.test_gen(x==2,"Constructor");
    {
      cout << "Point 2" << endl;
      o2_shared_ptr<obj>::type ptr2=ptr;
      cout << "Point 3" << endl;
    }
    cout << "Point 4" << endl;
  }
  t.test_gen(x==6,"Destructor");
  cout << "Point 5" << endl;
  
  t.report();
  return 0;
}

