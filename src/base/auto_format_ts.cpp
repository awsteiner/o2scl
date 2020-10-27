/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020, Andrew W. Steiner
  
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
#include <o2scl/auto_format.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_auto_format;

int main(void) {

  test_mgr t;
  t.set_output_level(2);

  auto_format at;
  at.verbose=0;
  double d=3.0;

  //at << d << o2scl_auto_format::endl << d*d << o2scl_auto_format::endl;

  at << d << endo;
  
  at.start_table();
  at << "column_1" << "column_2" << endo;
  //at.debug_table();
  at << 2*d << 'x' << endo;
  //at.debug_table();
  at << 34 << 1.0e-200 << endo;
  //at.debug_table();
  at << 8*d << "blah" << 4*d << 5*d << endo;
  //at.debug_table();
  at.end_table();

  at << endo;
  
  at.start_table();
  at.table_lines=1;
  at << 4.0 << "x" << -1.0 << endo;
  //at.debug_table();
  at << 5.0 << "a b" << -2.0e-220 << "foo" << endo;
  //at.debug_table();
  at.end_table();
  
  at << 62 << endo;
  at << endo;
  at << 6.0 << endo;
  vector<double> vd={3,1,4,1};
  at << vd << endo;
  vector<int> vi={3,1,4,1};
  at << vi << endo;
  vector<std::string> vs={"this","is","a  test"};
  at << vs << endo;
  at.done();
  
  t.report();
  return 0;
}

