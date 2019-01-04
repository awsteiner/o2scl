/*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2019, Andrew W. Steiner
  
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
#include <o2scl/reaction_lib.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  
  // Create a temporary file with data
  ofstream fout("reaction_temp.dat");
  fout << "1" << endl;
  fout << "         n    p                            wc07w     0.00000e+00" 
       << endl;
  fout << "-6.786180e+00 0.000000e+00 0.000000e+00 0.000000e+00" << endl;
  fout << " 0.000000e+00 0.000000e+00 0.000000e+00" << endl;
  fout << "1" << endl;
  fout << "         t  he3                            wc07w     0.00000e+00" 
       << endl;
  fout << "-2.014510e+01 0.000000e+00 0.000000e+00 0.000000e+00" << endl;
  fout << " 0.000000e+00 0.000000e+00 0.000000e+00" << endl;
  fout.close();

  reaction_lib r;
  int ret=r.read_file_reaclib2("reaction_temp.dat");
  cout << ret << " " << r.lib.size() << endl;
  cout << endl;
  
  t.test_gen(r.lib.size()==2,"size");

  t.report();
  return 0;
}
