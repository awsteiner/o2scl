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

// sphinx-example-start
/* Example: ex_string.cpp
   -------------------------------------------------------------------
   Demonstrate the string conversion classes and functions
*/

#include <o2scl/test_mgr.h>
#include <o2scl/format_float.h>
#include <o2scl/misc.h>
#include <o2scl/funct_to_fp.h>

using namespace std;
using namespace o2scl;

int main(void) {
  
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);
  
  format_float ff;
  ff.latex_mode();

  cout << "Output a LaTeX table of numbers: " << endl;
  cout << endl;
  cout << "\\begin{table}" << endl;
  cout << "\\begin{tabular}{cc}" << endl;
  cout << "\\hline" << endl;
  cout << "Column 1 & Column 2 \\\\" << endl;
  cout << "\\hline" << endl;
  cout << ff.convert(1.1e-6) << " & " 
       << ff.convert(-4.0) << " \\\\" << endl;
  cout << ff.convert(134.567) << " & " 
       << ff.convert(1568234034) << " \\\\" << endl;
  cout << "\\hline" << endl;
  cout << "\\end{tabular}" << endl;
  cout << "\\end{table}" << endl;
  cout << endl;

  cout << "Function function_to_double():" << endl;
  double x=function_to_double("cyl_bessel_i(2,pi)");
  cout << x << endl;
  t.test_rel(x,2.618495,1.0e-6,"BesselI(2,x)");

  cout << "Function string_to_uint_list():\n  " << endl;
  vector<size_t> ulist;
  string_to_uint_list("1,2-4,6,10-11",ulist);
  vector_out(cout,ulist,true);
  t.test_gen(ulist==vector<size_t>({1,2,3,4,6,10,11}),"uint list");

  t.report();
  return 0;
}

