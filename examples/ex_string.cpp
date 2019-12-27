/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

/* Example: ex_string_conv.cpp
   -------------------------------------------------------------------
   Demonstrate the string conversion classes and functions
*/

#include <o2scl/test_mgr.h>
#include <o2scl/format_float.h>
#include <o2scl/misc.h>

using namespace std;
using namespace o2scl;

int main(void) {
  
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  cout << "Convert several strings to numbers: " << endl;
  std::string d1="1.0e-6";
  std::string d2="1.0E-6";
  std::string d3="-0.001";
  std::string i1="1034";
  // We have to use the full o2scl::stod() specification to 
  // avoid confusion with std::stod()
  t.test_rel(o2scl::stod(d1),1.0e-6,1.0e-12,"stod 1");
  t.test_rel(o2scl::stod(d2),1.0e-6,1.0e-12,"stod 2");
  t.test_rel(o2scl::stod(d3),-1.0e-3,1.0e-12,"stod 3");
  t.test_gen(o2scl::stod(i1)==1034,"stod 4");
  cout << endl;

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

  t.report();
  return 0;
}
// End of example

