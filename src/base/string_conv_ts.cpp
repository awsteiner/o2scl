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
#include <o2scl/misc.h>
#include <o2scl/test_mgr.h>
#include <climits>

using namespace std;
using namespace o2scl;

int main(void) {
  double x=sqrt(2.0);
  int i=3;
  string s1="4", s2="1.7320508";
  test_mgr t;
  t.set_output_level(2);
  
  cout << "These are strings: ";
  cout << itos(i) << " " << dtos(x) << endl;
  cout << "These are numbers: ";
  cout << o2scl::stoi(s1) << " " << o2scl::stod(s2) << endl;
  cout << endl;

  cout << "Testing incorrect inputs: " << endl;
  cout << "stoi(\"5.6\"): \t\t" << o2scl::stoi("5.6") << endl;
  cout << "stoi(\"987.6\"): \t\t" << o2scl::stoi("987.6") << endl;
  cout << "stoi(\"234.43e1zxepq\"): \t" << o2scl::stoi("234.43e1zxepq") << endl;
  cout << "stod(\"234.43e1zxepq\"): \t" << o2scl::stod("234.43e1zxepq") << endl;

  cout << "This demonstrates has_minus_sign(): " << endl;
  cout.setf(ios::scientific);
  double z1=0.0;
  double z2=-0.0;
  cout << z1 << " " << z2 << endl;
  cout << (z1>=0.0) << " " << (z2>=0.0) << endl;
  cout << has_minus_sign(&z1) << " " << has_minus_sign(&z2) << endl;
  
  cout << "Here." << endl;
  cout << stob(" true") << endl;
  cout << stob(" false") << endl;
  cout << stob("True") << endl;
  cout << stob("1") << endl;
  cout << stob(" 1") << endl;
  cout << stob(" 0") << endl;
  cout << stob("-1") << endl;

  cout << size_of_exponent(1.0e-111) << " " << 1.0e-111 << endl;
  cout << size_of_exponent(1.0e-11) << " " << 1.0e-11 << endl;
  cout << size_of_exponent(1.0e-1) << " " << 1.0e-1 << endl;
  cout << size_of_exponent(1.0e1) << " " << 1.0e1 << endl;
  cout << size_of_exponent(1.0e11) << " " << 1.0e11 << endl;
  cout << size_of_exponent(1.0e111) << " " << 1.0e111 << endl;

  t.report();
  return 0;
}

