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

#include <o2scl/test_mgr.h>
#include <o2scl/boson.h>
#include <o2scl/boson_rel.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  boson b;
  b.g=2.0;
  boson_rel br;
  
  cout << "An alternative way of testing ndeg_terms(), by\n  "
       << "comparing with numerical values in bosons.ipynb." << endl;
      
  double tt, psi;
  b.ms=0.511/197.33;
  b.inc_rest_mass=true;
  tt=0.1/0.511;
  psi=(8.0e-12-b.ms)/0.1*197.33;
  for(size_t j=1;j<6;j++) {
    double pterm, nterm, enterm, edterm;
    br.ndeg_terms(j,tt,psi*tt,b.ms,b.inc_rest_mass,false,
                   pterm,nterm,enterm,edterm);
    std::cout.setf(std::ios::showpos);
    std::cout.precision(8);
    std::cout << j << " " << pterm << " " << nterm << " "
              << enterm << " " << edterm << std::endl;
  }
  std::cout << endl;
    
  for(size_t j=1;j<6;j++) {
    double pterm, nterm, enterm, edterm;
    br.ndeg_terms(j,tt,psi*tt,b.ms,b.inc_rest_mass,true,
                   pterm,nterm,enterm,edterm);
    std::cout.setf(std::ios::showpos);
    std::cout.precision(8);
    std::cout << j << " " << pterm << " " << nterm << " "
              << enterm << " " << edterm << std::endl;
  }
  
  t.report();
  return 0;
}
