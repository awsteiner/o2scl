 /*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/fermion_rel2.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/test_mgr.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/eos_sn.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  fermion f;
  f.g=2.0;
  fermion_rel fr;
  fermion_rel2 fr2;
  
  for(double m=1.0e-2;m<1.01e2;m*=10.0) {
    f.m=m;
    for(double T=1.0e-2;T<1.01e2;T*=10.0) {
      for(double mu=1.0e-2;mu<1.01e2;mu*=10.0) {
        f.mu=mu;
        fr.calc_mu(f,T);
        cout << f.n << " " << f.en << " ";
        fr2.calc_mu(f,T);
        cout << f.n << " " << f.en << endl;
      }
    }
  }
  
  t.report();

  return 0;
}

