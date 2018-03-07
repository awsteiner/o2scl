/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/test_mgr.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/eos_had_rmf_hyp.h>
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(void) {

  cout.setf(ios::scientific);
  t.set_output_level(2);

  test_mgr t;

  eos_had_rmf_hyp re;
  re.cs=9.927;
  re.cw=4.820;
  re.cr=4.791;
  re.b=0.008659;
  re.c=-0.00241;
  
  re.saturation();
  cout << "  Saturation density: " << re.n0 << endl;
  //t.test_rel(re.n0,0.148,1.0e-2,"sat density");
  cout << "  Effective mass: " << re.msom << " " << nferm.ms/nferm.m << endl;
  //t.test_rel(re.msom,nferm.ms/nferm.m,1.0e-6,"msom");
  cout << "  Zero pressure: " << th.pr << endl;
  //t.test_rel(th.pr,0.0,1.0e-8,"zero press");
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << " " 
       << (th.ed/re.n0-re.mnuc)*hc_mev_fm << endl;
  //t.test_rel(re.eoa,th.ed/re.n0-re.mnuc,1.0e-6,"eoa");
  cout << "  Thermodynamic identity: " 
       << th.ed+th.pr-nferm.n*nferm.mu-p.n*p.mu << endl;
  //t.test_rel(th.ed+th.pr-nferm.n*nferm.mu-p.n*p.mu,0.0,1.0e-9,"TI");
  cout << endl;

  t.report();

  return 0;
}


