/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#include <o2scl/nstar_cold.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(void) {

  cout.setf(ios::scientific);
  test_mgr t;
  t.set_output_level(1);

  eos_had_rmf_hyp re;
  re.cs=sqrt(9.927);
  re.cw=sqrt(4.820);
  re.cr=sqrt(4.791);
  re.b=0.008659;
  re.c=-0.00241;
  re.mnuc=(re.def_neutron.m+re.def_proton.m)/2.0;

  cout << "GM91: " << endl;
  re.saturation();
  cout << "  Saturation density: " << re.n0 << endl;
  t.test_rel(re.n0,0.153,1.0e-3,"sat density");
  cout << "  Effective mass: " << re.msom << " "
       << re.def_neutron.ms/re.def_neutron.m << endl;
  t.test_rel(re.msom,re.def_neutron.ms/re.def_neutron.m,1.0e-6,"msom");
  t.test_rel(re.msom,0.78,1.0e-2,"msom 2");
  cout << "  Symmetry energy: " << re.esym*hc_mev_fm << endl;
  t.test_rel(re.esym*hc_mev_fm,32.5,0.1,"esym");
  cout << "  Zero pressure: " << re.def_thermo.pr << endl;
  t.test_rel(re.def_thermo.pr,0.0,1.0e-8,"zero press");
  cout << "  Energy per baryon: " << re.eoa*hc_mev_fm << " " 
       << (re.def_thermo.ed/re.n0-re.mnuc)*hc_mev_fm << endl;
  t.test_rel(re.eoa,re.def_thermo.ed/re.n0-re.mnuc,1.0e-6,"eoa");
  t.test_rel(re.eoa,-16.3/hc_mev_fm,0.1/hc_mev_fm,"eoa");
  cout << "  Thermomodynamic identity: " 
       << re.def_thermo.ed+re.def_thermo.pr-
    re.def_neutron.n*re.def_neutron.mu-
    re.def_proton.n*re.def_proton.mu << endl;
  t.test_rel(re.def_thermo.ed+re.def_thermo.pr-
	     re.def_neutron.n*re.def_neutron.mu-
	     re.def_proton.n*re.def_proton.mu,0.0,1.0e-9,"TI");
  cout << endl;

  // This corrects an apparent typo, 568 -> 658 in Table I of GM91
  vector<double> xw_dat={0.091,0.233,0.375,0.517,0.658,0.800,
			 0.942,1.08,1.23};
  size_t j=0;
  for(re.xs=0.2;re.xs<1.0001;re.xs+=0.1,j++) {
    re.calc_xw(-28.0/hc_mev_fm);
    cout << "  " << j << " " << re.xs << " " << re.xw << endl;
    if (j<7) {
      t.test_rel(re.xw,xw_dat[j],0.002,"xw");
    } else {
      t.test_rel(re.xw,xw_dat[j],0.01,"xw");
    }
  }
  cout << endl;

  re.inc_cascade=false;
  nstar_cold nc;
  nc.set_eos(re);
  //nc.verbose=1;
  //nc.calc_eos();
  //nc.calc_nstar();
  
  t.report();

  return 0;
}


