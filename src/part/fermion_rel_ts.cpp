/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
#include <o2scl/fermion_rel.h>
#include <o2scl/test_mgr.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/eos_sn.h>

#ifdef O2SCL_SET_MULTIP
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  
#ifdef O2SCL_SET_MULTIP
  
  std::string arg;
  if (argc>=2) {
    arg=argv[1];
  }
  
  fermion f;
  f.g=2;

  fermion_rel fr;
  fr.verify_ti=true;
  fr.err_nonconv=false;
  fr.nit.err_nonconv=false;
  fr.dit.err_nonconv=false;
  fr.it_multip.err_nonconv=false;
  
  fermion_ld fld;
  fld.g=2;
  
  fermion_rel_ld frld;
  frld.verify_ti=true;
  frld.err_nonconv=false;
  frld.nit.err_nonconv=false;
  frld.dit.err_nonconv=false;
  frld.it_multip.err_nonconv=false;
  
  fermion_cdf25 f25;
  f25.g=2;
  
  fermion_rel_cdf25 fr25;
  fr25.verify_ti=true;
  fr25.err_nonconv=false;
  fr25.nit.err_nonconv=false;
  fr25.dit.err_nonconv=false;
  fr25.it_multip.err_nonconv=false;

  int first_test=0;
  int last_test=1e6;
  double test_shift=1.0;
  
  // An exhaustive comparison of the fermion_rel class at various
  // levels of precision

  if (arg=="1") {

    // This tests higher-accuracy settings. Smaller tolerances for the
    // integrals cause some of the integrations to fail. Note that a
    // larger value of 'test_shift' is also required, but it's not
    // clear to me if this really represents an accuracy issue.
    
    fr.dit.tol_abs=1.0e-11;
    fr.dit.tol_rel=1.0e-11;
    fr.nit.tol_abs=1.0e-11;
    fr.nit.tol_rel=1.0e-11;
    fr.upper_limit_fac=40.0;
    fr.density_root.tol_rel=1.0e-10;
    fr.def_massless_root.tol_rel=1.0e-10;
    
    frld.dit.tol_abs=1.0e-15;
    frld.dit.tol_rel=1.0e-15;
    frld.nit.tol_abs=1.0e-15;
    frld.nit.tol_rel=1.0e-15;
    frld.upper_limit_fac=60.0;
    frld.density_root.tol_rel=1.0e-14;
    frld.def_massless_root.tol_rel=1.0e-14;
    
    fr25.dit.tol_abs=1.0e-18;
    fr25.dit.tol_rel=1.0e-18;
    fr25.nit.tol_abs=1.0e-18;
    fr25.nit.tol_rel=1.0e-18;
    fr25.upper_limit_fac=80.0;
    fr25.density_root.tol_rel=1.0e-18;
    fr25.def_massless_root.tol_rel=1.0e-18;

    test_shift=10.0;
  }
  
  if (arg=="2") {

    last_test=0;
    
    // I think this runs without throwing any exceptions, but it needs
    // some work, especially to ensure the thermodynamic identity is
    // satisfied at higher precision. Either way, it is quite slow.
    
    fr.multip=true;
    fr.upper_limit_fac=40.0;
    fr.density_root.tol_rel=1.0e-10;
    fr.def_massless_root.tol_rel=1.0e-10;
    
    frld.multip=true;
    frld.upper_limit_fac=52.0;
    frld.deg_entropy_fac=52.0;
    frld.tol_expan=1.0e-18;
    frld.exp_limit=11400.0;
    frld.density_root.tol_rel=1.0e-18;
    frld.def_massless_root.tol_rel=1.0e-18;
    
    fr25.multip=true;
    fr25.upper_limit_fac=62.0;
    fr25.deg_entropy_fac=62.0;
    fr25.tol_expan=1.0e-23;
    fr25.exp_limit=6.7e7;
    fr25.density_root.tol_rel=1.0e-23;
    fr25.def_massless_root.tol_rel=1.0e-23;
    
  }

  /*
  if (arg=="3") {
    
    fr.dit.tol_abs=1.0e-13;
    fr.dit.tol_rel=1.0e-13;
    fr.nit.tol_abs=1.0e-13;
    fr.nit.tol_rel=1.0e-13;
    fr.upper_limit_fac=40.0;
    fr.density_root.tol_rel=1.0e-10;
    
    frld.dit.tol_abs=1.0e-18;
    frld.dit.tol_rel=1.0e-18;
    frld.nit.tol_abs=1.0e-18;
    frld.nit.tol_rel=1.0e-18;
    frld.density_root.tol_rel=1.0e-18;
    frld.def_massless_root.tol_rel=1.0e-18;
    frld.upper_limit_fac=52.0;
    frld.deg_entropy_fac=52.0;
    frld.tol_expan=1.0e-18;
    frld.exp_limit=11400.0;
    
    fr25.dit.tol_abs=1.0e-25;
    fr25.dit.tol_rel=1.0e-25;
    fr25.nit.tol_abs=1.0e-25;
    fr25.nit.tol_rel=1.0e-25;
    fr25.density_root.tol_rel=1.0e-23;
    fr25.def_massless_root.tol_rel=1.0e-23;
    fr25.upper_limit_fac=62.0;
    fr25.deg_entropy_fac=62.0;
    fr25.tol_expan=1.0e-23;
    fr25.exp_limit=6.7e7;
    
  }
  */
  
  int count=0;
  
  if (argc<2) {
    
    part_cal_new<> pcn;
    pcn.test_calc_mu(f,fld,f25,fr,frld,fr25,t,count,first_test,last_test,
                     1502,1479,2370,2134,1400,1971,2694);
    pcn.test_pair_mu(f,fld,f25,fr,frld,fr25,t,count,first_test,last_test,
                     1502,1479,2272,2132,1376,2006,2680);
    
    /*
    // Normal case, used in 'make o2scl-test'
    t.test_gen(cmu_n>=1543,"cmu_n");
    t.test_gen(cmu_en>=1544,"cmu_en");
    t.test_gen(cmu_ld_n>=2407,"cmu_ld_n");
    t.test_gen(cmu_ld_en>=2179,"cmu_ld_en");
    t.test_gen(cmu_ti>=1400,"cmu_ti");
    t.test_gen(cmu_ld_ti>=1906,"cmu_ld_ti");
    t.test_gen(cmu_25_ti>=2597,"cmu_25_ti");
    
    t.test_gen(pmu_n>=1582,"pmu_n");
    t.test_gen(pmu_en>=1604,"pmu_en");
    t.test_gen(pmu_ld_n>=2336,"pmu_ld_n");
    t.test_gen(pmu_ld_en>=2218,"pmu_ld_en");
    t.test_gen(pmu_ti>=1376,"pmu_ti");
    t.test_gen(pmu_ld_ti>=1626,"pmu_ld_ti");
    t.test_gen(pmu_25_ti>=1931,"pmu_25_ti");
    
    t.test_gen(cd_mu>=891,"cd_mu");
    t.test_gen(cd_en>=837,"cd_en");
    t.test_gen(cd_ld_mu>=1253,"cd_ld_mu");
    t.test_gen(cd_ld_en>=1117,"cd_ld_en");
    t.test_gen(cd_ti>=967,"cd_ti");
    t.test_gen(cd_ld_ti>=967,"cd_ld_ti");
    t.test_gen(cd_25_ti>=967,"cd_25_ti");
    
    t.test_gen(pd_mu>=850,"pd_mu");
    t.test_gen(pd_en>=814,"pd_en");
    t.test_gen(pd_ld_mu>=1190,"pd_ld_mu");
    t.test_gen(pd_ld_en>=1101,"pd_ld_en");
    t.test_gen(pd_ti>=676,"pd_ti");
    t.test_gen(pd_ld_ti>=817,"pd_ld_ti");
    t.test_gen(pd_25_ti>=1025,"pd_25_ti");
    */
    
  } else if (arg=="1") {

    part_cal_new<> pcn;
    pcn.test_calc_mu(f,fld,f25,fr,frld,fr25,t,count,first_test,last_test,
                     1738,1612,2371,2169,1540,1992,2702);
    
    /*
    // Improved accuracy case
    t.test_gen(cmu_n>=1738,"cmu_n");
    t.test_gen(cmu_en>=1612,"cmu_en");
    t.test_gen(cmu_ld_n>=2379,"cmu_ld_n");
    t.test_gen(cmu_ld_en>=2174,"cmu_ld_en");
    t.test_gen(cmu_ti>=1540,"cmu_ti");
    t.test_gen(cmu_ld_ti>=1992,"cmu_ld_ti");
    t.test_gen(cmu_25_ti>=2702,"cmu_25_ti");
    
    t.test_gen(pmu_n>=1756,"pmu_n");
    t.test_gen(pmu_en>=1670,"pmu_en");
    t.test_gen(pmu_ld_n>=2298,"pmu_ld_n");
    t.test_gen(pmu_ld_en>=2222,"pmu_ld_en");
    t.test_gen(pmu_ti>=1592,"pmu_ti");
    t.test_gen(pmu_ld_ti>=1821,"pmu_ld_ti");
    t.test_gen(pmu_25_ti>=2144,"pmu_25_ti");
    
    t.test_gen(cd_mu>=991,"cd_mu");
    t.test_gen(cd_en>=901,"cd_en");
    t.test_gen(cd_ld_mu>=1237,"cd_ld_mu");
    t.test_gen(cd_ld_en>=1104,"cd_ld_en");
    t.test_gen(cd_ti>=969,"cd_ti");
    t.test_gen(cd_ld_ti>=969,"cd_ld_ti");
    t.test_gen(cd_25_ti>=969,"cd_25_ti");
    
    t.test_gen(pd_mu>=939,"pd_mu");
    t.test_gen(pd_en>=888,"pd_en");
    t.test_gen(pd_ld_mu>=1190,"pd_ld_mu");
    t.test_gen(pd_ld_en>=1094,"pd_ld_en");
    t.test_gen(pd_ti>=819,"pd_ti");
    t.test_gen(pd_ld_ti>=997,"pd_ld_ti");
    t.test_gen(pd_25_ti>=1271,"pd_25_ti");
    */
    
  } else if (arg=="2") {

    part_cal_new<> pcn;
    pcn.test_calc_mu(f,fld,f25,fr,frld,fr25,t,count,first_test,last_test,
                     1826,1827,2420,2405,1754,2154,3311);
    
  }

#endif
  
  t.report();

  return 0;
}

