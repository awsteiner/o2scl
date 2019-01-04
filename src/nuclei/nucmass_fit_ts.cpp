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
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  
  cout.setf(ios::scientific);

  double res;
  nucmass_semi_empirical sem;
  nucmass_fit mf;

  nucmass_ame ame;
  o2scl_hdf::ame_load_ext(ame,"../../data/o2scl/nucmass/ame12.o2",
			  "ame12.o2");
  nucmass_ame amex;
  o2scl_hdf::ame_load_ext(amex,"../../data/o2scl/nucmass/ame12.o2",
			  "ame12.o2",true);
  nucmass_mnmsk_exp mexp;
  o2scl_hdf::mnmsk_load(mexp,"../../data/o2scl/nucmass/");
  nucmass_mnmsk mm;
  o2scl_hdf::mnmsk_load(mm,"../../data/o2scl/nucmass/");

  nucdist_set(mf.dist,ame);
  mf.fit(sem,res);
  cout << sem.B << " " << sem.Sv << " " << sem.Ss << " " 
       << sem.Ec << " " << sem.Epair << endl;
  cout << res << endl;
  t.test_gen(res<4.0,"Successful fit.");
  
  ubvector unc(1);
  unc[0]=3.0;
  mf.set_uncerts(unc);
  mf.fit_method=nucmass_fit::chi_squared_me;
  mf.fit(sem,res);
  cout << sem.B << " " << sem.Sv << " " << sem.Ss << " " 
       << sem.Ec << " " << sem.Epair << endl;

  /*
    int max_iso=30;
    ubvector qual, qual2;
    ubvector_int n_qual, n_qual2;
    
    mf.eval_isospin(sem,n_qual,qual);
    mf.eval_isospin(mm,n_qual2,qual2);
    for(size_t i=0;i<qual.size();i++) {
    cout << i << " " << n_qual[i] << " " << qual[i] << " ";
    cout << n_qual2[i] << " " << qual2[i] << endl;
    }
  */

  mf.fit_method=nucmass_fit::rms_mass_excess;
  nucdist_set(mf.dist,mexp);
  mf.eval(mm,res);
  cout << res << endl;
  t.test_rel(res,0.6806,1.0e-4,"Moller fit 1");

  nucdist_set(mf.dist,amex);
  mf.eval(mm,res);
  cout << res << endl;
  t.test_rel(res,0.6540311,1.0e-4,"Moller fit 2");

  nucdist_set(mf.dist,ame);
  mf.eval(mm,res);
  cout << res << endl;
  t.test_rel(res,0.894578,1.0e-4,"Moller fit 3");

  t.report();
  return 0;
}
  
