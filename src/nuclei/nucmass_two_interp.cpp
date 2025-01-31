/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2025, Andrew W. Steiner
  
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

#include <o2scl/nucmass_two_interp.h>
// For unit conversions
#include <o2scl/lib_settings.h>

#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;

nucmass_two_interp::nucmass_two_interp() {
  nfb=&def_fit;
  ib=&def_ib;
  nfit=def_fit.nfit;
}

nucmass_two_interp::~nucmass_two_interp() {
}

void nucmass_two_interp::set_default() {

  nucmass_ame ame;
  ame.load("20");

  nucmass_fit nf;
  nucdist_set(nf.dist,ame);
  double res;
  nf.def_mmin.verbose=2;
  nf.fit(*nfb,res);

  table<> tab;
  nf.eval_table(*nfb,res,true,tab);

  // Now use 'tab' to initialize interpm_idw

  const_matrix_view_table<> ix(tab,{"Z","N"});
  std::vector<std::string> out_cols={"me_th"};
  const_matrix_view_table<> iy(tab,out_cols);
  def_ib.verbose=2;
  ib->set_data(2,out_cols.size(),tab.get_nlines(),ix,iy);
  
  return;
}

double nucmass_two_interp::mass_excess_d(double Z, double N) {
  double ret=nfb->mass_excess_d(Z,N);
  ubvector x(2);
  ubvector y(1);
  ib->eval(x,y);
  ret+=y[0];
  return ret;
}

