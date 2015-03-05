/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#include <o2scl/eos_nse.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_nse::eos_nse() {
  root=&def_root;
  err_nonconv=true;
}

void eos_nse::calc_mu(double mun, double mup, double T,
		      double &nB, double &Ye, thermo &th, 
		      vector<nucleus> &nd) {
  
  nB=0.0;
  Ye=0.0;

  for(size_t i=0;i<nd.size();i++) {
    nd[i].mu=mun*nd[i].N+mup*nd[i].Z-nd[i].be;
    nd[i].non_interacting=true;
    nd[i].inc_rest_mass=false;
    cla.calc_mu(nd[i],T);
    nB+=nd[i].n*((double)nd[i].A);
    Ye+=nd[i].n*((double)nd[i].Z);
    th.ed+=nd[i].ed;
    th.pr+=nd[i].pr;
    th.en+=nd[i].en;
  }

  Ye/=nB;

  return;
}

int eos_nse::calc_density(double nB, double Ye, double T, 
			  double &mun, double &mup, thermo &th, 
			  vector<nucleus> &nd) {

  ubvector x(2);
  x[0]=mun/T;
  x[1]=mup/T;

  mm_funct11 mfm=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
		     double,double,vector<nucleus> &)>
    (&eos_nse::solve_fun),this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,nB,Ye,T,std::ref(nd));

  int ret=root->msolve(2,x,mfm);
  if (ret!=0) {
    O2SCL_CONV_RET("Solver failed in eos_nse::calc_density().",
		   exc_efailed,err_nonconv);
  }

  mun=x[0]*T;
  mup=x[1]*T;
  
  // Final evaluation given new chemical potentials
  calc_mu(mun,mup,T,nB,Ye,th,nd);

  return success;
}

int eos_nse::solve_fun(size_t nv, const ubvector &x, ubvector &y, 
		       double nB, double Ye, double T,
		       vector<nucleus> &nd) {

  double mun=x[0]*T;
  double mup=x[1]*T;

  double nB2, Ye2;
  thermo th;
  
  calc_mu(mun,mup,T,nB2,Ye2,th,nd);

  y[0]=(nB2-nB)/nB;
  y[1]=Ye2-Ye;
  
  if (nB2<=0.0 || y[0]==-1.0) {
    return exc_ebadfunc;
  }
  if (!std::isfinite(y[0]) || !std::isfinite(y[1])) {
    O2SCL_ERR("Function values not finite in eos_nse::solve_fun().",
	      exc_efailed);
  }

  return success;
}

