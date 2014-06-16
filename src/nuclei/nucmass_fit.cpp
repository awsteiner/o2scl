/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

#include <o2scl/nucmass_fit.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_fit::nucmass_fit() {
  exp=&def_dist;
  mm=&def_mmin;
  def_mmin.ntrial*=10;
  even_even=false;
  minZ=8;
  minN=8;
  fit_method=rms_mass_excess;
  uncs.resize(1);
  uncs[0]=1.0;
}

void nucmass_fit::set_exp_mass(nucmass &nm, int maxA, 
			       bool include_neutron) {
  def_dist.set_dist(nm,maxA,include_neutron);
  return;
}

double nucmass_fit::min_fun(size_t nv, const ubvector &x) {

  double y=0.0;
  
  nmf->fit_fun(nv,x);
  
  eval(*nmf,y);
  
  return y;
}

void nucmass_fit::fit(nucmass_fit_base &n, double &fmin) {
  
  if (exp->size()==0) {
    O2SCL_ERR("No experimental masses to fit to in nucmass_fit::fit().",
	      exc_efailed);
  }

  nmf=&n;
  size_t nv=nmf->nfit;
  ubvector mx(nv);
  nmf->guess_fun(nv,mx);
  
  multi_funct11 mfm=
    std::bind(std::mem_fn<double(size_t,const ubvector &)>
	      (&nucmass_fit::min_fun),
	      this,std::placeholders::_1,std::placeholders::_2);
  
  mm->mmin(nv,mx,fmin,mfm);
  fmin=mfm(nv,mx);

  return;
}

void nucmass_fit::eval(nucmass &n, double &fmin) {

  fmin=0.0;

  if (exp->size()==0) {
    O2SCL_ERR("No experimental masses to fit to in nucmass_fit::eval().",
	      exc_efailed);
  }

  if (fit_method==rms_mass_excess) {

    size_t nn=0;
    for(nucdist::iterator ndi=exp->begin();ndi!=exp->end();ndi++) {
      int Z=ndi->Z;
      int N=ndi->N;
      if (N>=minN && Z>=minZ && (even_even==false || (N%2==0 && Z%2==0))) {
	fmin+=pow(ndi->mex*hc_mev_fm-n.mass_excess(Z,N),2.0);
	if (!o2scl::is_finite(fmin)) {
	  std::string s=((std::string)"Non-finite value for nucleus with Z=")+
	    itos(Z)+" and N="+itos(N)+" in nucmass_fit::eval() (1).";
	  O2SCL_ERR(s.c_str(),exc_efailed);
	}
	nn++;
      }
    }
    fmin=sqrt(fmin/nn);

  } else if (fit_method==rms_binding_energy) {

    size_t nn=0;
    for(nucdist::iterator ndi=exp->begin();ndi!=exp->end();ndi++) {
      int Z=ndi->Z;
      int N=ndi->N;
      if (N>=minN && Z>=minZ && (even_even==false || (N%2==0 && Z%2==0))) {
	fmin+=pow(ndi->be*hc_mev_fm-n.binding_energy(Z,N),2.0);
	if (!o2scl::is_finite(fmin)) {
	  std::string s=((std::string)"Non-finite value for nucleus with Z=")+
	    itos(Z)+" and N="+itos(N)+" in nucmass_fit::eval() (2).";
	  O2SCL_ERR(s.c_str(),exc_efailed);
	}
	nn++;
      }
    }
    fmin=sqrt(fmin/nn);
    
  } else if (fit_method==chi_squared_me) {

    size_t unc_ix=0;
    for(nucdist::iterator ndi=exp->begin();ndi!=exp->end();ndi++) {
      int Z=ndi->Z;
      int N=ndi->N;
      if (N>=minN && Z>=minZ && (even_even==false || (N%2==0 && Z%2==0))) {
	if (unc_ix>=uncs.size()) unc_ix=0;
	fmin+=pow((ndi->mex*hc_mev_fm-n.mass_excess(Z,N))/
		  (uncs[unc_ix]),2.0);
	if (!o2scl::is_finite(fmin)) {
	  std::string s=((std::string)"Non-finite value for nucleus with Z=")+
	    itos(Z)+" and N="+itos(N)+" in nucmass_fit::eval() (3).";
	  O2SCL_ERR(s.c_str(),exc_efailed);
	}
	unc_ix++;
      }
    }

  } else if (fit_method==chi_squared_be) {

    size_t unc_ix=0;
    for(nucdist::iterator ndi=exp->begin();ndi!=exp->end();ndi++) {
      int Z=ndi->Z;
      int N=ndi->N;
      if (N>=minN && Z>=minZ && (even_even==false || (N%2==0 && Z%2==0))) {
	if (unc_ix>=uncs.size()) unc_ix=0;
	fmin+=pow((ndi->be*hc_mev_fm-n.binding_energy(Z,N))/
		  (uncs[unc_ix]),2.0);
	if (!o2scl::is_finite(fmin)) {
	  std::string s=((std::string)"Non-finite value for nucleus with Z=")+
	    itos(Z)+" and N="+itos(N)+" in nucmass_fit::eval() (4).";
	  O2SCL_ERR(s.c_str(),exc_efailed);
	}
	unc_ix++;
      }
    }

  } else {
    O2SCL_ERR("Unknown fit method in nucmass_fit::eval().",exc_einval);
  }
    
  return;
}

void nucmass_fit::eval_isospin_beta(nucmass &n, ubvector_int &n_qual,
				    ubvector &qual, int max_iso) {
  
  ubvector_size_t cnt(200);
  ubvector_int opt(200);
  ubvector min(200);
  for(size_t i=0;i<200;i++) cnt[i]=0;

  // First pass, count and find minimum for each isotopic chain
  for(nucdist::iterator ndi=exp->begin();ndi!=exp->end();ndi++) {
    int Z=ndi->Z;
    int N=ndi->N;
    if (Z>=minZ) {
      if (cnt[Z]==0) {
	min[Z]=ndi->be/(N+Z);
	opt[Z]=N;
      } else {
	if ((ndi->be/(N+Z))<min[Z]) {
	  min[Z]=ndi->be/(N+Z);
	  opt[Z]=N;
	}
      }
      cnt[Z]++;
    }
  }

  // Second pass, compute quality for each isospin number
  qual.resize(2*max_iso+1);
  n_qual.resize(2*max_iso+1);
  for(int i=0;i<2*max_iso+1;i++) {
    n_qual[i]=0;
    qual[i]=0.0;
  }

  for(nucdist::iterator ndi=exp->begin();ndi!=exp->end();ndi++) {
    int Z=ndi->Z;
    int N=ndi->N;
    if (Z>=minZ && n.is_included(Z,N)) {
      int dev=N-opt[Z];
      if (abs(dev)<=max_iso) {
	int index=dev+max_iso;
	n_qual[index]++;
	qual[index]+=pow(ndi->mex*hc_mev_fm-n.mass_excess(Z,N),2.0);
      }
    }
  }
  
  for(int i=0;i<2*max_iso+1;i++) {
    if (n_qual[i]>0) qual[i]=sqrt(qual[i]/n_qual[i]);
  }

  return;
}

void nucmass_fit::eval_isospin(nucmass &n, ubvector_int &n_qual,
			       ubvector &qual, int min_iso, int max_iso) {
  
  if (max_iso<min_iso) {
    O2SCL_ERR("Max must be less than min in eval_isospin().",exc_einval);
  }
  size_t nv=max_iso-min_iso+1;

  qual.resize(nv);
  n_qual.resize(nv);
  for(size_t i=0;i<nv;i++) {
    n_qual[i]=0;
    qual[i]=0.0;
  }

  // First pass, count and find minimum for each isotopic chain
  for(nucdist::iterator ndi=exp->begin();ndi!=exp->end();ndi++) {
    int Z=ndi->Z;
    int N=ndi->N;
    int A=N-Z;
    if (Z>=8 && N>=8) {
      if (N-Z>=min_iso && N-Z<=max_iso) {
	int ix=N-Z-min_iso;
	n_qual[ix]++;
	qual[ix]+=pow(ndi->mex*hc_mev_fm-n.mass_excess(Z,N),2.0);
      }
    }
  }
  
  for(int i=0;i<((int)nv);i++) {
    if (n_qual[i]>0) qual[i]=sqrt(qual[i]/n_qual[i]);
  }

  return;
}

