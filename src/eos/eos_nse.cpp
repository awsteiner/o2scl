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
  mroot_ptr=&def_mroot;
  err_nonconv=true;
  verbose=0;
  def_mroot.ntrial=1000;
  make_guess_iters=40;
}

void eos_nse::calc_mu(double mun, double mup, double T,
		      double &nn, double &np, thermo &th, 
		      vector<nucleus> &nd) {

  nn=0.0;
  np=0.0;

  for(size_t i=0;i<nd.size();i++) {
    nd[i].mu=mun*nd[i].N+mup*nd[i].Z-nd[i].be;
    nd[i].non_interacting=true;
    nd[i].inc_rest_mass=false;
    cla.calc_mu(nd[i],T);
    nn+=nd[i].n*((double)nd[i].N);
    np+=nd[i].n*((double)nd[i].Z);
    th.ed+=nd[i].ed;
    th.pr+=nd[i].pr;
    th.en+=nd[i].en;
  }

  return;
}

int eos_nse::calc_density(double nn, double np, double T, 
			  double &mun, double &mup, thermo &th, 
			  vector<nucleus> &nd) {
  
  int ret=make_guess(mun,mup,T,th,nd,nn*1.0e-4,nn*1.0e4,
		     np*1.0e-4,np*1.0e4);
  if (ret!=0) {
    O2SCL_CONV_RET("Function make_guess() failed in eos_nse::calc_density().",
		   exc_efailed,err_nonconv);
  }
  if (verbose>0) {
    cout << "calc_density(), nn, np, mun, mup: " << nn << " " << np << " "
	 << mun << " " << mup << endl;
  }
  
  ubvector x(2);
  x[0]=mun/T;
  x[1]=mup/T;

  mm_funct11 mfm=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
		     double,double,vector<nucleus> &)>
    (&eos_nse::solve_fun),this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,nn,np,T,std::ref(nd));

  ret=mroot_ptr->msolve(2,x,mfm);
  if (ret!=0) {
    O2SCL_CONV_RET("Solver failed in eos_nse::calc_density().",
		   exc_efailed,err_nonconv);
  }

  mun=x[0]*T;
  mup=x[1]*T;
  
  // Final evaluation given new chemical potentials
  calc_mu(mun,mup,T,nn,np,th,nd);

  return success;
}

int eos_nse::solve_fun(size_t nv, const ubvector &x, ubvector &y, 
		       double nn, double np, double T,
		       vector<nucleus> &nd) {

  double mun=x[0]*T;
  double mup=x[1]*T;

  double nn2, np2;
  thermo th;
  
  calc_mu(mun,mup,T,nn2,np2,th,nd);

  y[0]=(nn2-nn)/nn;
  y[1]=(np2-np)/np;
  
  if (nn2<=0.0 || np2<=0.0 || std::isinf(nn2) || std::isinf(np2)) {
    return exc_ebadfunc;
  }

  return success;
}

int eos_nse::make_guess(double &mun, double &mup, double T,
			 thermo &th, std::vector<nucleus> &nd,
			 double nn_min, double nn_max,
			 double np_min, double np_max) {
  
  double nn, np;
  
  // Initial result
  calc_mu(mun,mup,T,nn,np,th,nd);
  if (verbose>0) {
    cout << mun << " " << mup << " " << nn << " " << np << endl;
  }

  // If we're already done, return
  if (std::isfinite(nn) && std::isfinite(np) &&
      nn>0.0 && np>0.0) {
    return o2scl::success;
  }
  
  double mun_step=T, mup_step=T, nn2, np2, mun2, mup2;
  
  // If the densities are infinite, or larger than 10^8, decrease the
  // chemical potentials until they are not.
  if (std::isinf(nn) || std::isinf(np) || nn>nn_max || np>np_max) {
    if (verbose>0) {
      cout << mun << " " << mup << " " << nn << " " << np
	   << " Infinite." << endl;
    }

    bool done=false;
    size_t k=0;
    while (done==false && k<make_guess_iters) {
      mun2=mun-mun_step;
      mup2=mup-mup_step;
      calc_mu(mun2,mup2,T,nn2,np2,th,nd);
      if (verbose>0) {
	cout << mun2 << " " << mup2 << " " << nn2 << " " << np2
	     << " k=" << k << endl;
      }
      if (!std::isinf(nn2) && !std::isinf(np2) && nn2<nn_max && np2<np_max) {
	done=true;
	mun=mun2;
	mup=mup2;
      } else {
	mun_step*=2.0;
	mup_step*=2.0;
      }
      k++;
    }
    if (done==false) {
      O2SCL_CONV2_RET("Failed to make densities small enough ",
		      "in eos_nse::make_guess().",exc_efailed,err_nonconv);
    }
    calc_mu(mun,mup,T,nn,np,th,nd);
    if (verbose>0) {
      cout << mun << " " << mup << " " << nn << " " << np
	   << " Done." << endl;
    }
  }

  // Now, if one of the densities are zero, then increase the
  // chemical potentials until they are nonzero

  mun_step=T;
  mup_step=T;
  if (nn==0.0 || np==0.0) {
    if (nn==0.0) mun_step=1.0e5*T;
    if (np==0.0) mup_step=1.0e5*T;
    
    bool done=false;
    size_t k=0;
    while (done==false && k<make_guess_iters) {
      mun2=mun+mun_step;
      mup2=mup+mup_step;
      calc_mu(mun2,mup2,T,nn2,np2,th,nd);
      if (verbose>0) {
	cout << mun2 << " " << mup2 << " " << nn2 << " " << np2
	     << " k=" << k << endl;
      }

      done=true;
      // If the new point gives infinite densities, decrease the step
      // size by a factor of two.
      if (std::isinf(nn2) || nn2>nn_max) {
	mun_step/=2.0;
	done=false;
	if (verbose>0) {
	  cout << "nn too large." << endl;
	}
      }
      if (std::isinf(np2) || np2>np_max) {
	mup_step/=2.0;
	done=false;
	if (verbose>0) {
	  cout << "np too large." << endl;
	}
      }
      if (done) {
	mun=mun2;
	mup=mup2;
	nn=nn2;
	np=np2;
	// Otherwise if one of the densities is still zero, update the
	// point and increase the step size by a factor of 3/2.
	if (nn2<=nn_min || np2<=np_min) {
	  if (nn2<=nn_min) {
	    mun_step*=1.5;
	    if (verbose>0) {
	      cout << "nn too small." << endl;
	    }
	  }
	  if (np2<=np_min) {
	    mup_step*=1.5;
	    if (verbose>0) {
	      cout << "np too small." << endl;
	    }
	  }
	  done=false;
	  calc_mu(mun,mup,T,nn,np,th,nd);
	}
      }
      k++;
    }
    if (done==false) {
      O2SCL_CONV2_RET("Failed to get densities in specified range ",
		      "in eos_nse::make_guess().",exc_efailed,err_nonconv);
    }
    calc_mu(mun,mup,T,nn,np,th,nd);
    if (verbose>0) {
      cout << mun << " " << mup << " " << nn << " " << np
	   << " Done." << endl;
    }
  }

  return o2scl::success;
}

