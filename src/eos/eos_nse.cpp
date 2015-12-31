/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
  mmin_ptr=&def_mmin;
  err_nonconv=true;
  verbose=0;
  def_mroot.ntrial=1000;
  make_guess_iters=60;
  make_guess_init_step=1.0e5;
  def_mmin.ntrial=1000;
  def_mmin.err_nonconv=false;
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

  int mg_ret, ds_ret, dm_ret;
  bool err_nonconv_save;

  // First, use make_guess() to get close to the correct solution
  double fac=2.0;
  mg_ret=make_guess(mun,mup,T,th,nd,nn/fac,nn*fac,
		    np/fac,np*fac,false);
  if (verbose>0) {
    cout << "In calc_density(), initial make_guess() returns "
	 << mg_ret << endl;
  }
  
  if (mg_ret==success) {

    // If that worked, try using the solver directly, ignoring
    // convergence errors for now
    err_nonconv_save=mroot_ptr->err_nonconv;
    mroot_ptr->err_nonconv=false;
    ds_ret=direct_solve(nn,np,T,mun,mup,th,nd);
    mroot_ptr->err_nonconv=err_nonconv_save;
    if (verbose>0) {
      cout << "In calc_density(), initial direct_solve() returns "
	   << ds_ret << endl;
    }

    // If the solver found a solution, then return
    if (ds_ret==success) {
      return success;
    }
    
  } 

  // Try minimizing instead
  dm_ret=density_min(nn,np,T,mun,mup,th,nd);
  if (verbose>0) {
    cout << "In calc_density(), initial density_min() returns "
	 << dm_ret << endl;
  }

  // Now a final attempt using the solver, use variable err_nonconv to
  // decide how to handle convergence errors
  if (err_nonconv==false) {
    err_nonconv_save=mroot_ptr->err_nonconv;
    mroot_ptr->err_nonconv=false;
  } 
  ds_ret=direct_solve(nn,np,T,mun,mup,th,nd);
  if (err_nonconv==false) {
    mroot_ptr->err_nonconv=err_nonconv_save;
  }
  if (verbose>0) {
    cout << "In calc_density(), final direct_solve() returns "
	 << ds_ret << endl;
  }

  if (ds_ret==success) return success;
  return exc_efailed;
}

int eos_nse::direct_solve(double nn, double np, double T, 
			  double &mun, double &mup, thermo &th, 
			  vector<nucleus> &nd) {
  ubvector x(2);
  x[0]=mun/T;
  x[1]=mup/T;

  mm_funct11 mfm=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &,double,
		     double,double,vector<nucleus> &)>
     (&eos_nse::solve_fun),this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,nn,np,T,std::ref(nd));
  
  int ret=mroot_ptr->msolve(2,x,mfm);

  mun=x[0]*T;
  mup=x[1]*T;
  
  // Final evaluation given new chemical potentials
  calc_mu(mun,mup,T,nn,np,th,nd);

  return ret;
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
			o2scl::thermo &th, std::vector<o2scl::nucleus> &nd,
			double nn_min, double nn_max,
			double np_min, double np_max, bool err_on_fail) {
  
  double nn, np;
    
  // Initial result
  calc_mu(mun,mup,T,nn,np,th,nd);
  if (verbose>1) {
    cout << "In make_guess()." << endl;
    cout << mun << " " << mup << " " << nn << " " << np << endl;
  }

  // If we're already done, return
  if (std::isfinite(nn) && std::isfinite(np) &&
      nn>nn_min && np>np_min && nn<nn_max && np<np_max) {
    return o2scl::success;
  }

  bool mun_changing=true, mup_changing=true;
  double mun_step=0.0, mup_step=0.0, nn2, np2, mun2, mup2;

  if (std::isinf(nn) || nn>nn_max) {
    mun_step=-make_guess_init_step*T;
  } else if (nn<nn_min) {
    mun_step=make_guess_init_step*T;
  } else {
    mun_changing=false;
  }
  if (std::isinf(np) || np>np_max) {
    mup_step=-make_guess_init_step*T;
  } else if (np<np_min) {
    mup_step=make_guess_init_step*T;
  } else {
    mup_changing=false;
  }

  bool done=false;
  size_t k=0;
  while (done==false && k<make_guess_iters) {

    // Go to and evaluate new point
    if (mun_changing) {
      mun2=mun+mun_step;
    } else {
      mun2=mun;
    }
    if (mup_changing) {
      mup2=mup+mup_step;
    } else {
      mup2=mup;
    }
    calc_mu(mun2,mup2,T,nn2,np2,th,nd);
      
    if (verbose>1) {
      cout << "k=" << k << endl;
      cout.setf(ios::showpos);
      cout << "before   : " << mun << " " << mup << " ";
      cout.unsetf(ios::showpos);
      cout << nn << " " << np << endl;
      cout.setf(ios::showpos);
      cout << "step, min: " << mun_step << " " << mup_step << " ";
      cout.unsetf(ios::showpos);
      cout << nn_min << " " << np_min << endl;
      cout.setf(ios::showpos);
      cout << "after    : " << mun2 << " " << mup2 << " ";
      cout.unsetf(ios::showpos);
      cout << nn2 << " " << np2 << endl;
      cout.setf(ios::showpos);
      cout << "0,max    : " << 0.0 << " " << 0.0 << " ";
      cout.unsetf(ios::showpos);
      cout << nn_max << " " << np_max << endl;
    }

    bool accept=true;
    bool eval_protons=true;
    
    // Update based on neutron density
    if (mun_changing==false) {
      if (std::isinf(nn2) || nn2>nn_max || nn2<nn_min) {
	// The proton step put the neutron density out of range, so
	// shrink the proton step and try again
	mup_step/=1.5;
	accept=false;
	// No point in looking at protons
	eval_protons=false;
	if (verbose>1) {
	  cout << "Value of mun_changing false but neutrons "
	       << "out of range." << endl;
	}
      }
    } else if (mun_step<0.0) {
      // Neutron step is negative
      if (nn2<nn_min) {
	// If the step went too far, decrease the neutron step
	mun_step/=1.5;
	accept=false;
	if (verbose>1) {
	  cout << "Neutron step went too far (mun_step<0)." << endl;
	}
      } else if (std::isfinite(nn2) && nn2<nn_max) {
	// If the step is just right, stop changing it
	// and accept the step
	//mun_changing=false;
	mun_step/=2.0;
	if (verbose>1) {
	  cout << "Neutrons now in range (mun_step<0)." << endl;
	}
      }
      // Otherwise, the step didn't go far enough to we accept and try
      // another one.
    } else {
      // Neutron step is positive
      if (std::isinf(nn2) || nn2>nn_max) {
	// If the step went too far, decrease the neutron step
	mun_step/=1.5;
	accept=false;
	if (verbose>1) {
	  cout << "Neutron step went too far (mun_step>0)." << endl;
	}
      } else if (nn2>nn_min) {
	// If the step is just right, stop changing it
	// and accept the step
	//mun_changing=false;
	mun_step/=2.0;
	if (verbose>1) {
	  cout << "Neutrons now in range (mun_step>0)." << endl;
	}
      }
      // Otherwise, the step didn't go far enough to we accept and try
      // another one.
    }

    // Update based on proton density
    if (eval_protons) {
      if (mup_changing==false) {
	if (std::isinf(np2) || np2>np_max || np2<np_min) {
	  // The neutron step put the proton density out of range, so
	  // shrink the neutron step and try again
	  mun_step/=1.5;
	  accept=false;
	  if (verbose>1) {
	    cout << "Value of mup_changing false but protons "
		 << "out of range." << endl;
	  }
	}
      } else if (mup_step<0.0) {
	// Proton step is negative
	if (np2<np_min) {
	  // If the step went too far, decrease the proton step
	  mup_step/=1.5;
	  accept=false;
	  if (verbose>1) {
	    cout << "Proton step went too far (mup_step<0)." << endl;
	  }
	} else if (std::isfinite(np2) && np2<np_max) {
	  // If the step is just right, stop changing it
	  // and accept the step
	  //mup_changing=false;
	  mup_step/=2.0;
	  if (verbose>1) {
	    cout << "Protons now in range (mup_step<0)." << endl;
	  }
	}
	// Otherwise, the step didn't go far enough to we accept and try
	// another one.
      } else {
	// Proton step is positive
	if (std::isinf(np2) || np2>np_max) {
	  // If the step went too far, decrease the proton step
	  mup_step/=1.5;
	  accept=false;
	  if (verbose>1) {
	    cout << "Proton step went too far (mup_step>0)." << endl;
	  }
	} else if (np2>np_min) {
	  // If the step is just right, stop changing it
	  // and accept the step
	  //mup_changing=false;
	  mup_step/=2.0;
	  if (verbose>1) {
	    cout << "Protons now in range (mup_step>0)." << endl;
	  }
	}
	// Otherwise, the step didn't go far enough to we accept and try
	// another one.
      }
    }

    if (accept) {
      mun=mun2;
      mup=mup2;
      calc_mu(mun,mup,T,nn,np,th,nd);
      if (verbose>1) {
	cout << "Accept." << endl;
      }
      if (std::isfinite(nn) && std::isfinite(np) &&
	  nn>nn_min && np>np_min && nn<nn_max && np<np_max) {
	done=true;
	if (verbose>1) {
	  cout << "Done." << endl;
	}
      }
    }

    if (verbose>2) {
      char ch;
      cin >> ch;
    }

    // Go to next iteration
    k++;
  }

  if (done==false) {
    if (err_on_fail) {
      O2SCL_ERR("Failed in eos_nse::make_guess().",exc_efailed);
    } else {
      return exc_efailed;
    }
  }

  return o2scl::success;
}

int eos_nse::density_min(double nn, double np, double T, 
			 double &mun, double &mup, o2scl::thermo &th, 
			 std::vector<o2scl::nucleus> &nd) {

  o2scl::multi_funct11 mf=std::bind
    (std::mem_fn<double(size_t,const ubvector &, double, double,
			double, o2scl::thermo &,
			std::vector<o2scl::nucleus> &)>
     (&eos_nse::minimize_fun),this,std::placeholders::_1,
     std::placeholders::_2,T,nn,np,std::ref(th),std::ref(nd));
  
  ubvector x(2);
  x[0]=mun;
  x[1]=mup;
  
  double y=1.0;

  int ret;
  for(size_t i=0;i<5 && y>1.0e-4;i++) {
    ret=mmin_ptr->mmin(2,x,y,mf);
    if (verbose>0) {
      cout << "In density_min(), y=" << y << ", ret=" << ret << endl;
    }
     
  }

  mun=x[0];
  mup=x[1];
  
  return ret;
}
  
double eos_nse::minimize_fun(size_t nv, const ubvector &x, double T,
			     double nn, double np, o2scl::thermo &th,
			     std::vector<o2scl::nucleus> &nd) {
  double mun=x[0], mup=x[1], nn2, np2;
  calc_mu(mun,mup,T,nn2,np2,th,nd);
  if (std::isinf(nn2) || std::isinf(np2) || nn2>10.0 || np2>10.0) {
    return 1.0e100;
  }
  return pow((nn2-nn)/nn,2.0)+pow((np2-np)/np,2.0);
}
