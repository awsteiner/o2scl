/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016, Andrew W. Steiner
  
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
#include <o2scl/mcmc.h>
#include <o2scl/vec_stats.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef std::function<int(const ubvector &,double,
			  size_t,bool,std::array<double,1> &)> measure_funct;

typedef std::function<int(size_t,const ubvector &,double &,
			     std::array<double,1> &)> point_funct;

std::vector<double> arr_x;
std::vector<double> arr_x2;
mcmc_table<point_funct,measure_funct,std::array<double,1>,ubvector> mct;

int point(size_t nv, const ubvector &pars, double &ret,
	  std::array<double,1> &dat) {
  dat[0]=pars[0]*pars[0];
  ret=exp(-pars[0]*pars[0]/2.0);
  return o2scl::success;
}

int meas(const ubvector &pars, double weight, size_t ix, bool new_meas,
	 std::array<double,1> &dat) {
  arr_x.push_back(pars[0]);
  arr_x2.push_back(dat[0]);
  if ((pars[0]*pars[0]-dat[0])>1.0e-10) {
    cerr << "Failure." << endl;
    exit(-1);
  }
  if (arr_x.size()==100) {
    return mcmc_base<point_funct,measure_funct,int,ubvector>::mcmc_done;
  }
  return 0;
}

int meas2(const ubvector &pars, double weight, size_t ix, bool new_meas,
	  std::array<double,1> &dat) {
  if ((pars[0]*pars[0]-dat[0])>1.0e-10) {
    cerr << "Failure 2." << endl;
    exit(-1);
  }
  mct.add_line(pars,weight,ix,new_meas,dat);
  if (mct.get_table()->get_nlines()>=100) {
    return mcmc_base<point_funct,measure_funct,int,ubvector>::mcmc_done;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(2);

  point_funct mf=point;
  measure_funct mf2=meas;
  measure_funct mf3=meas2;
    
  mcmc_base<point_funct,measure_funct,std::array<double,1>,ubvector> mc;
  ubvector init(1);
  ubvector low(1);
  ubvector high(1);
  init[0]=-0.01;
  low[0]=-5.0;
  high[0]=5.0;
  mc.verbose=1;
  mc.user_seed=1;

  mc.mcmc(1,init,low,high,mf,mf2);

  cout << vector_mean(arr_x) << endl;
  cout << vector_stddev(arr_x) << endl;
  cout << vector_mean(arr_x2) << endl;
  cout << vector_stddev(arr_x2) << endl;

  arr_x.clear();
  arr_x2.clear();
  mc.aff_inv=true;
  mc.n_walk=10;
  mc.step_fac=2.0;

  mc.mcmc(1,init,low,high,mf,mf2);

  cout << vector_mean(arr_x) << endl;
  cout << vector_stddev(arr_x) << endl;
  cout << vector_mean(arr_x2) << endl;
  cout << vector_stddev(arr_x2) << endl;

  arr_x.clear();
  arr_x2.clear();
  ubmatrix covar(1,1);
  covar(0,0)=1.0;
  prob_cond_mdim_gaussian<ubvector> pdmg(1,covar);

  mc.set_proposal(pdmg);
  mc.mcmc(1,init,low,high,mf,mf2);

  cout << vector_mean(arr_x) << endl;
  cout << vector_stddev(arr_x) << endl;
  cout << vector_mean(arr_x2) << endl;
  cout << vector_stddev(arr_x2) << endl;

  arr_x.clear();
  arr_x2.clear();
  mct.verbose=1;
  mct.user_seed=1;

  vector<string> pnames={"x"};
  vector<string> punits={"MeV"};
  mct.set_names_units(pnames,punits);

  mct.mcmc(1,init,low,high,mf,mf3);

  cout << vector_mean(mct.get_table()->get_column("param_x")) << endl;

  shared_ptr<table_units<> > t=mct.get_table();
  cout << "n_accept, n_reject, table lines: "
       << mct.n_accept << " " << mct.n_reject << " "
       << t->get_nlines() << endl;
  for(size_t i=0;i<t->get_nlines();i+=t->get_nlines()/10) {
    cout << i << " " << t->get("mult",i) << " "
	 << t->get("weight",i) << " " << t->get("param_x",i) << endl;
  }

  tm.report();
  
  return 0;
}
