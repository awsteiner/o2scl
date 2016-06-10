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
			  size_t,bool,std::array<double,2> &)> measure_funct;

typedef std::function<double(size_t,const ubvector &,
			     std::array<double,2> &)> point_funct;

std::vector<double> arr_x;
std::vector<double> arr_x2;
mcmc_table<point_funct,measure_funct,std::array<double,2>,ubvector> mct;

double point(size_t nv, const ubvector &pars, std::array<double,2> &dat) {
  return exp(-pars[0]*pars[0]/2.0);
}

int meas(const ubvector &pars, double weight, size_t ix, bool new_meas,
	 std::array<double,2> &dat) {
  arr_x.push_back(pars[0]);
  arr_x2.push_back(pars[0]*pars[0]);
  if (arr_x.size()==100) return -1;
  return 0;
}

int meas2(const ubvector &pars, double weight, size_t ix, bool new_meas,
	  std::array<double,2> &dat) {
  mct.add_line(pars,weight,ix,new_meas);
  if (mct.get_table()->get_nlines()>=100) return -1;
  return 0;
}

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(2);

  point_funct mf=point;
  measure_funct mf2=meas;
  measure_funct mf3=meas2;
    
  mcmc_base<point_funct,measure_funct,std::array<double,2>,ubvector> mc;
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
  mc.nwalk=10;
  mc.step_fac=2.0;

  mc.mcmc(1,init,low,high,mf,mf2);

  cout << vector_mean(arr_x) << endl;
  cout << vector_stddev(arr_x) << endl;
  cout << vector_mean(arr_x2) << endl;
  cout << vector_stddev(arr_x2) << endl;

  arr_x.clear();
  arr_x2.clear();
  ubvector peak(1);
  peak[0]=0.0;
  ubmatrix covar(1,1);
  covar(0,0)=1.0;
  prob_dens_mdim_gaussian<ubvector> pdmg(1,peak,covar);

  mc.set_hastings(pdmg);
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

  cout << vector_mean(mct.get_table()->get_column("x")) << endl;

  shared_ptr<table_units<> > t=mct.get_table();
  for(size_t i=0;i<t->get_nlines();i+=t->get_nlines()/10) {
    cout << i << " " << t->get("mult",i) << " "
	 << t->get("weight",i) << " " << t->get("x",i) << endl;
  }

  tm.report();
  
  return 0;
}
