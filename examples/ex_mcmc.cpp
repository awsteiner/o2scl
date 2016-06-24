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

typedef mcmc_table<point_funct,measure_funct,std::array<double,2>,
		   ubvector> child_t;

class ex_mcmc : public child_t {
  
public:
  
  virtual void fill_line(const vec_t &pars, double weight, 
			 std::vector<double> &line, data_t &dat) {
    child_t::fill_line(pars,weight,line,dat);
    line.push_back(dat[0]);
    line.push_back(dat[1]);
  }

  double point(size_t nv, const ubvector &pars, std::array<double,2> &dat) {
    dat[0]=pars[0]*pars[0];
    dat[1]=pars[0]*pars[1];
    return exp(-((pars[0]-0.2)*(pars[0]-0.2)+
		 (pars[1]-0.5)*(pars[1]-0.5))/2.0);
  }
  
  int meas(const ubvector &pars, double weight, size_t ix, bool new_meas,
	   std::array<double,2> &dat) {
    mct.add_line(pars,weight,ix,new_meas,dat);
    if (mct.get_table()->get_nlines()>=10000) {
      return mcmc_base<point_funct,measure_funct,int,ubvector>::mcmc_done;
    }
    return 0;
  }

};
  
int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(2);

  point_funct pf=point;
  measure_funct mf=meas;
    
  mct.verbose=1;
  mct.user_seed=1;

  vector<string> pnames={"x"};
  vector<string> punits={"MeV"};
  mct.set_names_units(pnames,punits);
  ubvector low(2)={0.0,0.0};
  ubvector high(2)={1.0,1.0};

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
