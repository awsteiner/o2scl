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

typedef std::function<int(size_t,const ubvector &,double &,
			     std::array<double,2> &)> point_funct;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  std::array<double,2> &)> fill_funct;

int point(size_t nv, const ubvector &pars, double &weight,
	  std::array<double,2> &dat) {
  dat[0]=pars[0]*pars[0];
  dat[1]=pars[0]*pars[1];
  weight=exp(-((pars[0]-0.2)*(pars[0]-0.2)+
	       (pars[1]-0.5)*(pars[1]-0.5))/2.0);
  return 0;
}
  
int fill_line(const ubvector &pars, double weight, 
	      std::vector<double> &line, std::array<double,2> &dat) {
  line.push_back(dat[0]);
  line.push_back(dat[1]);
  return 0;
}

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(2);

  point_funct pf=point;
  fill_funct ff=fill_line;
    
  mcmc_table<point_funct,fill_funct,std::array<double,2>,ubvector> mct;
  mct.user_seed=1;

  // Table columns
  vector<string> pnames={"x0","x1","x0x0","x0x1"};
  vector<string> punits={"MeV","MeV","MeV^2","MeV^2"};
  mct.set_names_units(pnames,punits);

  // Parameter limits
  ubvector low(2), high(2), init(2);
  low[0]=0.0;
  low[1]=0.0;
  high[0]=1.0;
  high[1]=1.0;
  init[0]=0.5;
  init[1]=0.5;
  
  shared_ptr<table_units<> > t=mct.get_table();

  // MCMC with a random walk of a fixed length
  mct.step_fac=2.0;
  mct.mcmc(10000,2,init,low,high,pf,ff);

  cout << "n_accept, n_reject, table lines: "
       << mct.n_accept << " " << mct.n_reject << " "
       << t->get_nlines() << endl;
  cout.precision(4);
  for(size_t i=0;i<t->get_nlines();i+=t->get_nlines()/10) {
    cout.width(4);
    cout << i << " " << t->get("mult",i) << " "
	 << t->get("weight",i) << " " << t->get("x0",i) << " ";
    cout << t->get("x1",i) << " " << t->get("x0x0",i) << " "
	 << t->get("x0x1",i) << endl;
  }
  cout.precision(6);
  cout << endl;
  
  // MCMC with affine-invariant sampling
  mct.aff_inv=true;
  mct.n_walk=10;
  mct.step_fac=2.0;

  mct.mcmc(10000,2,init,low,high,pf,ff);

  cout << "n_accept, n_reject, table lines: "
       << mct.n_accept << " " << mct.n_reject << " "
       << t->get_nlines() << endl;
  cout.precision(4);
  for(size_t i=0;i<t->get_nlines();i+=t->get_nlines()/10) {
    cout.width(4);
    cout << i << " " << t->get("mult",i) << " "
	 << t->get("weight",i) << " " << t->get("x0",i) << " ";
    cout << t->get("x1",i) << " " << t->get("x0x0",i) << " "
	 << t->get("x0x1",i) << endl;
  }
  cout.precision(6);

  tm.report();
  
  return 0;
}
