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

/// The MCMC object
mcmc_table<point_funct,fill_funct,std::array<double,2>,ubvector> mct;

/** \brief The objective function for the MCMC

    Here, the variable 'weight' stores the objective function based on
    the parameters stored in \c pars. The object 'dat' stores any
    auxillary quantities which can be computed at every point in
    parameter space.
 */
int point(size_t nv, const ubvector &pars, double &weight,
	  std::array<double,2> &dat) {

  weight=exp(-((pars[0]-0.2)*(pars[0]-0.2)+
	       (pars[1]-0.5)*(pars[1]-0.5))/2.0);

  dat[0]=pars[0]*pars[0];
  dat[1]=pars[0]*pars[1];
  return 0;
}

/** \brief Add auxillary quantities to 'line' so they can be
    stored in the table
 */
int fill_line(const ubvector &pars, double weight, 
	      std::vector<double> &line, std::array<double,2> &dat) {
  line.push_back(dat[0]);
  line.push_back(dat[1]);
  if (mct.get_table()->get_nlines()==10000) {
    return mcmc_base<point_funct,fill_funct,int,ubvector>::mcmc_done;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(2);

  point_funct pf=point;
  fill_funct ff=fill_line;

  // Table column names and units. We must specify first the names for
  // the parameters first and then the names of the auxillary
  // parameters.
  vector<string> pnames={"x0","x1","x0x0","x0x1"};
  vector<string> punits={"MeV","MeV","MeV^2","MeV^2"};
  mct.set_names_units(pnames,punits);

  // Parameter limits and initial point
  ubvector low(2), high(2), init(2);
  low[0]=0.0;
  low[1]=0.0;
  high[0]=1.0;
  high[1]=1.0;
  init[0]=0.5;
  init[1]=0.5;

  // Get a pointer to the results table
  shared_ptr<table_units<> > t=mct.get_table();

  // MCMC with a random walk of a fixed length
  mct.step_fac=2.0;
  mct.mcmc(2,init,low,high,pf,ff);

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

  mct.mcmc(2,init,low,high,pf,ff);

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
