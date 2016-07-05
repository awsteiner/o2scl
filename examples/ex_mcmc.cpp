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
/* Example: ex_mcmc.cpp
   -------------------------------------------------------------------

   An example which demonstrates the generation of an arbitrary
   distribution through Markov chain Monte Carlo.
*/
#include <o2scl/mcmc.h>
#include <o2scl/vec_stats.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_io.h>

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
  
  weight=-((pars[0]-0.2)*(pars[0]-0.2)+
	   (pars[1]-0.5)*(pars[1]-0.5));
  
  dat[0]=pars[0]*pars[0];
  dat[1]=pars[0]*pars[0]*pars[1]*pars[1];
  return 0;
}

/** \brief Add auxillary quantities to 'line' so they can be
    stored in the table
*/
int fill_line(const ubvector &pars, double weight, 
	      std::vector<double> &line, std::array<double,2> &dat) {
  line.push_back(dat[0]);
  line.push_back(dat[1]);
  if (mct.get_table()->get_nlines()==100000) {
    return mcmc_base<point_funct,fill_funct,int,ubvector>::mcmc_done;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(1);

  // Exact results
  double res[3]={0.0,0.511498,0.344514};
  
  point_funct pf=point;
  fill_funct ff=fill_line;

  // Table column names and units. We must specify first the names for
  // the parameters first and then the names of the auxillary
  // parameters.
  vector<string> pnames={"x0","x1","x0sq","x0sq_x1sq"};
  vector<string> punits={"MeV","MeV","MeV^2","MeV^4"};
  mct.set_names_units(pnames,punits);

  // Parameter limits and initial point
  ubvector low(2), high(2), init(2);
  low[0]=-2.0;
  low[1]=-2.0;
  high[0]=2.0;
  high[1]=2.0;
  init[0]=0.2;
  init[1]=0.5;

  // Get a pointer to the results table
  shared_ptr<table_units<> > t=mct.get_table();

  // MCMC with a random walk of a fixed length
  mct.step_fac=2.0;
  //mct.verbose=2;
  mct.mcmc(2,init,low,high,pf,ff);

  cout << "n_accept, n_reject, table lines: "
       << mct.n_accept << " " << mct.n_reject << " "
       << t->get_nlines() << endl;
  cout.precision(4);
  for(size_t i=0;i<t->get_nlines();i+=t->get_nlines()/10) {
    cout.width(5);
    cout << i << " " << t->get("mult",i) << " ";
    cout.setf(ios::showpos);
    cout << t->get("weight",i) << " " << t->get("x0",i) << " ";
    cout << t->get("x1",i) << " " << t->get("x0sq",i) << " "
	 << t->get("x0sq_x1sq",i) << endl;
    cout.unsetf(ios::showpos);
  }
  cout.precision(6);
  cout << endl;
  
  size_t n=t->get_nlines();

  double t_avg=wvector_mean(n,t->get_column("x0sq"),
			    t->get_column("mult"));
  double t_stddev=wvector_stddev(n,t->get_column("x0sq"),
				 t->get_column("mult"));
  double t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,res[1],t_avgerr*10.0,"tab 1");

  t_avg=wvector_mean(n,t->get_column("x0sq_x1sq"),
		     t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x0sq_x1sq"),
			  t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,res[2],t_avgerr*10.0,"tab 2");
  cout << endl;

  // MCMC with affine-invariant sampling
  mct.aff_inv=true;
  mct.n_walk=10;
  mct.step_fac=20.0;

  mct.mcmc(2,init,low,high,pf,ff);

  cout << "n_accept, n_reject, table lines: "
       << mct.n_accept << " " << mct.n_reject << " "
       << t->get_nlines() << endl;
  cout.precision(4);
  for(size_t i=0;i<t->get_nlines();i+=t->get_nlines()/10) {
    cout.width(5);
    cout << i << " " << t->get("mult",i) << " "
	 << t->get("weight",i) << " " << t->get("x0",i) << " ";
    cout << t->get("x1",i) << " " << t->get("x0sq",i) << " "
	 << t->get("x0sq_x1sq",i) << endl;
  }
  cout.precision(6);
  cout << endl;

  n=t->get_nlines();

  t_avg=wvector_mean(n,t->get_column("x0sq"),
		     t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x0sq"),
			  t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,res[1],t_avgerr*10.0,"tab 1");

  t_avg=wvector_mean(n,t->get_column("x0sq_x1sq"),
		     t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x0sq_x1sq"),
			  t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,res[2],t_avgerr*10.0,"tab 2");
  cout << endl;
  
  tm.report();
  
  return 0;
}
// End of example
