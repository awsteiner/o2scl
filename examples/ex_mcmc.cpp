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
#include <o2scl/cubature.h>

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

/** 
 */
int f_cub(unsigned ndim, size_t npt, const double *x, unsigned fdim,
	  double *fval) {
  for (size_t i=0;i<npt;i++) {
    const double *x2=x+i*ndim;
    double *f2=fval+i*fdim;
    f2[0]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)));
    f2[1]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)))*x2[0]*x2[0];
    f2[2]=exp(-((x2[0]-0.2)*(x2[0]-0.2)+
		(x2[1]-0.5)*(x2[1]-0.5)))*x2[0]*x2[0]*x2[1]*x2[1];
  }
  return 0;
}

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(1);

  // Parameter limits and initial point
  ubvector low(2), high(2), init(2);
  low[0]=-2.0;
  low[1]=-2.0;
  high[0]=2.0;
  high[1]=2.0;
  init[0]=0.2;
  init[1]=0.5;

  // Use cubature to compute integrals
  double dlow[2]={-2.0,-2.0};
  double dhigh[2]={2.0,2.0};
  double dres[3], derr[3];
  typedef std::function<
    int(unsigned,size_t,const double *,unsigned,double *)> cub_funct_arr;
  cub_funct_arr cfa=f_cub;
  inte_hcubature<cub_funct_arr> hc;
  inte_hcubature<cub_funct_arr>::error_norm enh=
    inte_hcubature<cub_funct_arr>::ERROR_INDIVIDUAL;
  int ret=hc.integ(3,cfa,2,dlow,dhigh,10000,0.0,1.0e-4,enh,dres,derr);
  tm.test_gen(ret==0,"cubature success.");
  
  // Exact results
  double exact_res[2];
  exact_res[0]=dres[1]/dres[0];
  exact_res[1]=dres[2]/dres[0];

  point_funct pf=point;
  fill_funct ff=fill_line;

  // Table column names and units. We must specify first the names for
  // the parameters first and then the names of the auxillary
  // parameters.
  vector<string> pnames={"x0","x1","x0sq","x0sq_x1sq"};
  vector<string> punits={"MeV","MeV","MeV^2","MeV^4"};
  mct.set_names_units(pnames,punits);

  // Get a pointer to the results table
  shared_ptr<table_units<> > t=mct.get_table();

  // MCMC with a random walk of a fixed length
  mct.step_fac=3.0;
  mct.mcmc(2,init,low,high,pf,ff);

  cout << "n_accept, n_reject, table lines: "
       << mct.n_accept << " " << mct.n_reject << " "
       << t->get_nlines() << endl;
  cout.precision(4);
  for(size_t i=0;i<t->get_nlines();i+=t->get_nlines()/10) {
    cout.width(5);
    cout << i << " " << t->get("mult",i) << " ";
    cout.setf(ios::showpos);
    cout << t->get("log_wgt",i) << " " << t->get("x0",i) << " ";
    cout << t->get("x1",i) << " " << t->get("x0sq",i) << " "
	 << t->get("x0sq_x1sq",i) << endl;
    cout.unsetf(ios::showpos);
  }
  cout.precision(6);
  cout << endl;

  /*
    This function averages the table into 40 blocks to remove
    autocorrelations. The autocorrelation length is not computed here
    (but is small enough for this example).
  */
  mct.reblock(40);

  // Compute and test the average values
  size_t n=t->get_nlines();
  double t_avg=wvector_mean(n,t->get_column("x0sq"),
			    t->get_column("mult"));
  double t_stddev=wvector_stddev(n,t->get_column("x0sq"),
				 t->get_column("mult"));
  double t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);

  tm.test_abs(t_avg,exact_res[0],t_avgerr*10.0,"tab 1");

  t_avg=wvector_mean(n,t->get_column("x0sq_x1sq"),
		     t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x0sq_x1sq"),
			  t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,exact_res[1],t_avgerr*10.0,"tab 2");
  cout << endl;

  // MCMC with affine-invariant sampling
  mct.aff_inv=true;
  mct.n_walk=10;
  mct.step_fac=5.0;

  mct.mcmc(2,init,low,high,pf,ff);

  cout << "n_accept, n_reject, table lines: "
       << mct.n_accept << " " << mct.n_reject << " "
       << t->get_nlines() << endl;
  cout.precision(4);
  for(size_t i=0;i<t->get_nlines();i+=t->get_nlines()/10) {
    cout.width(5);
    cout << i << " " << t->get("mult",i) << " ";
    cout.setf(ios::showpos);
    cout << t->get("log_wgt",i) << " " << t->get("x0",i) << " ";
    cout << t->get("x1",i) << " " << t->get("x0sq",i) << " "
	 << t->get("x0sq_x1sq",i) << endl;
    cout.unsetf(ios::showpos);
  }
  cout.precision(6);
  cout << endl;

  mct.reblock(40);

  n=t->get_nlines();

  t_avg=wvector_mean(n,t->get_column("x0sq"),
		     t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x0sq"),
			  t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,exact_res[0],t_avgerr*10.0,"tab 1");

  t_avg=wvector_mean(n,t->get_column("x0sq_x1sq"),
		     t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x0sq_x1sq"),
			  t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,exact_res[1],t_avgerr*10.0,"tab 2");
  cout << endl;
  
  tm.report();
  
  return 0;
}
// End of example
