/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2017, Andrew W. Steiner
  
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
#include <o2scl/mcarlo_miser.h>
#include <o2scl/multi_funct.h>
#include <o2scl/expval.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;

#ifdef O2SCL_OLDER_COMPILER
int main(void) {
  test_mgr t;
  t.report();
  return 0;
}
#else

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

typedef std::function<int(size_t,const ubvector &,double &,
			  std::array<double,1> &)> point_funct;

typedef std::function<int(const ubvector &,double,size_t,bool,
			  std::array<double,1> &)> measure_funct;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  std::array<double,1> &)> fill_funct;

static const int niters=1000000;
int it_count;
expval_scalar sev_x, sev_x2;
mcmc_base<point_funct,measure_funct,std::array<double,1>,ubvector> mc;
mcmc_table<point_funct,fill_funct,std::array<double,1>,ubvector> mct;

int point(size_t nv, const ubvector &pars, double &ret,
	  std::array<double,1> &dat) {
  dat[0]=pars[0]*pars[0];
  ret=-pars[0]*pars[0]/2.0;
  return o2scl::success;
}

double f0(size_t nv, const ubvector &pars) {
  return exp(-pars[0]*pars[0]/2.0);
}

double f1(size_t nv, const ubvector &pars) {
  return exp(-pars[0]*pars[0]/2.0)*pars[0];
}

double f2(size_t nv, const ubvector &pars) {
  return exp(-pars[0]*pars[0]/2.0)*pars[0]*pars[0];
}

int measure(const ubvector &pars, double log_weight, size_t ix, bool new_meas,
	    std::array<double,1> &dat) {
  sev_x.add(pars[0]);
  sev_x2.add(dat[0]);
  // Check that the 'dat' object is correctly filled
  if ((pars[0]*pars[0]-dat[0])>1.0e-10) {
    cerr << "Failure." << endl;
    exit(-1);
  }
  if (it_count==niters-1) {
    return mcmc_base<point_funct,measure_funct,int,ubvector>::mcmc_done;
  }
  it_count++;
  return 0;
}

int fill_func(const ubvector &pars, double log_weight,
	      std::vector<double> &line, std::array<double,1> &dat) {
  line.push_back(dat[0]);
  // Check that the 'dat' object is correctly filled
  if ((pars[0]*pars[0]-dat[0])>1.0e-10) {
    cerr << "Failure 2." << endl;
    exit(-1);
  }
  if (mct.get_table()->get_nlines()==niters/3) {
    return mcmc_base<point_funct,measure_funct,int,ubvector>::mcmc_done;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(1);

  // Domain limits
  ubvector low(1);
  ubvector high(1);
  low[0]=-5.0;
  high[0]=2.0;

  // Compute exact results
  mcarlo_miser<> mm;
  mm.n_points=100000;
  multi_funct11 mf0=f0;
  multi_funct11 mf1=f1;
  multi_funct11 mf2=f2;

  double res[3], err[3];
  mm.minteg_err(mf0,1,low,high,res[0],err[0]);
  mm.minteg_err(mf1,1,low,high,res[1],err[1]);
  mm.minteg_err(mf2,1,low,high,res[2],err[2]);
  res[1]/=res[0];
  res[2]/=res[0];
  cout << "Exact results:" << endl;
  cout << res[1] << endl;
  cout << res[2] << endl;
  cout << endl;
  
  // Set up MCMC
  point_funct pf=point;
  measure_funct mf=measure;
  fill_funct ff=fill_func;
    
  ubvector init(1);
  init[0]=-0.01;

  // Set up expectation value objects
  sev_x.set_blocks(40,niters/40);
  sev_x2.set_blocks(40,niters/40);
  double avg, std_dev, avg_err;
  size_t m_block, m_per_block;

  // ----------------------------------------------------------------
  // Plain MCMC

  cout << "Plain MCMC: " << endl;
  
  it_count=0;
  mc.step_fac=2.0;
  mc.mcmc(1,init,low,high,pf,mf);

  sev_x.current_avg_stats(avg,std_dev,avg_err,m_block,m_per_block);

  // Output and test results
  cout.setf(ios::showpos);
  cout << avg << " " << std_dev << " " << avg_err << " ";
  cout.unsetf(ios::showpos);
  cout << m_block << " " << m_per_block << endl;
  tm.test_abs(avg,res[1],avg_err*10.0,"plain 1");
  sev_x2.current_avg_stats(avg,std_dev,avg_err,m_block,m_per_block);
  cout.setf(ios::showpos);
  cout << avg << " " << std_dev << " " << avg_err << " ";
  cout.unsetf(ios::showpos);
  cout << m_block << " " << m_per_block << endl;
  cout << mc.n_accept << " " << mc.n_reject << endl;
  tm.test_abs(avg,res[2],avg_err*10.0,"plain 2");
  cout << endl;

  // Clear expectation value objects for next iteration
  sev_x.free();
  sev_x2.free();

  // ----------------------------------------------------------------
  // Affine-invariant MCMC

  cout << "Affine-invariant MCMC: " << endl;
  
  mc.aff_inv=true;
  mc.n_walk=10;
  mc.step_fac=20.0;
  
  it_count=0;
  mc.mcmc(1,init,low,high,pf,mf);
  
  sev_x.current_avg_stats(avg,std_dev,avg_err,m_block,m_per_block);

  // Output and test results
  cout.setf(ios::showpos);
  cout << avg << " " << std_dev << " " << avg_err << " ";
  cout.unsetf(ios::showpos);
  cout << m_block << " " << m_per_block << endl;
  tm.test_abs(avg,res[1],avg_err*10.0,"ai 1");
  sev_x2.current_avg_stats(avg,std_dev,avg_err,m_block,m_per_block);
  cout.setf(ios::showpos);
  cout << avg << " " << std_dev << " " << avg_err << " ";
  cout.unsetf(ios::showpos);
  cout << m_block << " " << m_per_block << endl;
  cout << mc.n_accept << " " << mc.n_reject << endl;
  tm.test_abs(avg,res[2],avg_err*10.0,"ai 2");
  cout << endl;
  
  sev_x.free();
  sev_x2.free();

  // ----------------------------------------------------------------
  // Plain MCMC with table

  cout << "Plain MCMC with table:" << endl;

  vector<string> pnames={"x","x2"};
  vector<string> punits={"MeV","MeV^2"};
  mct.set_names_units(pnames,punits);

  it_count=0;
  mct.step_fac=2.0;
  mct.mcmc(1,init,low,high,pf,ff);

  shared_ptr<table_units<> > t=mct.get_table();
  
  mct.reblock(40);
  
  size_t n=t->get_nlines();

  double t_avg=wvector_mean(n,t->get_column("x"),t->get_column("mult"));
  double t_stddev=wvector_stddev(n,t->get_column("x"),t->get_column("mult"));
  double t_avgerr=t_stddev/sqrt((double)n);
  
  // Output and test results
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,res[1],t_avgerr*10.0,"tab 1");
  
  t_avg=wvector_mean(n,t->get_column("x2"),t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x2"),t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,res[2],t_avgerr*10.0,"tab 1");
  
  cout << mct.n_accept << " " << mct.n_reject << endl;
  cout << endl;
  
  // ----------------------------------------------------------------
  // Affine-invariant MCMC with table

  cout << "Affine-invariant MCMC with table: " << endl;

  mct.aff_inv=true;
  mct.n_walk=10;
  mct.step_fac=20.0;

  it_count=0;
  mct.mcmc(1,init,low,high,pf,ff);

  mct.reblock(40);
  
  n=t->get_nlines();

  t_avg=wvector_mean(n,t->get_column("x"),t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x"),t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  // Output and test results
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,res[1],t_avgerr*10.0,"tab 1");
  
  t_avg=wvector_mean(n,t->get_column("x2"),t->get_column("mult"));
  t_stddev=wvector_stddev(n,t->get_column("x2"),t->get_column("mult"));
  t_avgerr=t_stddev/sqrt((double)n);
  
  cout.setf(ios::showpos);
  cout << t_avg << " " << t_stddev << " " << t_avgerr << endl;
  cout.unsetf(ios::showpos);
  tm.test_abs(t_avg,res[2],t_avgerr*10.0,"tab 1");
  
  cout << mct.n_accept << " " << mct.n_reject << endl;

  // ----------------------------------------------------------------

  tm.report();
  
  return 0;
}
#endif
