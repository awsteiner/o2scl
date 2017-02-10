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
#include <o2scl/mcmc_omp.h>
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

int niters=1000000;
int it_count;
expval_scalar sev_x, sev_x2;
mcmc_omp_base<point_funct,measure_funct,std::array<double,1>,ubvector> mc;

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
    cout << "Returning mcmc_done." << endl;
    return mcmc_omp_base<point_funct,measure_funct,int,ubvector>::mcmc_done;
  }
  it_count++;
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

  size_t n_threads=2;

  vector<point_funct> vpf(n_threads);
  vector<measure_funct> vmf(n_threads);
  for(size_t i=0;i<n_threads;i++) {
    vpf[i]=pf;
    vmf[i]=mf;
  }
  
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
  mc.verbose=2;
  mc.n_threads=2;
  niters=20;
  cout << "n_threads: " << n_threads << endl;
  mc.mcmc(1,init,low,high,vpf,vmf);
  exit(-1);

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
  cout << mc.n_accept[0] << " " << mc.n_reject[0] << endl;
  tm.test_abs(avg,res[2],avg_err*10.0,"plain 2");
  cout << endl;

  // Clear expectation value objects for next iteration
  sev_x.free();
  sev_x2.free();
  exit(-1);

  // ----------------------------------------------------------------
  // Affine-invariant MCMC

  cout << "Affine-invariant MCMC: " << endl;
  
  mc.aff_inv=true;
  mc.n_walk=10;
  mc.step_fac=20.0;
  
  it_count=0;
  mc.mcmc(1,init,low,high,vpf,vmf);
  
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
  cout << mc.n_accept[0] << " " << mc.n_reject[0] << endl;
  tm.test_abs(avg,res[2],avg_err*10.0,"ai 2");
  cout << endl;
  
  sev_x.free();
  sev_x2.free();

  tm.report();
  
  return 0;
}
#endif
