/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2023, Andrew W. Steiner
  
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
#include <o2scl/mcmc_para.h>
#include <o2scl/vec_stats.h>
#include <o2scl/test_mgr.h>
#include <o2scl/mcarlo_miser.h>
#include <o2scl/multi_funct.h>
#include <o2scl/expval.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

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

typedef std::function<int(const ubvector &,double,size_t,int,bool,
			  std::array<double,1> &)> measure_funct;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  std::array<double,1> &)> fill_funct;

class mcmc_para_class {

public:

  int count;
  
  mcmc_para_base<point_funct,measure_funct,std::array<double,1>,
                     ubvector> mc;

  mcmc_para_table<point_funct,fill_funct,std::array<double,1>,
                      ubvector> mct;

  expval_scalar sev_x;
  
  expval_scalar sev_x2;

  double last_par;
  
  double last_dat;
  
  int gauss(size_t nv, const ubvector &pars, double &ret,
	    std::array<double,1> &dat) {
    dat[0]=pars[0]*pars[0];
    ret=-pars[0]*pars[0]/2.0;
    return o2scl::success;
  }

  int flat(size_t nv, const ubvector &pars, double &ret,
	   std::array<double,1> &dat) {
    dat[0]=pars[0]*pars[0];
    ret=0.0;
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

  int measure(const ubvector &pars, double log_weight, size_t ix,
	      int ret, bool new_meas, std::array<double,1> &dat) {
    if (new_meas) {
      sev_x.add(pars[0]);
      sev_x2.add(dat[0]);
      last_par=pars[0];
      last_dat=dat[0];
    } else {
      sev_x.add(last_par);
      sev_x2.add(last_dat);
    }
    return 0;
  }

  int fill_func(const ubvector &pars, double log_weight,
		std::vector<double> &line, std::array<double,1> &dat) {
    line.push_back(dat[0]);
    return 0;
  }
  
};

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(2);

  mcmc_para_class mpc;
  
  // Domain limits
  ubvector low(1);
  ubvector high(1);
  low[0]=-5.0;
  high[0]=2.0;

  // Compute exact results
  mcarlo_miser<> mm;
  mm.n_points=100000;
  double res[3], err[3];
  
  multi_funct mf0=std::bind
    (std::mem_fn<double(size_t,const ubvector &)>(&mcmc_para_class::f0),
     &mpc,std::placeholders::_1,std::placeholders::_2);
  multi_funct mf1=std::bind
    (std::mem_fn<double(size_t,const ubvector &)>(&mcmc_para_class::f1),
     &mpc,std::placeholders::_1,std::placeholders::_2);
  multi_funct mf2=std::bind
    (std::mem_fn<double(size_t,const ubvector &)>(&mcmc_para_class::f2),
     &mpc,std::placeholders::_1,std::placeholders::_2);
  
  mm.minteg_err(mf0,1,low,high,res[0],err[0]);
  mm.minteg_err(mf1,1,low,high,res[1],err[1]);
  mm.minteg_err(mf2,1,low,high,res[2],err[2]);
  err[1]=res[1]/res[0]*sqrt(pow(err[0]/res[0],2.0)+pow(err[1]/res[1],2.0));
  res[1]/=res[0];
  err[2]=res[2]/res[0]*sqrt(pow(err[0]/res[0],2.0)+pow(err[2]/res[2],2.0));
  res[2]/=res[0];
  cout << "Direct Monte Carlo results:" << endl;
  cout << res[1] << " " << err[1] << endl;
  cout << res[2] << " " << err[2] << endl;
  cout << endl;
  
  // Set up MCMC
  point_funct gauss_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,double &,
		     std::array<double,1> &)>(&mcmc_para_class::gauss),
     &mpc,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  point_funct flat_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,double &,
		     std::array<double,1> &)>(&mcmc_para_class::flat),
     &mpc,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  fill_funct ff=std::bind
    (std::mem_fn<int(const ubvector &,double,std::vector<double> &,
		     std::array<double,1> &)>(&mcmc_para_class::fill_func),
     &mpc,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  measure_funct mf=std::bind
    (std::mem_fn<int(const ubvector &,double,size_t,int,bool,
		     std::array<double,1> &)>
     (&mcmc_para_class::measure),&mpc,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4,std::placeholders::_5,std::placeholders::_6);

  size_t n_threads=1;
#ifdef O2SCL_OPENMP
  n_threads=2;
#endif

  vector<point_funct> gauss_vec(n_threads);
  vector<point_funct> flat_vec(n_threads);
  vector<measure_funct> meas_vec(n_threads);
  vector<fill_funct> fill_vec(n_threads);

  // Enough for 10 walkers, 2 threads, times 2
  vector<array<double,1> > data_vec(40);
  for(size_t i=0;i<n_threads;i++) {
    gauss_vec[i]=gauss_func;
    flat_vec[i]=flat_func;
    meas_vec[i]=mf;
    fill_vec[i]=ff;
  }
  
  cout << "n_threads: " << n_threads << endl;
  cout << endl;

  size_t N=20000;
  
  double avg, std, avg_err;
  size_t i1, i2;

  // ----------------------------------------------------------------
  // Plain MCMC

  if (true) {
    
    cout << "Plain MCMC: " << endl;

    mpc.mc.step_fac=10.0;
    mpc.mc.verbose=2;
    mpc.mc.n_threads=1;
    mpc.mc.max_iters=N;
    mpc.mc.prefix="mcmc";
    
    mpc.sev_x.set_blocks(20,N/20);
    mpc.sev_x2.set_blocks(20,N/20);

    mpc.mc.meas_for_initial=false;
    
    mpc.mc.mcmc(1,low,high,gauss_vec,meas_vec,data_vec);
    
    tm.test_gen(mpc.sev_x.finished(),"plain sev finished");
    
    mpc.sev_x.current_avg(avg,std,avg_err);
    cout << avg << " " << avg_err << endl;
    tm.test_rel(avg,res[1],100.0*sqrt(avg_err*avg_err+err[1]*err[1]),
		"plain mcmc 1");
    mpc.sev_x2.current_avg(avg,std,avg_err);
    cout << avg << " " << avg_err << endl;
    tm.test_rel(avg,res[2],4.0*sqrt(avg_err*avg_err+err[2]*err[2]),
		"plain mcmc 2");
    
    tm.test_gen(mpc.mc.n_accept[0]+mpc.mc.n_reject[0]==mpc.mc.max_iters,
		"plain n_iters 0");
    if (mpc.mc.n_threads>1) {
      tm.test_gen(mpc.mc.n_accept[1]+mpc.mc.n_reject[1]==mpc.mc.max_iters,
		  "plain n_iters 1");
    }
    cout << endl;

    mpc.sev_x.free();
    mpc.sev_x2.free();
    
  }

  // ----------------------------------------------------------------
  // Affine-invariant MCMC

  if (true) {
    cout << "Affine-invariant MCMC: " << endl;

    mpc.sev_x.set_blocks(40,1);
    mpc.sev_x2.set_blocks(40,1);
    
    mpc.sev_x.get_block_indices(i1,i2);
    cout << i1 << " " << i2 << endl;

    mpc.mc.aff_inv=true;
    mpc.mc.n_walk=10;
    mpc.mc.step_fac=2.0;
    mpc.mc.verbose=2;
    mpc.mc.n_threads=1;
    mpc.mc.max_iters=N*10;
    mpc.mc.prefix="mcmc_ai";
    
    mpc.mc.mcmc(1,low,high,gauss_vec,meas_vec,data_vec);

    mpc.sev_x.get_block_indices(i1,i2);
    cout << i1 << " " << i2 << endl;
    
    //tm.test_gen(mpc.sev_x.finished(),"aff_inv sev finished");
    
    mpc.sev_x.current_avg_stats(avg,std,avg_err,i1,i2);
    cout << avg << " " << avg_err << " " << i1 << " " << i2 << endl;
    //tm.test_rel(avg,res[1],8.0*sqrt(avg_err*avg_err+err[1]*err[1]),
    //"aff_inv mcmc 1");
    mpc.sev_x2.current_avg_stats(avg,std,avg_err,i1,i2);
    cout << avg << " " << avg_err << " " << i1 << " " << i2 << endl;
    //tm.test_rel(avg,res[2],4.0*sqrt(avg_err*avg_err+err[2]*err[2]),
    //"aff_inv mcmc 2");
    
    tm.test_gen(mpc.mc.n_accept[0]+mpc.mc.n_reject[0]==mpc.mc.max_iters,
		"aff_inv n_iters 0");
    if (mpc.mc.n_threads>1) {
      tm.test_gen(mpc.mc.n_accept[1]+mpc.mc.n_reject[1]==mpc.mc.max_iters,
		  "aff_inc n_iters 1");
    }
    cout << endl;

    mpc.sev_x.free();
    mpc.sev_x2.free();
    
  }

  // ----------------------------------------------------------------
  // Plain MCMC with a table

  cout << "Plain MCMC with a table: " << endl;
  
  vector<string> pnames={"x","x2"};
  vector<string> punits={"MeV","MeV^2"};
  mpc.mct.set_names_units(pnames,punits);

  mpc.mct.aff_inv=false;
  mpc.mct.step_fac=10.0;
  mpc.mct.verbose=2;
  mpc.mct.n_threads=n_threads;
  mpc.mct.max_iters=N;
  mpc.mct.prefix="mcmct";
  mpc.mct.table_prealloc=N*n_threads;

  mpc.mct.mcmc_fill(1,low,high,gauss_vec,fill_vec,data_vec);

  std::shared_ptr<o2scl::table_units<> > table=mpc.mct.get_table();
  
  mpc.sev_x.free();
  mpc.sev_x2.free();
  mpc.sev_x.set_blocks(40,1);
  mpc.sev_x2.set_blocks(40,1);
  for(size_t i=0;i<table->get_nlines();i++) {
    for(size_t j=0;j<((size_t)(table->get("mult",i)+1.0e-8));j++) {
      mpc.sev_x.add(table->get("x",i));
      mpc.sev_x2.add(table->get("x2",i));
    }
  }
  
  mpc.sev_x.current_avg_stats(avg,std,avg_err,i1,i2);
  cout << avg << " " << avg_err << " " << i1 << " " << i2 << endl;
  tm.test_rel(avg,res[1],100.0*sqrt(avg_err*avg_err+err[1]*err[1]),
	      "plain table mcmc 1");
  mpc.sev_x2.current_avg_stats(avg,std,avg_err,i1,i2);
  cout << avg << " " << avg_err << " " << i1 << " " << i2 << endl;
  tm.test_rel(avg,res[2],10.0*sqrt(avg_err*avg_err+err[2]*err[2]),
	      "plain table mcmc 2");
  
  std::vector<size_t> chain_sizes;
  mpc.mct.get_chain_sizes(chain_sizes);
  tm.test_gen(chain_sizes[0]==mpc.mct.n_accept[0]+1,"plain table size 0");
  if (mpc.mct.n_threads>1) {
    tm.test_gen(chain_sizes[1]==mpc.mct.n_accept[1]+1,"plain table size 1");
  }
  
  tm.test_gen(mpc.mct.n_accept[0]+mpc.mct.n_reject[0]==mpc.mct.max_iters,
              "plain table n_iters 0");
  if (mpc.mct.n_threads>1) {
    tm.test_gen(mpc.mct.n_accept[1]+mpc.mct.n_reject[1]==mpc.mct.max_iters,
		"plain table n_iters 1");
  }
  cout << endl;
    
  // ----------------------------------------------------------------
  // Affine-invariant MCMC with a table
  
  cout << "Affine-invariant MCMC with a table: " << endl;
  
  mpc.mct.aff_inv=true;
  mpc.mct.n_walk=10;
  mpc.mct.step_fac=2.0;
  mpc.mct.verbose=2;
  mpc.mct.n_threads=n_threads;
  mpc.mct.max_iters=N;
  mpc.mct.prefix="mcmct_ai";
  mpc.mct.table_prealloc=N*n_threads;

  mpc.mct.mcmc_fill(1,low,high,gauss_vec,fill_vec,data_vec);

  // Get results
  table=mpc.mct.get_table();

  // Test results
  mpc.sev_x.free();
  mpc.sev_x2.free();
  mpc.sev_x.set_blocks(40,1);
  mpc.sev_x2.set_blocks(40,1);
  for(size_t i=0;i<table->get_nlines();i++) {
    for(size_t j=0;j<((size_t)(table->get("mult",i)+1.0e-8));j++) {
      mpc.sev_x.add(table->get("x",i));
      mpc.sev_x2.add(table->get("x2",i));
    }
  }
  
  mpc.sev_x.current_avg_stats(avg,std,avg_err,i1,i2);
  cout << avg << " " << avg_err << " " << i1 << " " << i2 << endl;
  tm.test_rel(avg,res[1],100.0*sqrt(avg_err*avg_err+err[1]*err[1]),
	      "aff_inv table mcmc 1");
  mpc.sev_x2.current_avg_stats(avg,std,avg_err,i1,i2);
  cout << avg << " " << avg_err << " " << i1 << " " << i2 << endl;
  tm.test_rel(avg,res[2],4.0*sqrt(avg_err*avg_err+err[2]*err[2]),
	      "aff_inv table mcmc 2");
  
  mpc.mct.get_chain_sizes(chain_sizes);

  tm.test_gen(vector_sum_double(mpc.mct.n_walk,chain_sizes)-mpc.mct.n_walk==
	      mpc.mct.n_accept[0],"accept chain 0");
  if (n_threads>1) {
    tm.test_gen(vector_sum_double(mpc.mct.n_walk*2,chain_sizes)-
		vector_sum_double(mpc.mct.n_walk,chain_sizes)-mpc.mct.n_walk==
		mpc.mct.n_accept[1],"accept chain 1");
  }
  tm.test_gen(mpc.mct.n_accept[0]+mpc.mct.n_reject[0]==mpc.mct.max_iters,
	      "aff_inv table n_iters 0");
  if (n_threads>1) {
    tm.test_gen(mpc.mct.n_accept[1]+mpc.mct.n_reject[1]==mpc.mct.max_iters,
		"aff_inv table n_iters 1");
  }
  cout << endl;

  if (true) {
    
    // ----------------------------------------------------------------
    // Plain MCMC with a table and a flat distribution
  
    cout << "Plain MCMC with a table and a flat distribution: "
	 << endl;

    mpc.count=0;
    mpc.mct.aff_inv=false;
    mpc.mct.n_walk=1;
    mpc.mct.step_fac=10.0;
    mpc.mct.verbose=2;
    mpc.mct.max_iters=N;
    mpc.mct.n_threads=n_threads;
    mpc.mct.prefix="mcmct_flat";
    mpc.mct.table_prealloc=N*n_threads;
    
    mpc.mct.mcmc_fill(1,low,high,flat_vec,fill_vec,data_vec);

    uniform_grid_end<double> ug(low[0],high[0],20);
    size_t hist_ix=0;
    hist h;
    h.set_bin_edges(ug);
    expval_vector ev(20,1,20);
    size_t count=0;
    std::shared_ptr<o2scl::table_units<> > tablex=mpc.mct.get_table();
    bool done=false;
    for (size_t j=0;j<tablex->get_nlines() && done==false;j++) {
      for(size_t k=0;k<((size_t)(tablex->get("mult",j)+1.0e-8));k++) {
	if (count/(mpc.mct.max_iters*n_threads/20)>hist_ix) {
	  ev.add(h.get_wgts());
	  h.clear_wgts();
	  hist_ix++;
	}
	if (hist_ix==20) done=true;
	h.update(tablex->get("x",j));
	count++;
      }
    }
    if (done==false) {
      ev.add(h.get_wgts());
    }

    ubvector vavg(20), vstd(20), vavge(20);
    ev.current_avg_stats(vavg,vstd,vavge,i1,i2);
    cout << vavg[0] << endl;
    for(size_t j=1;j<20;j++) {
      cout << vavg[j] << endl;
      tm.test_rel(vavg[j],vavg[0],sqrt(vavg[0]+vavg[j]),"flat dist");
    }
    cout << endl;

  }

  if (true) {
    
    // ----------------------------------------------------------------
    // Affine-invariant MCMC with a table and a flat distribution
  
    cout << "Affine-invariant MCMC with a table and a flat distribution: "
	 << endl;
  
    mpc.mct.aff_inv=true;
    mpc.mct.n_walk=10;
    mpc.mct.step_fac=2.0;
    mpc.mct.verbose=2;
    mpc.mct.max_iters=N;
    mpc.mct.n_threads=n_threads;
    mpc.mct.prefix="mcmct_ai_flat";
    mpc.mct.table_prealloc=N*n_threads;
    
    mpc.mct.mcmc_fill(1,low,high,flat_vec,fill_vec,data_vec);

    uniform_grid_end<double> ug(low[0],high[0],20);
    size_t hist_ix=0;
    hist h;
    h.set_bin_edges(ug);
    expval_vector ev(20,1,20);
    size_t count=0;
    std::shared_ptr<o2scl::table_units<> > tabley=mpc.mct.get_table();
    bool done=false;
    for (size_t j=0;j<tabley->get_nlines() && done==false;j++) {
      for(size_t k=0;k<((size_t)(tabley->get("mult",j)+1.0e-8));k++) {
	if (count/(mpc.mct.max_iters*n_threads/20)>hist_ix) {
	  ev.add(h.get_wgts());
	  h.clear_wgts();
	  hist_ix++;
	}
	if (hist_ix==20) done=true;
	h.update(tabley->get("x",j));
	count++;
      }
    }
    if (done==false) {
      ev.add(h.get_wgts());
    }

    ubvector vavg(20), vstd(20), vavge(20);
    ev.current_avg_stats(vavg,vstd,vavge,i1,i2);
    cout << vavg[0] << endl;
    for(size_t j=1;j<20;j++) {
      cout << vavg[j] << endl;
      tm.test_rel(vavg[j],vavg[0],sqrt(vavg[0]+vavg[j]),"flat dist");
    }

  }

  if (false) {
    
    // ----------------------------------------------------------------
    // Affine-invariant MCMC with a table and previously read results
  
    cout << "Affine-invariant MCMC with a table and previous results: "
	 << endl;
  
    mpc.mct.aff_inv=true;
    mpc.mct.n_walk=10;
    mpc.mct.step_fac=-1.0;
  
    mpc.mct.max_iters=40;
    mpc.mct.prefix="mcmct_aiprev";

    hdf_file hf;
    string fname="mcmct_ai_0_out";
    hf.open(fname);
    mpc.mct.read_prev_results(hf,1);
    hf.close();

    cout << "Going to mcmc." << endl;
    mpc.mct.verbose=1;
    mpc.mct.mcmc_fill(1,low,high,gauss_vec,fill_vec,data_vec);

    // Get results
    table=mpc.mct.get_table();
    std::vector<size_t> chain_sizes2;
    mpc.mct.get_chain_sizes(chain_sizes2);

    // This testing code doesn't work yet, possibly because it is
    // not yet written correctly
    for(size_t it=0;it<n_threads;it++) {
      size_t sum1=0, sum2=0;
      for(size_t i=0;i<10;i++) {
	sum1+=chain_sizes[it*10+i];
      }
      for(size_t i=0;i<10;i++) {
	sum2+=chain_sizes2[it*10+i];
      }
      cout << sum1 << " " << sum2 << " " << sum2-sum1 << " " 
	   << mpc.mct.n_accept[it] << " "
	   << mpc.mct.n_reject[it] << endl;
      //tm.test_gen(sum2-sum1==mpc.mct.n_accept[it],"Test chain size");
    }
    cout << endl;
  }
  
  tm.report();
  
  return 0;
}

#endif
