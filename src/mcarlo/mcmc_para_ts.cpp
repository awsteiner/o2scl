/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2016-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/mcmc_para.h>
#include <o2scl/vec_stats.h>
#include <o2scl/test_mgr.h>
#include <o2scl/mcarlo_miser.h>
#include <o2scl/multi_funct.h>
#include <o2scl/expval.h>
#include <o2scl/hdf_io.h>
#include <o2scl/interpm_krige.h>
#include <o2scl/kde_python.h>

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
			  std::vector<double> &)> point_funct;

typedef std::function<int(const ubvector &,double,size_t,int,bool,
			  std::vector<double> &)> measure_funct;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  std::vector<double> &)> fill_funct;

typedef std::function<int(size_t,const ubvector &,double &,
			  std::vector<double> &)> point_hmc;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  std::vector<double> &)> fill_hmc;

class mcmc_para_class {

public:

  int count;
  
  mcmc_para_base<point_funct,measure_funct,std::vector<double>,
                 ubvector> mc;
  
  mcmc_para_table<point_funct,fill_funct,std::vector<double>,
                  ubvector> mct;
  
  mcmc_para_table<point_hmc,fill_hmc,std::vector<double>,ubvector>
  mct_hmc;
  
  int hmc_point(size_t n_params, const ubvector &u, double &log_wgt,
                std::vector<double> &dat) {
    
    static const double mean_x=0.0;
    static const double mean_y=0.0;
    static const double width_x=1.0;
    static const double width_y=1.0;
    static const double rho=0.0;
    
    double shift_x=u[0]-mean_x;
    double shift_y=u[1]-mean_y;
    double x=pow(shift_x/width_x,2.0);
    double y=pow(shift_y/width_y,2.0);
    double z=2.0*rho*shift_x*shift_y/(width_x*width_y);
    double r=1.0-rho*rho;
    double c=2.0*o2scl_const::pi*width_x*width_y*sqrt(r);

    /*
    cout.precision(10);
    cout << "bp: " << u[0] << " " << u[1] << " "
         << z << " " << r << " " << c << " "
         << 1.0/c*exp(-0.5*(x+y-z)/r) << std::endl;
    cout.precision(6);
    */
    
    log_wgt=log(1.0/c*exp(-0.5*(x+y-z)/r));

    dat[0]=x;
    dat[1]=y;

    if (false) {
      std::cout << "hmc_point: " << u[0] << " " << u[1] << " "
                << log_wgt << " " << dat[0] << " " << dat[1] << std::endl;
    }
    
    return 0;
  }
  
  int hmc_fill(const ubvector &pars, double log_weight,
               std::vector<double> &line, std::vector<double> &dat) {
    for(size_t j=0;j<dat.size();j++) {
      line.push_back(dat[j]);
    }
    return 0;
  }
  
  mcmc_para_class() {
  }
  
  expval_scalar sev_x;
  
  expval_scalar sev_x2;

  double last_par;
  
  double last_dat;
  
  int gauss(size_t nv, const ubvector &pars, double &ret,
	    std::vector<double> &dat) {
    dat[0]=pars[0]*pars[0];
    ret=-pars[0]*pars[0]/2.0;
    return o2scl::success;
  }

  int flat(size_t nv, const ubvector &pars, double &ret,
	   std::vector<double> &dat) {
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
	      int ret, bool new_meas, std::vector<double> &dat) {
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
		std::vector<double> &line, std::vector<double> &dat) {
    line.push_back(dat[0]);
    return 0;
  }
  
};

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);

  test_mgr tm;
  tm.set_output_level(1);

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
		     std::vector<double> &)>(&mcmc_para_class::gauss),
     &mpc,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  point_funct flat_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,double &,
		     std::vector<double> &)>(&mcmc_para_class::flat),
     &mpc,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  fill_funct ff=std::bind
    (std::mem_fn<int(const ubvector &,double,std::vector<double> &,
		     std::vector<double> &)>(&mcmc_para_class::fill_func),
     &mpc,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  measure_funct mf=std::bind
    (std::mem_fn<int(const ubvector &,double,size_t,int,bool,
		     std::vector<double> &)>
     (&mcmc_para_class::measure),&mpc,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4,std::placeholders::_5,std::placeholders::_6);

  size_t n_threads=1;
#ifdef O2SCL_SET_OPENMP
  n_threads=2;
#endif

  vector<point_funct> gauss_vec(n_threads);
  vector<point_funct> flat_vec(n_threads);
  vector<measure_funct> meas_vec(n_threads);
  vector<fill_funct> fill_vec(n_threads);

  // Enough memory for 10 walkers and 2 threads, times 2
  vector<std::vector<double> > data_vec(40);
  for(size_t j=0;j<40;j++) {
    // We need enough space for 2 outputs because the hmc_point()
    // function has two outputs
    data_vec[j].resize(2);
  }
  
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

  mpc.mc.def_stepper->step_fac.resize(1);
  
  // ----------------------------------------------------------------
  // Plain MCMC

  if (true) {
    
    cout << "Plain MCMC: " << endl;

    mpc.mc.step_fac=10.0;
    mpc.mc.verbose=2;
    mpc.mc.n_threads=1;
    mpc.mc.max_iters=N;
    mpc.mc.prefix="mcmc";
    mpc.mc.def_stepper->step_fac[0]=10.0;
    
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

  if (true) {
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
    mpc.mct.def_stepper->step_fac[0]=10.0;

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

  }

  if (true) {
  
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
    std::shared_ptr<o2scl::table_units<> > table=mpc.mct.get_table();
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
  
    std::vector<size_t> chain_sizes;
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

  }

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
    mpc.mct.def_stepper->step_fac[0]=10.0;
    
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
    
    cout << endl;
  }

  if (true) {
    
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
    std::shared_ptr<o2scl::table_units<> > table=mpc.mct.get_table();
    std::vector<size_t> chain_sizes2;
    mpc.mct.get_chain_sizes(chain_sizes2);

    cout << endl;
  }

#ifdef O2SCL_SET_PYTHON
  
  if (true) {

    //n_threads=1;
    
    // ----------------------------------------------------------------
    // Independence MH with a KDE
    
    cout << "Independence MH with a KDE:" << endl;

    vector<string> pnames_imh={"x","x2"};
    vector<string> punits_imh={"MeV","MeV^2"};
    
    mpc.mct.stepper=shared_ptr<mcmc_stepper_base<point_hmc,
                               std::vector<double>,ubvector>>
      (new mcmc_stepper_mh<point_hmc,std::vector<double>,
       ubvector,ubmatrix,
       prob_cond_mdim_indep<>>);
    mpc.mct.set_names_units(pnames_imh,punits_imh);

    vector<vector<double>> data_vec_imh(2);
    data_vec_imh[0].resize(2);
    data_vec_imh[1].resize(2);

    mpc.mct.aff_inv=false;
    mpc.mct.verbose=3;
    mpc.mct.n_threads=1;
    mpc.mct.max_iters=200;
    mpc.mct.prefix="mcmct_imh_kde";

    // Read the previous table
    table_units<> last;
    hdf_file hfx;
    hfx.open("mcmct_0_out");
    hdf_input(hfx,last);
    hfx.close();
    
    // Train the KDE with the file, creating a new KDE object for each
    // thread, and then setting that KDE as the base distribution for
    // the independent conditional probability. 
    std::shared_ptr<mcmc_stepper_mh<point_hmc,std::vector<double>,
       ubvector,ubmatrix,prob_cond_mdim_indep<>>> new_stepper
      (new mcmc_stepper_mh<point_hmc,std::vector<double>,
       ubvector,ubmatrix,prob_cond_mdim_indep<>>);
    mpc.mct.stepper=new_stepper;
    new_stepper->proposal.resize(1);
    
    // Copy the table data to a tensor for use in kde_python.
    // We need a copy for each thread because kde_python takes
    // over the tensor data.
    tensor<> tin;
    vector<size_t> in_size={last.get_nlines(),1};
    tin.resize(2,in_size);
    
    for(size_t i=0;i<last.get_nlines();i++) {
      vector<size_t> ix;
      ix={i,0};
      tin.get(ix)=last.get("x",i);
    }
    
    vector<double> weights;
    std::shared_ptr<kde_python<ubvector>> kp(new kde_python<ubvector>);
    kp->set_function("o2sclpy",tin,
                     weights,"verbose=0","kde_scipy");
    
    new_stepper->proposal[0].set_base(kp);
    
    // Read initial points from the file
    mpc.mct.initial_points_file_last("mcmct_0_out",1);
    
    // Run MCMC
    mpc.mct.mcmc_fill(1,low,high,gauss_vec,fill_vec,data_vec);

    cout << endl;
    
  }
  
#endif
  
  if (true) {
    
    // ----------------------------------------------------------------
    // HMC with a table
    
    cout << "HMC with a table: " << endl;
  
    vector<string> pnames_hmc={"x","y","d1","d2"};
    vector<string> punits_hmc={"km","fm","cm","mm"};
    vector<vector<double>> data_vec_hmc(2);
    data_vec_hmc[0].resize(2);
    data_vec_hmc[1].resize(2);
    ubvector low_hmc(2), high_hmc(2);
    low_hmc[0]=-10.0;
    low_hmc[1]=-10.0;
    high_hmc[0]=10.0;
    high_hmc[1]=10.0;

    std::shared_ptr<
      mcmc_stepper_hmc<
      point_hmc,std::vector<double>,ubvector>> new_stepper
      (new mcmc_stepper_hmc<
       point_hmc,std::vector<double>,ubvector>);
    mpc.mct_hmc.stepper=new_stepper;

    mpc.mct_hmc.set_names_units(pnames_hmc,punits_hmc);
    
    mpc.mct_hmc.aff_inv=false;
    mpc.mct_hmc.verbose=3;
    mpc.mct_hmc.n_threads=1;
    mpc.mct_hmc.max_iters=2000;
    mpc.mct_hmc.prefix="mcmct_hmc";
    
    point_hmc ph=std::bind
      (std::mem_fn<int(size_t,const ubvector &,double &,
                       std::vector<double> &)>(&mcmc_para_class::hmc_point),
       &mpc,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,std::placeholders::_4);
    fill_hmc fh=std::bind
      (std::mem_fn<int(const ubvector &,double,std::vector<double> &,
                       std::vector<double> &)>(&mcmc_para_class::hmc_fill),
       &mpc,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3,std::placeholders::_4);
    
    vector<point_hmc> hmc_point_vec(1);
    vector<fill_hmc> hmc_fill_vec(1);
    hmc_point_vec[0]=ph;
    hmc_fill_vec[0]=fh;
    
    mpc.mct_hmc.mcmc_fill(2,low_hmc,high_hmc,hmc_point_vec,
                       hmc_fill_vec,data_vec_hmc);
    
    std::shared_ptr<o2scl::table_units<> > hmc_table=mpc.mct_hmc.get_table();

    tm.test_rel(vector_mean(hmc_table->get_nlines(),
                           (*hmc_table)["x"]),0.0,0.2,"hmc mean");
    tm.test_rel(vector_stddev(hmc_table->get_nlines(),
                           (*hmc_table)["x"]),1.0,0.2,"hmc mean");
    cout << endl;
    
  }

  if (true) {
    
    mcmc_para_emu<point_funct,fill_funct,std::vector<double>,
                  ubvector> mpe;

    vector<string> pnames={"x","x2"};
    vector<string> punits={"MeV","MeV^2"};
    mpe.set_names_units(pnames,punits);

    mpe.aff_inv=false;
    mpe.verbose=2;
    mpe.n_threads=1;
    mpe.max_iters=N/5;
    mpe.prefix="mpe";
    mpe.show_emu=1;
    mpe.def_stepper->step_fac[0]=10.0;

    // Set up the shared pointer to the interpolation object
    std::shared_ptr<interpm_krige_optim<>> iko(new interpm_krige_optim<>);
    mpe.emu.resize(1);
    mpe.emu[0]=iko;

    typedef const const_matrix_row_gen
      <o2scl::const_matrix_view_table<>> mat_x_row_t;
    
    // Setup the multidimensional covariance object
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(1);
    
    iko->verbose=2;
    iko->def_mmin.verbose=1;
    vector<double> len_list={0.1,0.3,1.0,3.0};
    vector<double> l10_list={-15,-13,-11};
    vector<vector<double>> ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);
    
    iko->set_covar(vmfrn,param_lists);
    
    mpe.emu_file="mcmct_0_out";

    mpe.test_emu_file="mcmc_test_emu.o2";
    
    mpe.mcmc_emu(1,low,high,gauss_vec,fill_vec,data_vec);
    
  }
  
  tm.report();
  
  return 0;
}

#endif
