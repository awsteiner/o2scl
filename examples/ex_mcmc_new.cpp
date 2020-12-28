/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2021, Andrew W. Steiner
  
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
#include <o2scl/mcmc_para.h>
#include <o2scl/vec_stats.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_io.h>
#include <o2scl/cubature.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

/// Several convenient typedefs

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

typedef std::function<int(size_t,const ubvector &,double &,
			  std::array<double,2> &)> point_funct;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  std::array<double,2> &)> fill_funct;

/// The MCMC object
mcmc_para_table<point_funct,fill_funct,std::array<double,2>,ubvector> mct;

class exc {

public:
  
  /** \brief A simple two-dimensional Gaussian

      Here, the variable 'log_weight' stores the natural logarithm of
      the objective function based on the parameters stored in \c pars.
      The object 'dat' stores any auxillary quantities which can be
      computed at every point in parameter space.
  */
  int gauss2d(size_t nv, const ubvector &pars, double &log_weight,
	    std::array<double,2> &dat) {
  
    log_weight=-(pars[0]-0.2)*(pars[0]-0.2)-
      (pars[1]-0.5)*(pars[1]-0.5);
    
    dat[0]=pars[0]*pars[0];
    dat[1]=pars[0]*pars[0]*pars[1]*pars[1];
    return 0;
  }

  /** \brief A one-dimensional bimodal distribution

      Here, the variable 'log_weight' stores the natural logarithm of
      the objective function based on the parameters stored in \c pars.
      The object 'dat' stores any auxillary quantities which can be
      computed at every point in parameter space.
  */
  int bimodal(size_t nv, const ubvector &pars, double &log_weight,
	     std::array<double,2> &dat) {
    
    double x=pars[0];
    log_weight=log(exp(-x*x)*(sin(x-1.4)+1.0));
    dat[0]=0.0;
    dat[1]=1.0;
    return 0;
  }

  /** \brief Add auxillary quantities to 'line' so they can be
      stored in the table
  */
  int fill_line(const ubvector &pars, double log_weight, 
		std::vector<double> &line, std::array<double,2> &dat) {
    line.push_back(dat[0]);
    line.push_back(dat[1]);
    return 0;
  }

  exc() {
  }
  
};

/** Function for exact integration of the Gaussian
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
  
  exc e;
  
  test_mgr tm;
  tm.set_output_level(1);

  // Parameter limits and initial points

  ubvector low_gauss2d(2), high_gauss2d(2), init_gauss2d(2);
  low_gauss2d[0]=-1.5;
  low_gauss2d[1]=-1.5;
  high_gauss2d[0]=1.5;
  high_gauss2d[1]=1.5;
  init_gauss2d[0]=0.2;
  init_gauss2d[1]=0.5;
  
  ubvector low_bimodal(1), high_bimodal(1);
  low_bimodal[0]=-5.0;
  high_bimodal[0]=5.0;

  // Use cubature to compute integrals
  std::vector<double> dlow(2), dhigh(2);
  dlow[0]=-1.5;
  dlow[1]=-1.5;
  dhigh[0]=1.5;
  dhigh[1]=1.5;
  std::vector<double> dres(3), derr(3);
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

  // Function objects for the MCMC object
  point_funct gauss2d_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,double &,
		     std::array<double,2> &)>(&exc::gauss2d),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  point_funct bimodal_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,double &,
		     std::array<double,2> &)>(&exc::bimodal),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  fill_funct fill_func=std::bind
    (std::mem_fn<int(const ubvector &,double,std::vector<double> &,
		     std::array<double,2> &)>(&exc::fill_line),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);

  vector<point_funct> gauss2d_vec;
  gauss2d_vec.push_back(gauss2d_func);
  vector<point_funct> bimodal_vec;
  bimodal_vec.push_back(bimodal_func);
  vector<fill_funct> fill_vec;
  fill_vec.push_back(fill_func);

  if (true) {

    cout << "----------------------------------------------------------"
	 << endl;
    cout << "Plain MCMC example with a bimodal distribution:" << endl;
    
    // Set parameter names and units
    vector<string> pnames={"x","0","1"};
    vector<string> punits={"","",""};
    mct.set_names_units(pnames,punits);

    // Set MCMC parameters
    mct.step_fac=2.0;
    mct.max_iters=100000;

    // Perform MCMC
    mct.mcmc(1,low_bimodal,high_bimodal,bimodal_vec,fill_vec);

    // Compute autocorrelation length and effective sample size
    shared_ptr<table_units<> > t=mct.get_table();
    std::vector<double> ac, ftom;
    o2scl::vector_autocorr_vector((*t)["log_wgt"],ac);
    size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
    cout << "Autocorrelation length, effective sample size: "
	 << ac_len << " " << t->get_nlines()/ac_len << endl;
    //mct.reblock(t->get_nlines()/ac_len);
    hist h;
    h.from_table(*t,"x",40);
    hdf_file hf;
    hf.open_or_create("ex_mcmc_bimodal.o2");
    hdf_output(hf,h,"hist");
    hdf_output(hf,*t,"table");
    hf.close();

  }
  
  if (true) {
    
    cout << "----------------------------------------------------------"
	 << endl;
    cout << "Plain MCMC example with Gaussian distribution:" << endl;

    // Parameter names and units
    vector<string> pnames={"x0","x1","0","1"};
    vector<string> punits={"","","",""};
    mct.set_names_units(pnames,punits);

    // MCMC parameters
    mct.step_fac=2.0;
    mct.max_iters=100000;

    mct.verbose=2;
    mct.mcmc(2,low_gauss2d,high_gauss2d,gauss2d_vec,fill_vec);

    // Analyze results
    cout << "Steps accepted, rejected: " << mct.n_accept[0] << " "
	 << mct.n_reject[0] << endl;
    std::vector<double> x02, x12;
    shared_ptr<table_units<> > t;
    t=mct.get_table();
    t->new_column("weight");
    t->delete_rows_func("mult<0.5");
    for(size_t j=0;j<t->get_nlines();j++) {
      t->set("weight",j,exp(t->get("log_wgt",j)));
      for(size_t k=0;k<t->get("mult",j);k++) {
	x02.push_back(t->get("x0",j));
	x12.push_back(t->get("x1",j));
      }
    }
    size_t ss;
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector(x02,ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x02.size()/ac_len << endl;
      ss=x02.size()/ac_len;
    }
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector_mult((*t)["x0"],(*t)["mult"],ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x02.size()/ac_len << endl;
      ss=x02.size()/ac_len;
    }
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector(x12,ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x12.size()/ac_len << endl;
      if (x12.size()/ac_len<ss) ss=x12.size()/ac_len;
    }
    table_units<> t2(ss);
    t2.line_of_names("x0 x1 x0b x1b");
    size_t ac_len=x02.size()/ss;
    if (true) {
      // Block averaging
      for(size_t j=0;j<ss;j++) {
	t2.set("x0b",j,0.0);
	t2.set("x1b",j,0.0);
	for(size_t k=0;k<ac_len;k++) {
	  t2.set("x0b",j,t2.get("x0b",j)+x02[j*ac_len+k]/ac_len);
	  t2.set("x1b",j,t2.get("x1b",j)+x12[j*ac_len+k]/ac_len);
	}
      }
    }
    if (true) {
      // Thinning
      for(size_t j=0;j<ss;j++) {
	t2.set("x0",j,x02[j*ac_len]);
	t2.set("x1",j,x12[j*ac_len]);
      }
    }
    hdf_file hf2;
    hf2.open_or_create("ex_mcmc_gaussian.o2");
    hf2.setd_vec("x0",x02);
    hdf_output(hf2,t2,"table2");
    hf2.close();
    if (false) {
      {
	std::vector<double> ac, ftom;
	o2scl::vector_autocorr_vector((*t)["log_wgt"],ac);
	size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
	cout << "ac_len,samp_size: " << ac_len << " "
	     << t->get_nlines()/ac_len << endl;
      }
      {
	std::vector<double> ac, ftom;
	o2scl::vector_autocorr_vector((*t)["weight"],ac);
	size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
	cout << "ac_len,samp_size: " << ac_len << " "
	     << t->get_nlines()/ac_len << endl;
      }
      //mct.reblock(t->get_nlines()/ac_len);
      hist h;
      h.from_table(*t,"x0",40);
      hdf_file hf;
      hf.open_or_create("ex_mcmc_gaussian.o2");
      hf.setd_vec("x0",x02);
      hdf_output(hf,h,"hist");
      hdf_output(hf,*t,"table");
      hf.close();
    }
  }

  if (true) {
    cout << "----------------------------------------------------------"
	 << endl;
    cout << "MCMC with affine-invariant sampling:" << endl;
    vector<string> pnames={"x0","x1","0","1"};
    vector<string> punits={"","","",""};
    mct.set_names_units(pnames,punits);
    mct.step_fac=2.0;
    mct.aff_inv=true;
    mct.n_walk=5;
    mct.max_iters=100000;
    std::vector<double> x0[5], x1[5];
    shared_ptr<table_units<> > t;
    mct.mcmc(2,low_gauss2d,high_gauss2d,gauss2d_vec,fill_vec);
    cout << "Steps accepted, rejected: " << mct.n_accept[0] << " "
	 << mct.n_reject[0] << endl;
    t=mct.get_table();
    for(size_t j=0;j<t->get_nlines();j++) {
      size_t walk=((size_t)(t->get("walker",j)+1.0e-6));
      for(size_t k=0;k<t->get("mult",j);k++) {
	x0[walk].push_back(t->get("x0",j));
	x1[walk].push_back(t->get("x1",j));
      }
    }
    size_t ac_len_all;
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector(x0[0],ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      ac_len_all=ac_len;
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x0[0].size()/ac_len << endl;
    }
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector(x1[0],ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      if (ac_len_all<ac_len) ac_len_all=ac_len;
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x1[0].size()/ac_len << endl;
    }
    table_units<> t2;
    t2.line_of_names("x0 x1");
    for(size_t k=0;k<mct.n_walk;k++) {
      for(size_t j=0;j<x0[k].size();j+=ac_len_all) {
	double line[2]={x0[k][j],x0[k][j]};
	t2.line_of_data(2,line);
      }
    }
    hdf_file hf2;
    hf2.open_or_create("ex_mcmc_affinv.o2");
    hdf_output(hf2,t2,"table2");
    hf2.close();
  }
  
  if (true) {
    cout << "----------------------------------------------------------"
	 << endl;
    cout << "MCMC with perfect Gaussian proposal distribution:"
	 << endl;
    vector<string> pnames={"x0","x1","0","1"};
    vector<string> punits={"","","",""};

    mct.set_names_units(pnames,punits);
    mct.step_fac=2.0;
    mct.max_iters=40000;

    vector<prob_cond_mdim_indep<ubvector> *> pcmi;
    ubvector peak(2);
    peak[0]=0.2;
    peak[1]=0.5;
    ubmatrix covar(2,2);
    covar(0,0)=0.75;
    covar(1,1)=0.75;
    covar(0,1)=0.0;
    covar(1,0)=0.0;
    prob_dens_mdim_gaussian<ubvector> pdmg(2,peak,covar);
    
    pcmi.resize(1);
    pcmi[0]=new prob_cond_mdim_indep<ubvector>(pdmg);
    ubvector step(2);
    mct.set_proposal_ptrs(pcmi);
    
    std::vector<double> x02, x12;
    shared_ptr<table_units<> > t;
    low_gauss2d[0]=-5.0;
    low_gauss2d[1]=-5.0;
    high_gauss2d[0]=5.0;
    high_gauss2d[1]=5.0;
    std::cout << "Going to mcmc." << std::endl;
    mct.mcmc(2,low_gauss2d,high_gauss2d,gauss2d_vec,fill_vec);
    // Almost all steps should be accepted because the proposal
    // distribution is so close to the true distribution.
    // Only steps near the boundaries are rejected.
    cout << "Steps accepted, rejected: " << mct.n_accept[0] << " "
	 << mct.n_reject[0] << endl;
    t=mct.get_table();
    t->new_column("weight");
    for(size_t j=0;j<t->get_nlines();j++) {
      t->set("weight",j,exp(t->get("log_wgt",j)));
      for(size_t k=0;k<t->get("mult",j);k++) {
	x02.push_back(t->get("x0",j));
	x12.push_back(t->get("x1",j));
      }
    }
    size_t ss;
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector(x02,ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x02.size()/ac_len << endl;
      ss=x02.size()/ac_len;
    }
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector(x12,ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x12.size()/ac_len << endl;
      if (x12.size()/ac_len<ss) ss=x12.size()/ac_len;
    }
    table_units<> t2(ss);
    t2.line_of_names("x0 x1");
    size_t ac_len=x02.size()/ss;
    // Thinning
    for(size_t j=0;j<ss;j++) {
      t2.set("x0",j,x02[j*ac_len]);
      t2.set("x1",j,x12[j*ac_len]);
    }
    hdf_file hf2;
    hf2.open_or_create("ex_mcmc_pdgauss.o2");
    hdf_output(hf2,t2,"table2");
    hf2.close();
    exit(-1);
  }

  if (true) {
    // Try the fixed_step conditional probability
    vector<string> pnames={"x0","x1","0","1"};
    vector<string> punits={"","","",""};

    mct.set_names_units(pnames,punits);
    mct.step_fac=2.0;
    mct.max_iters=100000;

    vector<prob_cond_mdim_fixed_step<ubvector> > pcmrw;
    pcmrw.resize(1);
    ubvector step(2);
    step[0]=(high_gauss2d[0]-low_gauss2d[0])/mct.step_fac;
    step[1]=(high_gauss2d[1]-low_gauss2d[1])/mct.step_fac;
    pcmrw[0].set(step,low_gauss2d,high_gauss2d);
    mct.set_proposal(pcmrw);
    
    std::vector<double> x02, x12;
    shared_ptr<table_units<> > t;
    mct.mcmc(2,low_gauss2d,high_gauss2d,gauss2d_vec,fill_vec);
    cout << "Steps accepted, rejected: " << mct.n_accept[0] << " "
	 << mct.n_reject[0] << endl;
    t=mct.get_table();
    t->new_column("weight");
    for(size_t j=0;j<t->get_nlines();j++) {
      t->set("weight",j,exp(t->get("log_wgt",j)));
      for(size_t k=0;k<t->get("mult",j);k++) {
	x02.push_back(t->get("x0",j));
	x12.push_back(t->get("x1",j));
      }
    }
    size_t ss;
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector(x02,ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x02.size()/ac_len << endl;
      ss=x02.size()/ac_len;
    }
    {
      std::vector<double> ac, ftom;
      o2scl::vector_autocorr_vector(x12,ac);
      size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
      cout << "ac_len,samp_size: " << ac_len << " "
	   << x12.size()/ac_len << endl;
      if (x12.size()/ac_len<ss) ss=x12.size()/ac_len;
    }
    table_units<> t2(ss);
    t2.line_of_names("x0 x1");
    size_t ac_len=x02.size()/ss;
    // Thinning
    for(size_t j=0;j<ss;j++) {
      t2.set("x0",j,x02[j*ac_len]);
      t2.set("x1",j,x12[j*ac_len]);
    }
    hdf_file hf2;
    hf2.open_or_create("ex_mcmc_propdist.o2");
    hdf_output(hf2,t2,"table2");
    hf2.close();
  }

  if (false) {
    
    // Table column names and units. We must specify first the names for
    // the parameters first and then the names of the auxillary
    // parameters.
    vector<string> pnames={"x0","x1","x0sq","x0sq_x1sq"};
    vector<string> punits={"MeV","MeV","MeV^2","MeV^4"};
    mct.set_names_units(pnames,punits);

    // MCMC with a random walk of a fixed length
    cout << "MCMC with random walk:\n" << endl;
    // This step factor is chosen to give approximately equal number of
    // accept/reject steps but will be different for different problems
    mct.step_fac=3.0;
    mct.max_iters=1000;
    mct.mcmc(2,low_gauss2d,high_gauss2d,gauss2d_vec,fill_vec);

    // Get a pointer to the results table
    shared_ptr<table_units<> > t=mct.get_table();

    // Output table and other information
    cout << "n_accept, n_reject, table lines: "
	 << mct.n_accept[0] << " " << mct.n_reject[0] << " "
	 << t->get_nlines() << endl;
    cout << "    i mult        log_wgt     x0          x1          "
	 << "x0sq        x0sq_x1sq" << endl;
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
    cout << "MCMC with affine-invariant sampling:\n" << endl;

    mct.aff_inv=true;
    mct.n_walk=10;
    // This step factor is chosen to give approximately equal number
    // of accept/reject steps but will be different for
    // different problems
    mct.step_fac=5.0;
    mct.mcmc(2,low_gauss2d,high_gauss2d,gauss2d_vec,fill_vec);

    // Output table and other information
    cout << "n_accept, n_reject, table lines: "
	 << mct.n_accept[0] << " " << mct.n_reject[0] << " "
	 << t->get_nlines() << endl;
    cout << "    i mult        log_wgt     x0          x1          "
	 << "x0sq        x0sq_x1sq" << endl;
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

    // Perform block averaging
    mct.reblock(40);

    // Compute and test the average values
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

  }
  
  tm.report();
  
  return 0;
}
// End of example
