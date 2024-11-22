/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2024, Andrew W. Steiner
  
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
// sphinx-example-start
/* Example: ex_mcmc_nn.cpp
   ───────────────────────────────────────────────────────────────────
   An example which demonstrates ..
*/
#include <o2scl/mcmc_para.h>
#include <o2scl/vec_stats.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_io.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/set_openmp.h>
#include <o2scl/nflows_python.h>
#include <o2scl/interpm_python.h>
#include <o2scl/classify_python.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

/// Convenient typedefs

typedef boost::numeric::ublas::vector<double> ubvector;

typedef boost::numeric::ublas::matrix<double> ubmatrix;

typedef ubvector data_t;

typedef std::function<int(size_t,const ubvector &,double &,
			  data_t &)> point_funct;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  data_t &)> fill_funct;

class mcmc_stepper_mh_record :
  public mcmc_stepper_mh<point_funct,data_t,ubvector,
                         ubmatrix,prob_cond_mdim_indep<>> {
  
  public:

  std::vector<double> vw_next, vq_next;

  virtual ~mcmc_stepper_mh_record() {
  }
  
  virtual void step(size_t i_thread, size_t n_params, point_funct &f,
                    const ubvector &current, ubvector &next,
                    double w_current, double &w_next,
                    const ubvector &low, const ubvector &high,
                    int &func_ret, bool &accept, data_t &dat,
                    rng<> &r, int verbose) {
    
    // Use proposal distribution and compute associated weight
    double q_prop=proposal[i_thread % proposal.size()].log_metrop_hast
      (current,next);
    
    accept=false;
    
    func_ret=success;
    this->check_bounds(i_thread,n_params,next,low,high,
                       func_ret,verbose);
    if (func_ret!=this->mcmc_skip) {
      func_ret=f(n_params,next,w_next,dat);
    } 

    vw_next.push_back(w_next);
    double q_next=proposal[i_thread %
                           proposal.size()].log_pdf(current,next);
    vq_next.push_back(q_next);
    
    if (func_ret==success) {
      double rand=r.random();
      
      // Metropolis-Hastings algorithm
      if (rand<exp(w_next-w_current+q_prop)) {
        accept=true;
      }
    }
    
    return;
  }
  
};

/// The MCMC object
mcmc_para_emu<point_funct,fill_funct,data_t,ubvector> mct;

/** \brief A demonstration class for the MCMC example. This example
    could have been written with global functions, but we put them
    in a class to show how it would work in that case.
 */
class exc {

public:
  
  /** \brief A two-dimesional probability distribution

      Here, the variable 'log_weight' stores the natural logarithm of
      the objective function based on the parameters stored in \c pars.
      The object 'dat' stores any auxillary quantities which can be
      computed at every point in parameter space.
  */
  int test_func(size_t nv, const ubvector &pars, double &log_weight,
                data_t &dat) {
    
    double x=pars[0];
    double y=pars[1];

    /*
      plot with:
      o2graph -set fig_dict "dpi=250" -create table3d \
      x grid:0,3,0.02 y grid:0,4,0.02 z \
      "if(10*(x-1)*(x-1.5)*(x-0.5)>y-2,-(x-1)^2-(y-2)^2,0)" \
      -den-plot z -show
    */
    
    // Add a complicated nonlinear constraint
    double cubic=10.0*(x-1.0)*(x-1.5)*(x-0.5);
    
    // A simple two-dimensional Gaussian
    log_weight=-pow(x-1.0,2.0)-pow(y-2.0,2.0);
    
    dat[0]=cubic;
    // The constraint for the classifier
    dat[1]=1.0;

    if (cubic<y-2.0) {
      dat[1]=0.0;
      return 1;
    }
    
    return 0;
  }

  /** \brief Add auxillary quantities to 'line' so they can be
      stored in the table
  */
  int fill_line(const ubvector &pars, double log_weight, 
		std::vector<double> &line, data_t &dat) {
    line.push_back(dat[0]);
    line.push_back(dat[1]);
    return 0;
  }

  exc() {
  }
  
};

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);
  
  exc e;
  
  test_mgr tm;
  tm.set_output_level(1);

  // Parameter limits and initial points

  ubvector low_tf(2), high_tf(2);
  low_tf[0]=-5.0;
  high_tf[0]=10.0;
  low_tf[1]=-5.0;
  high_tf[1]=10.0;

  // Function objects for the MCMC object
  point_funct tf_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,double &,
		     data_t &)>(&exc::test_func),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  fill_funct fill_func=std::bind
    (std::mem_fn<int(const ubvector &,double,std::vector<double> &,
		     data_t &)>(&exc::fill_line),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);

  // Create function object vectors
  vector<point_funct> tf_vec;
  tf_vec.push_back(tf_func);
  vector<fill_funct> fill_vec;
  fill_vec.push_back(fill_func);

  // Create and allocate data objects
  vector<data_t> data_vec(2);
  data_vec.resize(2);
  data_vec[0].resize(2);
  data_vec[1].resize(2);

  cout << "──────────────────────────────────────────────────────────"
       << endl;
    
  // Set parameter names and units
  vector<string> pnames={"x","y","cubic","const"};
  vector<string> punits={"","","",""};
  mct.set_names_units(pnames,punits);

  // Run an initial HMC simulation
  
  shared_ptr<mcmc_stepper_hmc<point_funct,data_t,ubvector>> hmc_stepper
    (new mcmc_stepper_hmc<point_funct,data_t,ubvector>);
  mct.stepper=hmc_stepper;
  hmc_stepper->mom_step.resize(1);
  hmc_stepper->mom_step[0]=0.5;
  hmc_stepper->epsilon=0.02;

  mct.store_pos_rets=true;
  mct.store_rejects=true;
  mct.n_retrain=0;
  mct.verbose=3;
  mct.n_threads=1;
  mct.max_iters=500;
  mct.prefix="ex_mcmc_nn1";

  // Run MCMC
  mct.mcmc_fill(2,low_tf,high_tf,tf_vec,fill_vec,data_vec);
  cout << endl;

  // Copy the table data to a tensor for use in the proposal
  // distribution
  tensor<> ten_in;
  vector<size_t> in_size={mct.get_table()->get_nlines(),2};
  ten_in.resize(2,in_size);
  vector<size_t> class_size={mct.get_table()->get_nlines(),1};
  ten_in.resize(2,in_size);
  
  for(size_t i=0;i<mct.get_table()->get_nlines();i++) {
    // Only select lines which have mult greater than 0.5
    // for training the proposal distribution.
    if (mct.get_table()->get("mult",i)>0.5) {
      vector<size_t> ix;
      ix={i,0};
      ten_in.get(ix)=mct.get_table()->get("x",i);
      ix={i,1};
      ten_in.get(ix)=mct.get_table()->get("y",i);
    }
  }

  // Train the normalizing flows probability distribution
  cout << "Training the proposal distribution." << endl;
  std::shared_ptr<nflows_python<>> np(new nflows_python<>);
  np->set_function("o2sclpy",ten_in,"max_iter=200,verbose=1",
                   "nflows_nsf",1);
  cout << "Done training the proposal distribution.\n" << endl;

  // Sample the proposal distribution and store to a file
  table<> tprop;
  tprop.line_of_names("x y");
  for(size_t i=0;i<200;i++) {
    vector<double> p(2);
    (*np)(p);
    tprop.line_of_data(2,p);
  }
  hdf_file hf;
  hf.open_or_create("ex_mcmc_nn2.o2");
  hdf_output(hf,tprop,"tprop");
  hf.close();

  // Set up the classifier
  cout << "Training the classifier." << endl;
  std::shared_ptr<classify_python
                  <ubvector,ubvector_int,
                   o2scl::const_matrix_view_table<>,
                   o2scl::matrix_view_table<>>> cp
    (new classify_python
     <ubvector,ubvector_int,
     o2scl::const_matrix_view_table<>,
     o2scl::matrix_view_table<>>);
  cp->set_functions("classify_sklearn_mlpc",
                    ((std::string)"hlayers=[100,100],activation=")+
                    "relu,verbose=1,max_iter=2000,"+
                    "n_iter_no_change=40,tol=1.0e-5",1);
  
  o2scl::const_matrix_view_table<> cp_in(*(mct.get_table()),{"x","y"});
  o2scl::matrix_view_table<> cp_out(*(mct.get_table()),{"const"});
  
  cp->set_data(2,1,mct.get_table()->get_nlines(),cp_in,cp_out);
  cout << "Done training the classifier.\n" << endl;

  // Test the classifier by giving it random points
  table<> tclass;
  rng<> r;
  tclass.line_of_names("x y c");
  for(size_t i=0;i<200;i++) {
    vector<double> in(2);
    vector<int> out(1);
    in[0]=r.random()*3.0;
    in[1]=r.random()*4.0;
    cp->eval_std_vec(in,out);
    vector<double> line={in[0],in[1],(double)out[0]};
    tclass.line_of_data(3,line);
  }
  hf.open("ex_mcmc_nn2.o2",true);
  hdf_output(hf,tclass,"tclass");
  hf.close();
  
  // Before we train the emulator, remove the points which don't satisfy
  // the constraint
  mct.get_table()->delete_rows_func("const<0.5");
  
  // Set up the emulator
  cout << "Training the emulator." << endl;
  std::shared_ptr<interpm_python
                  <ubvector,
                   o2scl::const_matrix_view_table<>,
                   o2scl::matrix_view_table<>>> ip
    (new interpm_python
     <ubvector,
     o2scl::const_matrix_view_table<>,
     o2scl::matrix_view_table<>>);
  ip->set_functions("interpm_tf_dnn",
                    ((std::string)"hlayers=[100,100],activations=")+
                    "['relu','relu'],verbose=1,epochs=2000",1);
  
  o2scl::const_matrix_view_table<> ip_in(*(mct.get_table()),{"x","y"});
  o2scl::matrix_view_table<> ip_out(*(mct.get_table()),{"log_wgt"});
  
  ip->set_data(2,1,mct.get_table()->get_nlines(),ip_in,ip_out);
  cout << "Done training the emulator\n." << endl;
  
  // Test the emulator by giving it random points
  table<> temu;
  temu.line_of_names("x y lw");
  for(size_t i=0;i<200;i++) {
    vector<double> in(2);
    vector<double> out(1);
    in[0]=r.random()*3.0;
    in[1]=r.random()*4.0;
    ip->eval_std_vec(in,out);
    vector<double> line={in[0],in[1],out[0]};
    temu.line_of_data(3,line);
  }
  hf.open("ex_mcmc_nn2.o2",true);
  hdf_output(hf,temu,"temu");
  hf.close();
  exit(-1);
  
#ifdef O2SCL_NEVER_DEFINED
  
  // Use independence Metropolis-Hastings
  mct.stepper=shared_ptr<mcmc_stepper_base<point_funct,data_t,ubvector>>
    (new mcmc_stepper_mh<point_funct,data_t,ubvector,ubmatrix,
     prob_cond_mdim_indep<>>);
  mct.stepper->proposal.resize(1);
  mct.stepper->proposal[0].set_base(np);
  
  // Set the MCMC parameters
  mct.use_classifier=true;
  mct.max_iters=20000;
  mct.prefix="ex_mcmc_nn2";
  mct.n_threads=1;
  mct.verbose=3;
  
  // Perform MCMC
  mct.mcmc_fill(1,low_tf,high_tf,tf_vec,fill_vec,
                data_vec);

  // Output acceptance and rejection rate
  cout << "n_accept, n_reject: " << mct.n_accept[0] << " "
       << mct.n_reject[0] << endl;

  // Get results
  shared_ptr<table_units<> > t=mct.get_table();

  // Remove empty rows from table.
  t->delete_rows_func("mult<0.5");

  // Compute the autocorrelation length
  ubvector ac, ftom;
  o2scl::vector_autocorr_vector_mult(t->get_nlines(),
                                     (*t)["x2"],(*t)["mult"],ac);
  size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
  
  // Create a separate table of statistically independent samples
  table_units<> indep;
  copy_table_thin_mcmc(ac_len,*t,indep,"mult");
  
  cout << "Autocorrelation length, effective sample size: "
       << ac_len << " " << indep.get_nlines() << endl;
  
  // Write these samples to a file
  hf.open_or_create("ex_mcmc_kde.o2");
  hdf_output(hf,*t,"mcmc");
  hdf_output(hf,indep,"indep");
  hf.setd_vec("q_next",local_stepper->vq_next);
  hf.setd_vec("w_next",local_stepper->vw_next);
  hf.close();

  // Compute the average of the correlated samples for comparison
  double avg2=vector_mean(t->get_nlines(),(*t)["x2"]);
  cout << "Average of correlated samples: " << avg2 << endl;
  
  // Use the independent samples to compute the final integral and
  // compare to the exact result. Note that we must specify the
  // number of elements in the vector, indep["x2"], because the
  // table_units object often has space at the end to add extra rows.
  
  double avg=vector_mean(indep.get_nlines(),indep["x2"]);
  double std=vector_stddev(indep.get_nlines(),indep["x2"]);
  cout << "Average and std. dev. of uncorrelated samples: "
       << avg << " " << std << endl;
  cout << "Absolute difference: " << fabs(avg-exact) << endl;
  cout << "Uncertainty in the average: "
       << std/sqrt(indep.get_nlines()) << endl;
  tm.test_rel(avg,exact,10.0*std/sqrt(indep.get_nlines()),"ex_mcmc");
  
  tm.report();

#endif
  
  return 0;
}
