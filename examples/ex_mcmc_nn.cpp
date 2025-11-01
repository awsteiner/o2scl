/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2024-2025, Andrew W. Steiner
  
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
#include <o2scl/kde_python.h>
#include <o2scl/gmm_python.h>
#include <o2scl/interpm_python.h>
#include <o2scl/interpm_krige.h>
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
      o2graph -set colbar 1 -set fig_dict "dpi=250" -create table3d \
      x grid:0,3,0.02 y grid:0,4,0.02 z \
      "if(10*(x-1)*(x-1.5)*(x-0.5)>y-2,exp(-(x-1)^2-(y-2)^2),-1)" \
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

  if (argc>=2 && ((std::string)argv[1]=="nn")) {
  
    // ─────────────────────────────────────────────────────────────────
    // Initial setup
  
    // Parameter limits
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

    // Set parameter names and units
    vector<string> pnames={"x","y"};
    vector<string> punits={"",""};
    vector<string> dnames={"cubic","class"};
    vector<string> dunits={"",""};
    mct.set_names_units(pnames,punits,dnames,dunits);

    // Set up HMC stepper
    shared_ptr<mcmc_stepper_hmc<point_funct,data_t,ubvector>> hmc_stepper
      (new mcmc_stepper_hmc<point_funct,data_t,ubvector>);
    mct.stepper=hmc_stepper;
    hmc_stepper->mom_step.resize(1);
    hmc_stepper->mom_step[0]=0.5;
    hmc_stepper->epsilon=0.04;

    /// MCMC parameters
    mct.store_pos_rets=true;
    mct.store_rejects=true;
    mct.n_retrain=0;
    mct.verbose=3;
    mct.n_threads=1;
    mct.max_iters=10000;
    mct.prefix="data/ex_mcmc_nn1";

    cout << "──────────────────────────────────────────────────────────"
         << endl;
    
    // ─────────────────────────────────────────────────────────────────
    // Create an initial guess using HMC
  
    mct.mcmc_fill(2,low_tf,high_tf,tf_vec,fill_vec,data_vec);
    cout << endl;

    // ─────────────────────────────────────────────────────────────────
    // Proposal distributions
  
    cout << "Training the proposal distributions." << endl;
  
    // Copy the table data to a tensor
    tensor<> ten_in, ten_in_copy;
    vector<size_t> in_size={mct.get_table()->get_nlines(),2};
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

    // Try three different kinds of proposal distributions

    // Normalizing flows
    ten_in_copy=ten_in;
    std::shared_ptr<nflows_python<ubvector>> nflows
      (new nflows_python<ubvector>);
    nflows->set_function("o2sclpy",ten_in_copy,"max_iter=20,verbose=1",
                         "nflows_nsf",1);
    cout << "Done training the nflows proposal distribution.\n" << endl;

    // KDE
    ten_in_copy=ten_in;
    std::shared_ptr<kde_python<ubvector>> kde
      (new kde_python<ubvector>);
    uniform_grid_log_end<double> ug(1.0e-3,1.0e3,99);
    vector<double> bw_array;
    ug.vector(bw_array);
    kde->set_function("o2sclpy","verbose=0","kde_sklearn",0);
    kde->set_data(ten_in_copy,bw_array);
                
    cout << "Done training the KDE proposal distribution.\n" << endl;

    // GMM
    ten_in_copy=ten_in;
    std::shared_ptr<gmm_python> gmm(new gmm_python);
    gmm->set_function("o2sclpy",4,ten_in_copy,"verbose=0","gmm_sklearn",0);
    gmm->get_python();
    std::shared_ptr<prob_dens_mdim_gmm<>> pdmg=gmm->get_gmm();
    cout << "Done training the GMM proposal distribution.\n" << endl;

    // Setup an array of shared pointers to the proposal distributions

    std::shared_ptr<mcmc_stepper_mh<point_funct,data_t,ubvector,ubmatrix,
                                    prob_cond_mdim_indep<>>> indep_stepper[3];
    for(size_t k=0;k<3;k++) {
      std::shared_ptr<mcmc_stepper_mh<point_funct,data_t,ubvector,ubmatrix,
                                      prob_cond_mdim_indep<>>> snew
        (new mcmc_stepper_mh<point_funct,data_t,ubvector,ubmatrix,
         prob_cond_mdim_indep<>>);
      indep_stepper[k]=snew;
    }
  
    indep_stepper[0]->proposal.resize(1);
    indep_stepper[0]->proposal[0].set_base(nflows);
    indep_stepper[1]->proposal.resize(1);
    indep_stepper[1]->proposal[0].set_base(kde);
    indep_stepper[2]->proposal.resize(1);
    indep_stepper[2]->proposal[0].set_base(pdmg);

    // ─────────────────────────────────────────────────────────────────
    // Test each proposal distribution using several MH steps.

    // First, sample the initial point from all three proposal
    // distributions.
    ubvector p[3], q[3];
    for(size_t k=0;k<3;k++) {
      p[k].resize(2);
      q[k].resize(2);
    }
    (*nflows)(p[0]);
    (*kde)(p[1]);
    (*pdmg)(p[2]);

    // Setup the testing tables
    table<> tprop[3];
    for(size_t k=0;k<3;k++) {
      tprop[k].line_of_names("px py qx qy lwp lwq q_prop qual");
    }
  
    // Compute the average deviation in the log_wgt for MH steps
    vector<double> qual(3);
    static const size_t N_test=2000;
    for(size_t i=0;i<N_test;i++) {
      (*nflows)(q[0]);
      (*kde)(q[1]);
      (*pdmg)(q[2]);
      for(size_t k=0;k<3;k++) {
        double q_prop=(indep_stepper[k])->proposal[0].log_metrop_hast
          (p[k],q[k]);

        double lw_p, lw_q;
        e.test_func(2,p[k],lw_p,data_vec[0]);
        e.test_func(2,q[k],lw_q,data_vec[0]);
      
        qual[k]+=fabs(lw_q-lw_p+q_prop);
        vector<double> line={p[k][0],p[k][1],q[k][0],q[k][1],
                             lw_p,lw_q,q_prop,fabs(lw_q-lw_p+q_prop)};
        tprop[k].line_of_data(line.size(),line);
      
        // "Accept" the step
        p[k][0]=q[k][0];
        p[k][1]=q[k][1];
      }
    
    }
    for(size_t k=0;k<3;k++) qual[k]/=N_test;
    std::cout << "Quality of proposal distributions: nflows: "
              << qual[0] << " KDE: " << qual[1] << " GMM: "
              << qual[2] << std::endl;

    // Select best proposal distribution and put it in index 0
    if (qual[1]<qual[2] && qual[1]<qual[0]) {
      cout << "Selecting KDE." << endl;
      std::swap(indep_stepper[0],indep_stepper[1]);
    } else if (qual[2]<qual[0] && qual[2]<qual[1]) {
      cout << "Selecting GMM." << endl;
      std::swap(indep_stepper[0],indep_stepper[2]);
    } else {
      cout << "Selecting nflows." << endl;
    }

    // Store proposal distribution tests to file
    hdf_file hf;
    hf.open_or_create("data/ex_mcmc_nn2.o2");
    hdf_output(hf,tprop[0],"prop_nflow");
    hdf_output(hf,tprop[1],"prop_KDE");
    hdf_output(hf,tprop[2],"prop_GMM");
    hf.close();
  
    // ─────────────────────────────────────────────────────────────────
    // Setup the three classifiers
  
    for(size_t k=0;k<3;k++) {
      std::shared_ptr<classify_python
                      <ubvector,ubvector_int,
                       o2scl::const_matrix_view_table<>,
                       o2scl::matrix_view_table<>>> cpnew
        (new classify_python
         <ubvector,ubvector_int,
         o2scl::const_matrix_view_table<>,
         o2scl::matrix_view_table<>>);
      mct.cl_list.push_back(cpnew);
    }
    mct.cl_list[0]->set_function("classify_sklearn_dtc","verbose=1",0);
    mct.cl_list[1]->set_function
      ("classify_sklearn_mlpc",
       ((std::string)"hlayers=[100,100],activation=")+
       "relu,verbose=1,max_iter=2000,"+
       "n_iter_no_change=40,tol=1.0e-5",0);
    mct.cl_list[2]->set_function("classify_sklearn_gnb","verbose=1",0);

    // ─────────────────────────────────────────────────────────────────
    // Setup the three emulators
  
    // TensorFlow neural network
    std::shared_ptr<interpm_python
                    <ubvector,
                     o2scl::const_matrix_view_table<>,
                     o2scl::matrix_view_table<>>> emu0
      (new interpm_python
       <ubvector,
       o2scl::const_matrix_view_table<>,
       o2scl::matrix_view_table<>>);
    emu0->set_function("interpm_tf_dnn",
                       ((std::string)"hlayers=[100,100],activations=")+
                       "[relu,relu],verbose=1,epochs=200",0);
    mct.emu.push_back(emu0);

    // Scikit-learn neural network
    std::shared_ptr<interpm_python
                    <ubvector,
                     o2scl::const_matrix_view_table<>,
                     o2scl::matrix_view_table<>>> emu1
      (new interpm_python
       <ubvector,
       o2scl::const_matrix_view_table<>,
       o2scl::matrix_view_table<>>);
    emu1->set_function("interpm_sklearn_mlpr",
                       ((std::string)"hlayers=[100,100],activation=")+
                       "relu,verbose=1,max_iter=400",0);
    mct.emu.push_back(emu1);

    // Gaussian process
    std::shared_ptr<interpm_krige_optim<>> iko(new interpm_krige_optim<>);
    typedef const const_matrix_row_gen
      <o2scl::const_matrix_view_table<>> mat_x_row_t;
  
    // Setup the covariance object
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(1);

    // Parameters for the Gaussian process optimization
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

    // Send the GP emulator to the MCMC object
    mct.emu.push_back(iko);
  
    // ─────────────────────────────────────────────────────────────────
    // Run the final MCMC
  
    // Set the MCMC parameters
    mct.stepper=indep_stepper[0];
    mct.use_emulator=true;
    mct.use_classifier=true;
    mct.n_retrain=0;
    mct.max_iters=200;
    //mct.max_iters=5000;
    mct.prefix="data/ex_mcmc_nn3";
    mct.n_threads=1;
    mct.verbose=3;
    mct.show_emu=1;
    mct.test_size=0.1;
    mct.max_emu_size=1000;
    mct.max_class_size=1000;
  
    // Set up the file for the emulator and classifier input. In this
    // example, they're the same, but they need not be.
    mct.emu_file="data/ex_mcmc_nn1_0_out";
    mct.emu_tname="markov_chain_0";
    mct.class_file="data/ex_mcmc_nn1_0_out";
    mct.class_tname="markov_chain_0";

    cout << "Calling mcmc_emu()" << endl;
    mct.mcmc_emu(2,low_tf,high_tf,tf_vec,fill_vec,data_vec);

  }

  tm.report();
  
  return 0;
}
