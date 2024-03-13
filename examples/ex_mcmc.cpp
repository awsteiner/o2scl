/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2017-2024, Andrew W. Steiner
  
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
/* Example: ex_mcmc.cpp
   ───────────────────────────────────────────────────────────────────
   An example which demonstrates the generation of an arbitrary
   distribution through Markov chain Monte Carlo. See "License 
   Information" section of the documentation for license information.
*/
#include <o2scl/mcmc_para.h>
#include <o2scl/vec_stats.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_io.h>
#include <o2scl/inte_qag_gsl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

/// Convenient typedefs

typedef boost::numeric::ublas::vector<double> ubvector;

typedef boost::numeric::ublas::matrix<double> ubmatrix;

typedef std::function<int(size_t,const ubvector &,double &,
			  std::vector<double> &)> point_funct;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  std::vector<double> &)> fill_funct;
/// The MCMC object
mcmc_para_table<point_funct,fill_funct,std::vector<double>,ubvector> mct;

/** \brief A demonstration class for the MCMC example. This example
    could have been written with global functions, but we put them
    in a class to show how it would work in that case.
 */
class exc {

public:
  
  /** \brief A one-dimensional bimodal distribution

      Here, the variable 'log_weight' stores the natural logarithm of
      the objective function based on the parameters stored in \c pars.
      The object 'dat' stores any auxillary quantities which can be
      computed at every point in parameter space.
  */
  int bimodal(size_t nv, const ubvector &pars, double &log_weight,
	      std::vector<double> &dat) {
    
    double x=pars[0];
    log_weight=log(exp(-x*x)*(sin(x-1.4)+1.0));
    dat[0]=x*x;
    return 0;
  }

  /** \brief Add auxillary quantities to 'line' so they can be
      stored in the table
  */
  int fill_line(const ubvector &pars, double log_weight, 
		std::vector<double> &line, std::vector<double> &dat) {
    line.push_back(dat[0]);
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

  ubvector low_bimodal(1), high_bimodal(1);
  low_bimodal[0]=-5.0;
  high_bimodal[0]=5.0;

  // Function objects for the MCMC object
  point_funct bimodal_func=std::bind
    (std::mem_fn<int(size_t,const ubvector &,double &,
		     std::vector<double> &)>(&exc::bimodal),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  fill_funct fill_func=std::bind
    (std::mem_fn<int(const ubvector &,double,std::vector<double> &,
		     std::vector<double> &)>(&exc::fill_line),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);

  // Create function object vectors
  vector<point_funct> bimodal_vec;
  bimodal_vec.push_back(bimodal_func);
  vector<fill_funct> fill_vec;
  fill_vec.push_back(fill_func);

  // Create and allocate data objects
  vector<std::vector<double> > data_vec(2);
  data_vec[0].resize(1);
  data_vec[1].resize(1);

  // Compute exact value of <x^2>. The function format is a bit
  // different so we use lambda expressions to construct the
  // functions for the integrators. 
  inte_qag_gsl<> iqg;
  funct f=[e,data_vec](double x) mutable -> double {
    ubvector u(1); double lw; u[0]=x; 
    e.bimodal(1,u,lw,data_vec[0]); return exp(lw); };
  funct fx2=[e,data_vec](double x) mutable -> double {
    ubvector u(1); double lw; u[0]=x; 
    e.bimodal(1,u,lw,data_vec[0]); return data_vec[0][0]*exp(lw); };
  cout << iqg.integ(fx2,low_bimodal[0],high_bimodal[0]) << " "
       << iqg.integ(f,low_bimodal[0],high_bimodal[0]) << endl;
  double exact=iqg.integ(fx2,low_bimodal[0],high_bimodal[0])/
    iqg.integ(f,low_bimodal[0],high_bimodal[0]);
  cout << "exact: " << exact << endl;
  
  cout << "──────────────────────────────────────────────────────────"
       << endl;
  cout << "Plain MCMC example with a bimodal distribution:" << endl;
    
  // Set parameter names and units
  vector<string> pnames={"x","x2"};
  vector<string> punits={"",""};
  mct.set_names_units(pnames,punits);

  // Set MCMC parameters
  mct.step_fac=2.0;
  mct.max_iters=20000;
  mct.prefix="ex_mcmc";
  mct.table_prealloc=mct.max_iters/3;
  
  // Perform MCMC
  mct.mcmc_fill(1,low_bimodal,high_bimodal,bimodal_vec,fill_vec,
                data_vec);

  // Output acceptance and rejection rate
  cout << "n_accept, n_reject: " << mct.n_accept[0] << " "
       << mct.n_reject[0] << endl;

  // Get results
  shared_ptr<table_units<> > t=mct.get_table();

  // Remove empty rows from table
  t->delete_rows_func("mult<0.5");
  
  // Compute autocorrelation length and effective sample size
  std::vector<double> ac, ftom;
  o2scl::vector_autocorr_vector_mult(t->get_nlines(),
				     (*t)["x2"],(*t)["mult"],ac);
  size_t ac_len=o2scl::vector_autocorr_tau(ac,ftom);
  cout << "Autocorrelation length, effective sample size: "
       << ac_len << " " << t->get_nlines()/ac_len << endl;

  // Create a set of fully independent samples
  t->new_column("N");
  for(size_t j=0;j<t->get_nlines();j++) {
    t->set("N",j,j);
  }
  std::string func=((std::string)"N%")+szttos(ac_len)+">0.5";
  t->delete_rows_func(func);
  for(size_t j=0;j<t->get_nlines();j++) {
    t->set("N",j,j);
  }

  // Write these samples to a file
  hdf_file hf;
  hf.open_or_create("ex_mcmc.o2");
  hdf_output(hf,*t,"mcmc");
  hf.close();

  // Use the independent samples to compute the final integral and
  // compare to the exact result
  double avg=vector_mean((*t)["x2"]);
  double std=vector_stddev((*t)["x2"]);
  tm.test_rel(avg,exact,10.0*std/sqrt(t->get_nlines()),"ex_mcmc");
  cout << "avg. from MCMC, exact avg., diff., std. from MCMC, "
       << "unc. in mean times 10:\n  "
       << avg << " " << exact << " " << fabs(avg-exact) << " "
       << std << " " << 10.0*std/sqrt(t->get_nlines()) << endl;
  
  tm.report();
  
  return 0;
}
