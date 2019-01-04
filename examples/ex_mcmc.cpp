/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2019, Andrew W. Steiner
  
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

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

/// Convenient typedefs

typedef boost::numeric::ublas::vector<double> ubvector;

typedef boost::numeric::ublas::matrix<double> ubmatrix;

typedef std::function<int(size_t,const ubvector &,double &,
			  std::array<double,1> &)> point_funct;

typedef std::function<int(const ubvector &,double,std::vector<double> &,
			  std::array<double,1> &)> fill_funct;

/// The MCMC object
mcmc_para_table<point_funct,fill_funct,std::array<double,1>,ubvector> mct;

class exc {

public:
  
  /** \brief A one-dimensional bimodal distribution

      Here, the variable 'log_weight' stores the natural logarithm of
      the objective function based on the parameters stored in \c pars.
      The object 'dat' stores any auxillary quantities which can be
      computed at every point in parameter space.
  */
  int bimodal(size_t nv, const ubvector &pars, double &log_weight,
	      std::array<double,1> &dat) {
    
    double x=pars[0];
    log_weight=log(exp(-x*x)*(sin(x-1.4)+1.0));
    dat[0]=x*x;
    return 0;
  }

  /** \brief Add auxillary quantities to 'line' so they can be
      stored in the table
  */
  int fill_line(const ubvector &pars, double log_weight, 
		std::vector<double> &line, std::array<double,1> &dat) {
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
		     std::array<double,1> &)>(&exc::bimodal),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);
  fill_funct fill_func=std::bind
    (std::mem_fn<int(const ubvector &,double,std::vector<double> &,
		     std::array<double,1> &)>(&exc::fill_line),&e,
     std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,
     std::placeholders::_4);

  vector<point_funct> bimodal_vec;
  bimodal_vec.push_back(bimodal_func);
  vector<fill_funct> fill_vec;
  fill_vec.push_back(fill_func);

  cout << "----------------------------------------------------------"
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
  mct.mcmc(1,low_bimodal,high_bimodal,bimodal_vec,fill_vec);

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
  std::vector<double> indep;
  size_t count=0;
  for(size_t j=0;j<t->get_nlines();j++) {
    for(size_t k=0;k<((size_t)(t->get("mult",j)+1.0e-8));k++) {
      if (count==ac_len) {
	indep.push_back(t->get("x2",j));
	count=0;
      }
      count++;
    }
  }

  // Use the independent samples to compute the final integral and
  // compare to the exact result
  double avg=vector_mean(indep);
  double std=vector_stddev(indep);
  tm.test_rel(avg,1.32513,10.0*std/sqrt(indep.size()),"ex_mcmc");
  cout << avg << " " << 1.32513 << " " << fabs(avg-1.32513) << " "
       << std << " " << 10.0*std/sqrt(indep.size()) << endl;
  
  tm.report();
  
  return 0;
}
// End of example
