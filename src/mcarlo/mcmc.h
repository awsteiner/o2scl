/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2016, Andrew W. Steiner
  
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
/** \file mcmc.h
    \brief File for definition of \ref o2scl::mcmc_base and 
    \ref o2scl::mcmc_table
*/
#ifndef O2SCL_MCMC_H
#define O2SCL_MCMC_H

#include <iostream>
#include <random>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/uniform_grid.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_file.h>
#include <o2scl/exception.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/cholesky.h>
#include <o2scl/vector.h>
#include <o2scl/multi_funct.h>

/** \brief Main namespace
    
    This file is documented in mcmc.h .
*/
namespace o2scl {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  
  /** \brief A generic MCMC simulation class

      This class generates a Markov chain of a user-specified
      function. The chain can be generated using the
      Metropolis-Hastings algorithm with a user-specified proposal
      distribution or using the affine-invariant sampling method of
      Goodman and Weare.

      By default, the Metropolis-Hastings algorithm is executed with a
      simple walk, with steps in each dimension of size \f$
      (\mathrm{high} - \mathrm{low})/\mathrm{step\_fac} \f$ with the
      denominator specified in \ref step_fac.

      The function type is a template type, \c func_t, which should
      be of the form 
      \code
      int f(size_t num_of_parameters, const vec_t &parameters,
            double &log_pdf, data_t &dat)
      \endcode
      which computes \c log_pdf, the natural logarithm of the function
      value, for any point in parameter space (any point between \c
      low and \c high ).

      If the function being simulated returns \ref mcmc_skip then the
      point is automatically rejected. After each acceptance or
      rejection, a user-specified "measurement" function (of type \c
      measure_t ) is called, which can be used to store the results.
      In order to stop the simulation, either this function or the
      probability distribution being simulated should return the value
      \ref mcmc_done .
      
      A generic proposal distribution can be specified in \ref
      set_proposal(). To go back to the default random walk method,
      one can call the function \ref unset_proposal().

      If \ref aff_inv is set to true and the number of walkers, \ref
      n_walk is set to a number larger than 1, then affine-invariant
      sampling is used. For affine-invariant sampling, the variable 
      \ref step_fac represents the value of \f$ a \f$, the 
      limits of the distribution for \f$ z \f$.

      In order to store data at each point, the user can store this
      data in any object of type \c data_t . If affine-invariant
      sampling is used, then each chain has it's own data object. The
      class keeps twice as many copies of these data object as would
      otherwise be required, in order to avoid copying of data objects
      in the case that the steps are accepted or rejected.

      \note This class is experimental.
  */
  template<class func_t, class measure_t,
    class data_t, class vec_t=ubvector> class mcmc_base {
    
  protected:
  
  /// Random number generator
  std::uniform_real_distribution<double> unif;
  
  /// Proposal distribution
  o2scl::prob_cond_mdim<vec_t> *prop_dist;
  
  /// If true, then use the user-specified proposal distribution
  bool pd_mode;

  /// If true, we are in the warm up phase
  bool warm_up;

  /// Current points in parameter space
  std::vector<vec_t> current;

  /// Data array
  std::vector<data_t> data_arr;

  /// Data switch array
  std::vector<bool> switch_arr;

  /// \name Interface customization
  //@{
  /** \brief Initializations before the MCMC 
   */
  virtual int mcmc_init() {
    return 0;
  }

  /** \brief Cleanup after the MCMC
   */
  virtual void mcmc_cleanup() {
    return;
  }

  /** \brief Function to run for the best point 
   */
  virtual void best_point(ubvector &best, double w_best, data_t &dat) {
    return;
  }
  //@}

  /// Index of the current walker
  size_t curr_walker;
  
  public:

  /// Integer to indicate completion
  static const int mcmc_done=-10;

  /// Integer to indicate rejection
  static const int mcmc_skip=-20;

  /// \name Output quantities
  //@{
  /// The number of Metropolis steps which were accepted
  size_t n_accept;

  /// The number of Metropolis steps which were rejected
  size_t n_reject;
  //@}

  /// \name Settings
  //@{
  /// If true, use affine-invariant Monte Carlo
  bool aff_inv;
  
  /// Stepsize factor (default 10.0)
  double step_fac;

  /** \brief Number of warm up steps (successful steps not
      iterations) (default 0)
	
      \note Not to be confused with <tt>warm_up</tt>, which is 
      a boolean local variable in some functions not an int.
  */
  size_t n_warm_up;

  /** \brief If non-zero, use as the seed for the random number 
      generator (default 0)
  */
  int user_seed;

  /// Output control (default 0)
  int verbose;

  /** \brief Maximum number of failed steps when generating initial points
      with affine-invariant sampling (default 1000)
  */
  size_t max_bad_steps;

  /** \brief Number of walkers for affine-invariant MC or 1 
      otherwise (default 1)
  */
  size_t n_walk;

  /** \brief If true, call the error handler if msolve() or
      msolve_de() does not converge (default true)
  */
  bool err_nonconv;
  //@}
  
  mcmc_base() : unif(0.0,1.0) {
    user_seed=0;
    n_warm_up=0;

    // MC step parameters
    aff_inv=false;
    pd_mode=false;
    step_fac=10.0;
    n_walk=1;
    err_nonconv=true;
    verbose=0;
    warm_up=false;
    max_bad_steps=1000;

    n_accept=0;
    n_reject=0;
    prop_dist=0;
  }

  /// \name Basic usage
  //@{
  /** \brief Perform an MCMC simulation

      Perform an MCMC simulation over \c nparams parameters starting
      at initial point \c init, limiting the parameters to be between
      \c low and \c high, using \c func as the objective function and
      calling the measurement function \c meas at each MC point.
  */
  virtual int mcmc(size_t nparams, vec_t &init,
		   vec_t &low, vec_t &high, func_t &func,
		   measure_t &meas) {

    // Fix the number of walkers if it is too small
    if (aff_inv && n_walk<=1) n_walk=2;
    
    // Fix 'step_fac' if it's less than or equal to zero
    if (step_fac<=0.0) {
      step_fac=10.0;
    }
    
    // Set RNG seed
    unsigned long int seed=time(0);
    if (user_seed!=0) {
      seed=user_seed;
    }
    std::mt19937 rd(seed);

    // Keep track of successful and failed MH moves
    n_accept=0;
    n_reject=0;

    // Warm-up flag, not to be confused with 'n_warm_up', i.e. the
    // number of warm_up iterations.
    warm_up=true;
    if (n_warm_up==0) warm_up=false;

    // Allocate current point and current weight
    current.resize(n_walk);
    std::vector<double> w_current(n_walk);
    w_current[0]=0.0;
    for(size_t i=0;i<n_walk;i++) {
      current[i].resize(nparams);
      w_current[i]=0.0;
    }

    // Initialize data and switch arrays
    data_arr.resize(2*n_walk);
    switch_arr.resize(n_walk);
    for(size_t i=0;i<switch_arr.size();i++) switch_arr[i]=false;
    
    // Next and best points in parameter space
    vec_t next(nparams), best(nparams);
    double w_next=0.0, w_best=0.0, q_prop=0.0;
    
    // Set current to initial point
    for(size_t i=0;i<nparams;i++) {
      current[0][i]=init[i];
    }

    // Return value of the measurement function
    int meas_ret=0;

    // Run init() function
    int iret=mcmc_init();
    if (iret!=0) {
      O2SCL_ERR("Function mcmc_init() failed in mcmc_base::mcmc().",
		o2scl::exc_einval);
      return iret;
    }

    if (verbose>=1) {
      if (aff_inv) {
	std::cout << "mcmc: Affine-invariant step, n_parameters="
	<< nparams << " n_walk=" << n_walk << std::endl;
      } else if (pd_mode==true) {
	std::cout << "mcmc: With proposal distribution, n_parameters="
	<< nparams << std::endl;
      } else {
	std::cout << "mcmc: Random-walk with uniform dist., n_parameters="
	<< nparams << std::endl;
      }
    }
    
    // ---------------------------------------------------
    // Compute initial point and initial weights
    
    if (aff_inv) {

      // Stretch-move steps

      // Initialize each walker in turn
      for(curr_walker=0;curr_walker<n_walk;curr_walker++) {

	size_t init_iters=0;
	bool done=false;

	while (!done) {

	  // Make a perturbation from the initial point
	  for(size_t ipar=0;ipar<nparams;ipar++) {
	    if (init[ipar]<low[ipar] || init[ipar]>high[ipar]) {
	      O2SCL_ERR((((std::string)"Parameter ")+std::to_string(ipar)+
			 " of "+std::to_string(nparams)+
			 " out of range (value="+std::to_string(init[ipar])+
			 ") in mcmc_base::mcmc().").c_str(),
			o2scl::exc_einval);
	    }
	    do {
	      current[curr_walker][ipar]=init[ipar]+(unif(rd)*2.0-1.0)*
		(high[ipar]-low[ipar])/100.0;
	    } while (current[curr_walker][ipar]>=high[ipar] ||
		     current[curr_walker][ipar]<=low[ipar]);
	  }
	
	  // Compute the weight
	  int iret=func(nparams,current[curr_walker],
			w_current[curr_walker],data_arr[curr_walker]);
	  
	  if (verbose>=1) {
	    std::cout.precision(4);
	    std::cout << "mcmc: ";
	    std::cout.width((int)(1.0+log10((double)(n_walk-1))));
	    std::cout << curr_walker << " " << w_current[curr_walker]
		      << " (initial; ai)" << std::endl;
	    std::cout.precision(6);
	  }

	  // ------------------------------------------------
	  
	  // Increment iteration count
	  init_iters++;

	  // If we have a good point, call the measurement function and
	  // stop the loop
	  if (iret==o2scl::success) {
	    meas_ret=meas(current[curr_walker],w_current[curr_walker],
			  curr_walker,true,data_arr[curr_walker]);
	    done=true;
	  } else if (init_iters>max_bad_steps) {
	    if (err_nonconv) {
	      O2SCL_ERR2("Initial walkers failed in ",
			 "mcmc_base::mcmc().",o2scl::exc_einval);
	    }
	    return 1;
	  }
	  
	  if (verbose>=2) {
	    std::cout << "Press a key and type enter to continue. ";
	    char ch;
	    std::cin >> ch;
	  }
	  
	}

      }

    } else {

      curr_walker=0;
      
      // Uniform random-walk steps

      // Compute weight for initial point
      int iret=func(nparams,current[0],w_current[0],data_arr[0]);

      if (verbose>=1) {
	std::cout.precision(4);
	std::cout << "mcmc: " << w_current[0] << " (initial)" << std::endl;
	std::cout.precision(6);
      }

      if (iret!=o2scl::success) {
	if (err_nonconv) {
	  O2SCL_ERR("Initial weight vanished in mcmc_base::mcmc().",
		    o2scl::exc_einval);
	}
	return 2;
      }

      meas_ret=meas(current[0],w_current[0],0,true,data_arr[0]);

      best=current[0];
      w_best=w_current[0];
      best_point(best,w_best,data_arr[0]);

      if (meas_ret==mcmc_done) {
	mcmc_cleanup();
	return 0;
      }

    }
    
    // ---------------------------------------------------

    bool main_done=false;
    size_t mcmc_iters=0;
    
    while (!main_done) {

      // Walker to move (or zero when aff_inv is false)
      curr_walker=0;
      double smove_z=0.0;
    
      // ---------------------------------------------------
      // Select next point
    
      if (aff_inv) {

	// Choose walker to move
	curr_walker=mcmc_iters % n_walk;
	
	// Choose jth walker
	size_t ij;
	do {
	  ij=((size_t)(unif(rd)*((double)n_walk)));
	} while (ij==curr_walker || ij>=n_walk);
	
	// Select z 
	double p=unif(rd);
	double a=step_fac;
	smove_z=(1.0-2.0*p+2.0*a*p+p*p-2.0*a*p*p+a*a*p*p)/a;
	
	// Create new trial point
	for(size_t i=0;i<nparams;i++) {
	  next[i]=current[ij][i]+smove_z*(current[curr_walker][i]-
					  current[ij][i]);
	}
	
      } else if (pd_mode) {
	
	// Use proposal distribution and compute associated weight
	(*prop_dist)(current[0],next);
	q_prop=prop_dist->log_pdf(current[0],next)-
	  prop_dist->log_pdf(next,current[0]);
	if (!std::isfinite(q_prop)) {
	  O2SCL_ERR2("Proposal distribution not finite in ",
		     "mcmc_base::mcmc().",o2scl::exc_efailed);
	}

      } else {

	// Uniform random-walk step
	for(size_t k=0;k<nparams;k++) {
	  next[k]=current[0][k]+(unif(rd)*2.0-1.0)*
	    (high[k]-low[k])/step_fac;
	}
      
      }

      // ---------------------------------------------------
      // Compute next weight
      
      int iret=o2scl::success;
      // If the next point out of bounds, ensure that the point is rejected
      for(size_t k=0;k<nparams;k++) {
	if (next[k]<=low[k] || next[k]>=high[k]) iret=mcmc_skip;
      }
      if (iret!=mcmc_skip) {
	if (switch_arr[curr_walker]==false) {
	  iret=func(nparams,next,w_next,data_arr[curr_walker+n_walk]);
	} else {
	  iret=func(nparams,next,w_next,data_arr[curr_walker]);
	}
      }
      if (iret==o2scl::success && w_best>w_next) {
	best=next;
	w_best=w_next;
	if (switch_arr[curr_walker]==false) {
	  best_point(best,w_best,data_arr[curr_walker+n_walk]);
	} else {
	  best_point(best,w_best,data_arr[curr_walker]);
	}
      }
      
      // ---------------------------------------------------
    
      bool accept=false;

      if (iret==o2scl::success) {
	double r=unif(rd);
	
	if (aff_inv) {
	  if (r<pow(smove_z,((double)nparams)-1.0)*
	      exp(w_next-w_current[curr_walker])) {
	    accept=true;
	  }
	  if (verbose>=1) {
	    std::cout.precision(4);
	    std::cout << "mcmc: ";
	    std::cout.width((int)(1.0+log10((double)(nparams-1))));
	    std::cout << curr_walker << " "
		      << w_current[curr_walker] << " " << w_next << " "
		      << smove_z << " ratio: "
		      << pow(smove_z,((double)nparams)-1.0)*
	      exp(w_next-w_current[curr_walker])
		      << " accept: " << accept << std::endl;
	    std::cout.precision(6);
	  }
	} else if (pd_mode) {
	  if (r<exp(w_next-w_current[0]+q_prop)) {
	    accept=true;
	  }
	  if (verbose>=1) {
	    std::cout.precision(4);
	    std::cout << "mcmc: " << w_current[0] 
		      << " " << w_next << " " << q_prop << " ratio: "
		      << exp(w_next-w_current[0]+q_prop)
		      << " accept: " << accept << std::endl;
	    std::cout.precision(6);
	  }
	} else {
	  // Metropolis algorithm
	  if (r<exp(w_next-w_current[0])) {
	    accept=true;
	  }
	  if (verbose>=1) {
	    std::cout.precision(4);
	    std::cout << "mcmc: " << w_current[0] << " " << w_next
		      << " ratio: "
		      << exp(w_next-w_current[0])
		      << " accept: " << accept << std::endl;
	    std::cout.precision(6);
	  }
	}

	// End of 'if (iret==o2scl::success)'
      }

      if (accept) {
	  
	n_accept++;
	  
	// Store results from new point
	if (!warm_up) {
	  if (switch_arr[curr_walker]==false) {
	    meas_ret=meas(next,w_next,curr_walker,true,
			  data_arr[curr_walker+n_walk]);
	  } else {
	    meas_ret=meas(next,w_next,curr_walker,true,
			  data_arr[curr_walker]);
	  }
	}

	// Prepare for next point
	current[curr_walker]=next;
	w_current[curr_walker]=w_next;
	switch_arr[curr_walker]=!(switch_arr[curr_walker]);
	  
      } else {
	    
	// Point was rejected
	n_reject++;

	// Repeat measurement of old point
	if (!warm_up) {
	  if (switch_arr[curr_walker]==false) {
	    meas_ret=meas(current[curr_walker],w_current[curr_walker],
			  curr_walker,false,data_arr[curr_walker]);
	  } else {
	    meas_ret=meas(current[curr_walker],w_current[curr_walker],
			  curr_walker,false,data_arr[curr_walker+n_walk]);
	  }
	}

      }

      if (iret==mcmc_done) {
	main_done=true;
      } else {
	if (meas_ret!=0) {
	  main_done=true;
	  if (meas_ret!=mcmc_done && err_nonconv) {
	    O2SCL_ERR((((std::string)"Measurement function returned ")+
		       std::to_string(meas_ret)+
		       " in mcmc_base::mcmc().").c_str(),
		      o2scl::exc_efailed);
	  }
	}
	
	mcmc_iters++;
	
	if (warm_up && mcmc_iters==n_warm_up) {
	  warm_up=false;
	  mcmc_iters=0;
	  n_accept=0;
	  n_reject=0;
	  if (verbose>=1) {
	    std::cout << "Finished warmup." << std::endl;
	  }
	}
      }

      if (verbose>=2) {
	std::cout << "Press a key and type enter to continue. ";
	char ch;
	std::cin >> ch;
      }

      // --------------------------------------------------------------
      // End of main loop
    }
    
    // --------------------------------------------------------------
    
    mcmc_cleanup();
    
    return 0;
  }
  //@}

  /// \name Proposal distribution
  //@{
  /** \brief Set the proposal distribution
   */
  virtual void set_proposal(o2scl::prob_cond_mdim<vec_t> &p) {
    prop_dist=&p;
    pd_mode=true;
    aff_inv=false;
    n_walk=1;
    return;
  }

  /** \brief Go back to random-walk Metropolis with a uniform distribution
   */
  virtual void unset_proposal() {
    if (pd_mode) {
      prop_dist=0;
      pd_mode=false;
    }
    aff_inv=false;
    n_walk=1;
    return;
  }
  //@}

  };

  /** \brief A generic MCMC simulation class writing data to a 
      \ref o2scl::table_units object

      This class performs a MCMC simulation and stores the 
      results in a \ref o2scl::table_units object. The
      user must specify the column names and units in 
      \ref set_names_units before \ref mcmc() is called.

      The function \ref add_line is the measurement function of type
      \c measure_t in the parent. The overloaded function \ref mcmc()
      in this class works a bit differently in that it takes a
      function object (type \c fill_t) of the form
      \code
      int fill_func(const vec_t &pars, double log_weight, 
      std::vector<double> &line, data_t &dat);
      \endcode
      which should store any auxillary values stored in the data
      object to \c line, in order to be added to the table.

      The output table will contain the parameters, the logarithm of
      the function (called "log_wgt") and a multiplying factor called
      "mult". This "fill" function is called only when a step is
      accepted and the multiplier for that row is set to 1. If a
      future step is rejected, then the multiplier is increased by
      one, rather than adding the same row to the table again.

      This class forms the basis of the MCMC used in the Bayesian
      analysis of neutron star mass and radius in
      http://github.com/awsteiner/bamr .

      \note This class is experimental.
  */
  template<class func_t, class fill_t, class data_t, class vec_t=ubvector>
    class mcmc_table : public mcmc_base<func_t,
    std::function<int(const vec_t &,double,size_t,bool,data_t &)>,
    data_t,vec_t> {

  protected:

  /// Measurement functor type for the parent
  typedef std::function<int(const vec_t &,double,size_t,bool,data_t &)>
  internal_measure_t;
  
  /// Type of parent class
  typedef mcmc_base<func_t,internal_measure_t,data_t,vec_t> parent_t;

  /// Parameter names
  std::vector<std::string> col_names;
    
  /// Parameter units
  std::vector<std::string> col_units;
    
  /// Main data table for Markov chain
  std::shared_ptr<o2scl::table_units<> > tab;

  /** \brief MCMC initialization function

      This function sets the column names and units.
  */
  virtual int mcmc_init() {

    if (this->verbose>=2) {
      std::cout << "Start mcmc_table::mcmc_init()." << std::endl;
    }

    // -----------------------------------------------------------
    // Init table

    std::string s, u;
    tab->clear_table();
    tab->new_column("mult");
    tab->new_column("log_wgt");
    for(size_t i=0;i<col_names.size();i++) {
      tab->new_column(col_names[i]);
      if (col_units[i].length()>0) {
	tab->set_unit(col_names[i],col_units[i]);
      }
    }

    walker_rows.resize(this->n_walk);
    for(size_t i=0;i<this->n_walk;i++) {
      walker_rows[i]=-1;
    }
    
    if (this->verbose>=2) {
      std::cout << "mcmc: Table column names and units: " << std::endl;
      for(size_t i=0;i<tab->get_ncolumns();i++) {
	std::cout << tab->get_column_name(i) << " "
		  << tab->get_unit(tab->get_column_name(i)) << std::endl;
      }
    }
    
    if (this->verbose>=2) {
      std::cout << "End mcmc_table::mcmc_init()." << std::endl;
    }

    return 0;
  }

  /** \brief Fill \c line with data for insertion into the table
   */
  virtual int fill_line(const vec_t &pars, double log_weight, 
			 std::vector<double> &line, data_t &dat,
			 fill_t &fill) {

    // Initial multiplier
    line.push_back(1.0);
    line.push_back(log_weight);
    for(size_t i=0;i<pars.size();i++) {
      line.push_back(pars[i]);
    }
    return fill(pars,log_weight,line,dat);
  }

  /** \brief Record the last row in the table which corresponds
      to each walker
   */
  std::vector<int> walker_rows;
  
  public:

  mcmc_table() : tab(new o2scl::table_units<>) {
  }

  /// \name Basic usage
  //@{
  /** \brief Set the table names and units
   */
  virtual void set_names_units(std::vector<std::string> names,
			       std::vector<std::string> units) {
    col_names=names;
    col_units=units;
    return;
  }
  
  /** \brief Perform an MCMC simulation

      Perform an MCMC simulation over \c nparams parameters starting
      at initial point \c init, limiting the parameters to be between
      \c low and \c high, using \c func as the objective function and
      calling the measurement function \c meas at each MC point.
  */
  virtual int mcmc(size_t nparams, vec_t &init,
		   vec_t &low, vec_t &high, func_t &func,
		   fill_t &fill) {

    internal_measure_t meas=std::bind
    (std::mem_fn<int(const vec_t &,double,size_t,bool,data_t &,fill_t &)>
     (&mcmc_table::add_line),this,std::placeholders::_1,
     std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
     std::placeholders::_5,std::ref(fill));
    
    return parent_t::mcmc(nparams,init,low,high,func,meas);
  }

  /** \brief Get the output table
   */
  std::shared_ptr<o2scl::table_units<> > get_table() {
    return tab;
  }
  
  /** \brief Set the output table
   */
  void set_table(std::shared_ptr<o2scl::table_units<> > &t) {
    tab=t;
    return;
  }
  
  /** \brief A measurement function which adds the point to the
      table
  */
  virtual int add_line(const vec_t &pars, double log_weight,
		       size_t ix, bool new_meas, data_t &dat,
		       fill_t &fill) {

    // Test to see if we need to add a new line of data or increment
    // the weight on the previous line
    if (new_meas==true) {

      std::vector<double> line;
      int fret=fill_line(pars,log_weight,line,dat,fill);
      
      if (fret!=o2scl::success) {
	// If we're done, we stop before adding the last point to the
	// table. This is important because otherwise the last line in
	// the table will always only have unit multiplicity, which
	// may or may not be correct.
	if (fret==this->mcmc_done) {
	  if (this->verbose>=1) {
	    std::cout << "Fill function returned mcmc_done. " 
		      << "Stopping run." << std::endl;
	  }
	}
	else {
	  if (this->verbose>=1) {
	    std::cout << "Fill function returned " << fret
	    << ". Stopping run." << std::endl;
	  }
	}
	return this->mcmc_done;
      }
      
      if (line.size()!=tab->get_ncolumns()) {
	std::cout << "line: " << line.size() << " columns: "
	<< tab->get_ncolumns() << std::endl;
	O2SCL_ERR("Table misalignment in mcmc_table::add_line().",
		  exc_einval);
      }
      
      if (this->verbose>=2) {
	std::cout << "mcmc: Adding line:" << std::endl;
	std::vector<std::string> sc_in, sc_out;
	for(size_t k=0;k<line.size();k++) {
	  sc_in.push_back(tab->get_column_name(k)+": "+
			  o2scl::dtos(line[k]));
	}
	o2scl::screenify(line.size(),sc_in,sc_out);
	for(size_t k=0;k<sc_out.size();k++) {
	  std::cout << sc_out[k] << std::endl;
	}
      }

      walker_rows[this->curr_walker]=tab->get_nlines();
      tab->line_of_data(line.size(),line);

    } else if (tab->get_nlines()>0) {
	
      // Otherwise, just increment the multiplier on the previous line
      if (walker_rows[this->curr_walker]<0 ||
	  walker_rows[this->curr_walker]>=((int)(tab->get_nlines()))) {
	std::cout << "nlines: " << tab->get_nlines() << std::endl;
	std::cout << "walker: " << this->curr_walker << std::endl;
	std::cout << "row: " << walker_rows[this->curr_walker]
	<< std::endl;
	O2SCL_ERR2("Sanity in row counting in ",
		   "mcmc_table::add_line().",o2scl::exc_esanity);
      }

      double mult_old=tab->get("mult",walker_rows[this->curr_walker]);
      tab->set("mult",walker_rows[this->curr_walker],mult_old+1.0);
      
      if (this->verbose>=2) {
	std::cout << "mcmc: Updating line:" << std::endl;
	std::vector<std::string> sc_in, sc_out;
	for(size_t k=0;k<tab->get_ncolumns();k++) {
	  sc_in.push_back
	    (tab->get_column_name(k)+": "+
	     o2scl::dtos(tab->get(tab->get_column_name(k),
				  walker_rows[this->curr_walker])));
	}
	o2scl::screenify(tab->get_ncolumns(),sc_in,sc_out);
	for(size_t k=0;k<sc_out.size();k++) {
	  std::cout << sc_out[k] << std::endl;
	}
      }
      
    }
      
    return 0;
  }
  //@}

  /** \brief Reaverage the data into blocks of a fixed
      size in order to avoid autocorrelations
      
      \note The number of blocks \c n_blocks must be larger than the
      current table size. This function expects to find a column named
      "mult" which contains the multiplicity of each column, as is the
      case after a call to \ref mcmc_base::mcmc().
      
      This function is useful to remove autocorrelations to the table
      so long as the autocorrelation length is shorter than the block
      size. This function does not compute the autocorrelation length
      to check that this is the case.
  */
  void reblock(size_t n_blocks) {
    size_t n=tab->get_nlines();
    if (n_blocks>n) {
      O2SCL_ERR2("Cannot reblock. Not enough data in ",
		"mcmc_table::reblock().",o2scl::exc_einval);
    }
    size_t n_block=n/n_blocks;
    size_t m=tab->get_ncolumns();
    for(size_t j=0;j<n_blocks;j++) {
      double mult=0.0;
      ubvector dat(m);
      for(size_t i=0;i<m;i++) {
	dat[i]=0.0;
      }
      for(size_t k=j*n_block;k<(j+1)*n_block;k++) {
	mult+=(*tab)["mult"][k];
	for(size_t i=1;i<m;i++) {
	  dat[i]+=(*tab)[i][k]*(*tab)["mult"][k];
	}
      }
      tab->set("mult",j,mult);
      for(size_t i=1;i<m;i++) {
	dat[i]/=mult;
	tab->set(i,j,dat[i]);
      }
    }
    tab->set_nlines(n_blocks);
    return;
  }
  
  };
  
  // End of namespace
}

#endif
