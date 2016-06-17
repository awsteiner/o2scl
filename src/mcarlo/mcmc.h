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

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/rng_gsl.h>
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

      \note This class is experimental.
      
      \todo Fix Hastings steps and add better testing
   */
  template<class func_t, class measure_t,
    class data_t, class vec_t=ubvector> class mcmc_base {
    
  protected:
  
  /// Random number generator
  o2scl::rng_gsl gr;
  
  //// A Hastings distribution
  o2scl::prob_dens_mdim<vec_t> hast;
  
  /// If true, then use Metropolis-Hastings with a multivariate Gaussian
  bool hg_mode;

  /// If true, we are in the warm up phase
  bool warm_up;

  /// Current points in parameter space
  std::vector<vec_t> current;

  /// Data array
  std::vector<data_t> data_arr;

  /// Data switch array
  std::vector<bool> switch_arr;

  public:

  /// Integer to indicate completion
  static const int mcmc_done=-10;
  
  /// The number of Metropolis steps which were accepted
  size_t n_accept;

  /// The number of Metropolis steps which were rejected
  size_t n_reject;

  /// \name settings
  //@{
  /// If true, use affine-invariant Monte Carlo
  bool aff_inv;
  
  /// MCMC stepsize factor (default 10.0)
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

  /// (default 1000)
  size_t max_bad_steps;

  /** \brief Number of walkers for affine-invariant MC or 1 
      otherwise (default 1)
  */
  size_t nwalk;
  //@}

  /** \brief If true, call the error handler if msolve() or
      msolve_de() does not converge (default true)
  */
  bool err_nonconv;
  
  mcmc_base() {
    user_seed=0;
    n_warm_up=0;

    // MC step parameters
    aff_inv=false;
    hg_mode=false;
    step_fac=10.0;
    nwalk=1;
    err_nonconv=true;
    verbose=0;
    warm_up=false;
    max_bad_steps=1000;

    n_accept=0;
    n_reject=0;
  }
  
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

  /** \brief Desc
   */
  virtual void best_point(ubvector &best, double w_best, data_t &dat) {
    return;
  }
  
  /** \brief Set the distribution for Metropolis-Hastings
   */
  virtual void set_hastings(o2scl::prob_dens_mdim<vec_t> &p) {
    hast=p;
    hg_mode=true;
    aff_inv=false;
    return;
  }
  
  /** \brief The main MCMC function
   */
  virtual int mcmc(size_t nparams, vec_t &init,
		   vec_t &low, vec_t &high, func_t &func,
		   measure_t &meas) {

    // Fix the number of walkers if it is too small
    if (aff_inv && nwalk<=1) nwalk=2;
    
    // Fix 'step_fac' if it's less than or equal to zero
    if (step_fac<=0.0) {
      step_fac=10.0;
    }
    
    // Set RNG seed
    unsigned long int seed=time(0);
    if (user_seed!=0) {
      seed=user_seed;
    }
    gr.set_seed(seed);

    // Keep track of successful and failed MH moves
    n_accept=0;
    n_reject=0;

    // Warm-up flag, not to be confused with 'n_warm_up', i.e. the
    // number of warm_up iterations.
    warm_up=true;
    if (n_warm_up==0) warm_up=false;

    // Allocate current point and current weight
    current.resize(nwalk);
    std::vector<double> w_current(nwalk);
    w_current[0]=0.0;
    for(size_t i=0;i<nwalk;i++) {
      current[i].resize(nparams);
      w_current[i]=0.0;
    }
    double q_current=0.0;

    // Initialize data and switch arrays
    data_arr.resize(2*nwalk);
    switch_arr.resize(nwalk);
    for(size_t i=0;i<switch_arr.size();i++) switch_arr[i]=false;
    
    // Next and best points in parameter space
    vec_t next(nparams), best(nparams);
    double w_next=0.0, w_best=0.0, q_next=0.0;
    
    // Set current to initial point
    for(size_t i=0;i<nparams;i++) {
      current[0][i]=init[i];
    }

    // Return value of the measurement function
    int meas_ret=0;

    // Run init() function
    int iret=mcmc_init();
    if (iret!=0) {
      O2SCL_ERR("Function mcmc_init failed in mcmc_base::mcmc().",
		o2scl::exc_einval);
      return iret;
    }

    if (verbose>=1) {
      if (aff_inv) {
	std::cout << "mcmc: Affine-invariant step, n_parameters="
	<< nparams << " nwalk=" << nwalk << std::endl;
      } else if (hg_mode==true) {
	std::cout << "mcmc: Metropolis-Hastings, n_parameters="
	<< nparams << std::endl;
      } else {
	std::cout << "mcmc: Simple metropolis, n_parameters="
	<< nparams << std::endl;
      }
    }
    
    // ---------------------------------------------------
    // Compute initial point and initial weights
    
    if (aff_inv) {

      // Stretch-move steps
      size_t ij_best=0;

      // Initialize each walker in turn
      for(size_t ij=0;ij<nwalk;ij++) {

	size_t init_iters=0;
	bool done=false;

	while (!done) {

	  // Make a perturbation from the initial point
	  for(size_t ik=0;ik<nparams;ik++) {
	    if (init[ik]<low[ik] || init[ik]>high[ik]) {
	      O2SCL_ERR((((std::string)"Parameter ")+std::to_string(ik)+
			 " of "+std::to_string(nparams)+
			 " out of range (value="+std::to_string(init[ik])+
			 ") in mcmc_base::mcmc().").c_str(),
			o2scl::exc_einval);
	    }
	    do {
	      current[ij][ik]=init[ik]+(gr.random()*2.0-1.0)*
		(high[ik]-low[ik])/100.0;
	    } while (current[ij][ik]>=high[ik] || current[ij][ik]<=low[ik]);
	  }
	
	  // Compute the weight
	  w_current[ij]=func(nparams,current[ij],data_arr[ij]);
	  
	  if (verbose>=1) {
	    std::cout.precision(4);
	    std::cout << "mcmc: ";
	    std::cout.width((int)(1.0+log10((double)(nwalk-1))));
	    std::cout << ij << " " << w_current[ij]
		      << " (initial; ai)" << std::endl;
	    std::cout.precision(6);
	  }

	  // ------------------------------------------------
	  
	  // Increment iteration count
	  init_iters++;

	  // If we have a good point, call the measurement function and
	  // stop the loop
	  if (w_current[ij]>0.0) {
	    meas_ret=meas(current[ij],w_current[ij],ij,true,data_arr[ij]);
	    done=true;
	  } else if (init_iters>max_bad_steps) {
	    if (err_nonconv) {
	      O2SCL_ERR("Initial walkers failed.",o2scl::exc_einval);
	    }
	    return 1;
	  }
	}

      }

    } else {

      // Normal or Metropolis-Hastings steps

      // Compute weight for initial point
      w_current[0]=func(nparams,current[0],data_arr[0]);

      if (verbose>=1) {
	std::cout.precision(4);
	std::cout << "mcmc: " << w_current[0] << " (initial)" << std::endl;
	std::cout.precision(6);
      }

      if (w_current[0]<=0.0) {
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

      // Compute the initial Hastings proposal weight
      if (hg_mode>0) {
	q_current=hast.pdf(current[0]);
      }
    
    }
    
    // ---------------------------------------------------

    bool main_done=false;
    size_t mcmc_iters=0;
    
    while (!main_done) {

      // Walker to move (or zero when aff_inv is false)
      size_t ik=0;
      double smove_z=0.0;
    
      // ---------------------------------------------------
      // Select next point
    
      if (aff_inv) {

	// Choose walker to move
	ik=mcmc_iters % nwalk;
      
	bool in_bounds;
	size_t step_iters=0;
      
	do {

	  in_bounds=true;
	
	  // Choose jth walker
	  size_t ij;
	  do {
	    ij=((size_t)(gr.random()*((double)nwalk)));
	  } while (ij==ik || ij>=nwalk);
	
	  // Select z 
	  double p=gr.random();
	  double a=step_fac;
	  smove_z=(1.0-2.0*p+2.0*a*p+p*p-2.0*a*p*p+a*a*p*p)/a;
	
	  // Create new trial point
	  for(size_t i=0;i<nparams;i++) {
	    next[i]=current[ij][i]+smove_z*(current[ik][i]-current[ij][i]);
	    if (next[i]>=high[i] || next[i]<=low[i]) {
	      in_bounds=false;
	    }
	  }
	
	  step_iters++;
	  if (step_iters==1000) {
	    if (err_nonconv) {
	      O2SCL_ERR("Failed to find suitable step in mcmc_base::mcmc().",
			o2scl::exc_einval);
	    }
	    return 2;
	  }

	} while (in_bounds==false);	

      } else if (hg_mode>0) {
      
	// Make a Metropolis-Hastings step based on previous data
	hast(next);
	q_next=hast.pdf(next);

      } else {

	// Make a step, ensure that we're in bounds and that
	// the masses are not too large
	for(size_t k=0;k<nparams;k++) {
	  
	  next[k]=current[0][k]+(gr.random()*2.0-1.0)*
	    (high[k]-low[k])/step_fac;
	
	  // If it's out of range, redo step near boundary
	  if (next[k]<low[k]) {
	    next[k]=low[k]+gr.random()*(high[k]-low[k])/step_fac;
	  } else if (next[k]>high[k]) {
	    next[k]=high[k]-gr.random()*(high[k]-low[k])/step_fac;
	  }
	  
	  if (next[k]<low[k] || next[k]>high[k]) {
	    O2SCL_ERR("Sanity check in parameter step in mcmc_base::mcmc().",
		      o2scl::exc_esanity);
	  }
	}
      
      }

      // ---------------------------------------------------
      // Compute next weight

      if (switch_arr[ik]==false) {
	w_next=func(nparams,next,data_arr[ik+nwalk]);
      } else {
	w_next=func(nparams,next,data_arr[ik]);
      }
      if (w_next>w_best) {
	best=next;
	w_best=w_next;
	if (switch_arr[ik]==false) {
	  best_point(best,w_best,data_arr[ik+nwalk]);
	} else {
	  best_point(best,w_best,data_arr[ik]);
	}
      }
      
      // ---------------------------------------------------
    
      bool accept=false;
      double r=gr.random();

      // Metropolis algorithm
      if (aff_inv) {
	if (r<pow(smove_z,((double)nwalk)-1.0)*w_next/w_current[ik]) {
	  accept=true;
	}
	if (verbose>=1) {
	  std::cout.precision(4);
	  std::cout << "mcmc: ";
	  std::cout.width((int)(1.0+log10((double)(nwalk-1))));
	  std::cout << ik << " " << w_current[ik] << " " << w_next << " "
		    << pow(smove_z,((double)nwalk)-1.0) << " ratio: "
		    << pow(smove_z,((double)nwalk)-1.0)*w_next/w_current[ik]
		    << " accept: " << accept << std::endl;
	  std::cout.precision(6);
	}
      } else if (hg_mode>0) {
	if (r<w_next*q_current/w_current[0]/q_next) {
	  accept=true;
	}
	if (verbose>=1) {
	  std::cout.precision(4);
	  std::cout << "mcmc: " << w_current[0] << " " << q_current
		    << " " << w_next << " " << q_next << " ratio: "
		    << w_next*q_current/w_current[0]/q_next
		    << " accept: " << accept << std::endl;
	  std::cout.precision(6);
	}
      } else {
	if (r<w_next/w_current[0]) {
	  accept=true;
	}
	if (verbose>=1) {
	  std::cout.precision(4);
	  std::cout << "mcmc: " << w_current[0] << " " << w_next
		    << " ratio: " << w_next/w_current[0]
		    << " accept: " << accept << std::endl;
	  std::cout.precision(6);
	}
      }

      if (accept) {
	  
	n_accept++;
	  
	// Store results from new point
	if (!warm_up) {
	  if (switch_arr[ik]==false) {
	    meas_ret=meas(next,w_next,ik,true,data_arr[ik+nwalk]);
	  } else {
	    meas_ret=meas(next,w_next,ik,true,data_arr[ik]);
	  }
	}

	// Prepare for next point
	current[ik]=next;
	w_current[ik]=w_next;
	switch_arr[ik]=!(switch_arr[ik]);
	  
      } else {
	    
	// Point was rejected
	n_reject++;

	// Repeat measurement of old point
	if (!warm_up) {
	  if (switch_arr[ik]==false) {
	    meas_ret=meas(current[ik],w_current[ik],ik,false,
			  data_arr[ik]);
	  } else {
	    meas_ret=meas(current[ik],w_current[ik],ik,false,
			  data_arr[ik+nwalk]);
	  }
	}

      }

      if (meas_ret!=0) {
	main_done=true;
	if (meas_ret!=mcmc_done && err_nonconv) {
	  O2SCL_ERR((((std::string)"Measurement function returned ")+
		     std::to_string(meas_ret)+" in mcmc_base::mcmc().").c_str(),
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
      
      // --------------------------------------------------------------
      // End of main loop
    }
    
    // --------------------------------------------------------------
    
    mcmc_cleanup();
    
    return 0;
  }

  };

  /** \brief A generic MCMC simulation class
      
      \note This class is experimental.
   */
  template<class func_t, class measure_t, class data_t, class vec_t=ubvector>
    class mcmc_table : public mcmc_base<func_t,measure_t,data_t,vec_t> {
    
  protected:
    
  /// Parameter names
  std::vector<std::string> param_names;
    
  /// Parameter units
  std::vector<std::string> param_units;
  //@}
    
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
    tab->new_column("mult");
    tab->new_column("weight");
    for(size_t i=0;i<param_names.size();i++) {
      tab->new_column(((std::string)"param_")+param_names[i]);
      if (param_units[i].length()>0) {
	tab->set_unit(((std::string)"param_")+param_names[i],
		      param_units[i]);
      }
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
  virtual void fill_line(const vec_t &pars, double weight, 
			 std::vector<double> &line, data_t &dat) {

    // Initial multiplier
    line.push_back(1.0);
    line.push_back(weight);
    for(size_t i=0;i<pars.size();i++) {
      line.push_back(pars[i]);
    }
      
    return;
  }
    
  public:

  mcmc_table() : tab(new o2scl::table_units<>) {
  }
    
  /** \brief A measurement function which adds the point to the
      table
   */
  virtual int add_line(const vec_t &pars, double weight,
		       size_t ix, bool new_meas, data_t &dat) {

    // Test to see if we need to add a new line of data or
    // increment the weight on the previous line
    if (tab->get_nlines()<=(this->nwalk-1) || new_meas==true) {
	
      std::vector<double> line;
      fill_line(pars,weight,line,dat);
      
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
	std::cout << "Press a key and enter to continue." << std::endl;
	char ch;
	std::cin >> ch;
      }
      
      tab->line_of_data(line.size(),line);
	
    } else if (tab->get_nlines()>0) {
	
      // Otherwise, just increment the multiplier on the previous line
      tab->set("mult",tab->get_nlines()-this->nwalk,
	       tab->get("mult",tab->get_nlines()-this->nwalk)+1.0);

      if (this->verbose>=2) {
	std::cout << "mcmc: Updating line:" << std::endl;
	std::vector<std::string> sc_in, sc_out;
	for(size_t k=0;k<tab->get_ncolumns();k++) {
	  sc_in.push_back
	    (tab->get_column_name(k)+": "+
	     o2scl::dtos(tab->get(tab->get_column_name(k),
				  tab->get_nlines()-this->nwalk)));
	}
	o2scl::screenify(tab->get_ncolumns(),sc_in,sc_out);
	for(size_t k=0;k<sc_out.size();k++) {
	  std::cout << sc_out[k] << std::endl;
	}
	std::cout << "Press a key and enter to continue." << std::endl;
	char ch;
	std::cin >> ch;
      }
      
    }
      
    return 0;
  }
    
  /** \brief Set the table names and units
   */
  virtual void set_names_units(std::vector<std::string> names,
			       std::vector<std::string> units) {
    param_names=names;
    param_units=units;
    return;
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
    return tab;
  }
  
  };
  
  // End of namespace
}

#endif
