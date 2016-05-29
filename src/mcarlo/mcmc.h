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
  
  typedef std::function<int(const ubvector &,double,
			    size_t,bool)> measure_funct;
  
  /** \brief A generic MCMC simulation class

      \note This class is experimental.
   */
  template<class func_t=o2scl::multi_funct11,
    class measure_t=measure_funct, class vec_t=ubvector> class mcmc_base {
    
  protected:
  
  /// Random number generator
  o2scl::rng_gsl gr;
  
  /// A Gaussian probability distribution
  o2scl::prob_dens_gaussian pdg;

  //// Desc
  o2scl::prob_dens_mdim<vec_t> hast;
  
  /// If true, then use Metropolis-Hastings with a multivariate Gaussian
  int hg_mode;

  /// If true, we are in the warm up phase
  bool warm_up;

  /// Current points in parameter space
  std::vector<vec_t> current;

  public:

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
  int n_warm_up;

  /** \brief If non-zero, use as the seed for the random number 
      generator (default 0)
  */
  int user_seed;

  /// Output control (default 0)
  int verbose;

  /// (default 1000)
  int max_bad_steps;

  /** \brief Number of walkers for affine-invariant MC or 1 
      otherwise (default 1)
  */
  int nwalk;

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
    hg_mode=0;
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
  virtual void best_point(ubvector &best, double w_best) {
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

    if (step_fac<1.0) {
      step_fac=1.0;
    }
    
    // Set RNG seed
    unsigned long int seed=time(0);
    if (user_seed!=0) {
      seed=user_seed;
    }
    gr.set_seed(seed);
    pdg.set_seed(seed);

    // Keep track of successful and failed MH moves
    n_accept=0;
    n_reject=0;

    // Warm-up flag, not to be confused with 'n_warm_up', i.e. the
    // number of warm_up iterations.
    warm_up=true;
    if (n_warm_up==0) warm_up=false;

    current.resize(1);
    std::vector<double> w_current(1);
    w_current[0]=0.0;
    
    // For stretch-moves, allocate for each walker
    if (aff_inv) {
      current.resize(nwalk);
      w_current.resize(nwalk);
      for(size_t i=0;i<nwalk;i++) {
	w_current[i]=0.0;
      }
    }
    
    // Allocate memory for in all points
    for(size_t i=0;i<nwalk;i++) {
      current[i].resize(nparams);
    }

    // Entry objects (Must be after read_input() since nsources is set
    // in that function.)
    vec_t next(nparams), best(nparams);

    // Weights for each entry
    double w_next=0.0, w_best=0.0;

    for(size_t i=0;i<nparams;i++) {
      current[0][i]=init[i];
    }
    
    // Run init() function
    int iret=mcmc_init();
    if (iret!=0) {
      O2SCL_ERR("init function failed.",o2scl::exc_einval);
      return iret;
    }
      
    // ---------------------------------------------------
    // Compute initial point and initial weights
    
    int meas_ret=0;
      
    double q_current=0.0, q_next=0.0;

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
	    do {
	      current[ij][ik]=init[ik]+(gr.random()*2.0-1.0)*
		(high[ik]-low[ik])/100.0;
	    } while (current[ij][ik]>=high[ik] || current[ij][ik]<=low[ik]);
	  }
	
	  // Compute the weight
	  w_current[ij]=func(nparams,current[ij]);
	  meas_ret=meas(current[ij],w_current[ij],ij,true);
	  
	  if (verbose>=1) {
	    std::cout << "mcmc: " << ij << " " << w_current[ij] << std::endl;
	  }

	  // ------------------------------------------------
	  
	  // Increment iteration count
	  init_iters++;

	  // If we have a good point, stop the loop
	  if (w_current[ij]>0.0) {
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
      w_current[0]=func(nparams,current[0]);
      meas_ret=meas(current[0],w_current[0],0,true);
      if (verbose>=1) {
	std::cout << "mcmc: " << w_current[0] << std::endl;
      }
      if (w_current[0]<=0.0) {
	if (err_nonconv) {
	  O2SCL_ERR("Initial weight vanished.",o2scl::exc_einval);
	}
	return 2;
      }

      best=current[0];
      w_best=w_current[0];
      best_point(best,w_best);

      // Compute the initial Hastings proposal weight
      if (hg_mode>0) {
	q_current=hast.pdf(current[0]);
      }
    
    }

    if (meas_ret!=0) {
      O2SCL_ERR2("Measurement function returned non-zero value ",
		"for initial point",o2scl::exc_einval);
    }
    
    // ---------------------------------------------------

    // The MCMC is arbitrarily broken up into 20 'blocks', making
    // it easier to keep track of progress and ensure file updates
    size_t block_counter=0;

    bool main_done=false;
    size_t mcmc_iters=0;
    
    while (!main_done) {

      // Walker to move for smove
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
	      O2SCL_ERR("Failed to find suitable step.",
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
	    O2SCL_ERR("Sanity check in parameter step.",o2scl::exc_esanity);
	  }
	}
      
      }

      // ---------------------------------------------------
      // Compute next weight

      w_next=func(nparams,next);
      if (w_next>w_best) {
	best=next;
	w_best=w_next;
	best_point(best,w_best);
      }

      // ---------------------------------------------------
    
      // Test to ensure new point is good
      if (w_next<=0.0) {
	if (err_nonconv) {
	  O2SCL_ERR("Zero weight.",o2scl::exc_einval);
	}
	return 3;
      }

      bool accept=false;
      double r=gr.random();

      // Metropolis algorithm
      if (aff_inv) {
	if (r<pow(smove_z,((double)nwalk)-1.0)*w_next/w_current[ik]) {
	  accept=true;
	}
	if (verbose>=1) {
	  std::cout << "mcmc: " << w_current[ik] << " " << w_next << " "
		    << pow(smove_z,((double)nwalk)-1.0) << " ratio: "
		    << pow(smove_z,((double)nwalk)-1.0)*w_next/w_current[ik]
		    << " accept: " << accept << std::endl;
	}
      } else if (hg_mode>0) {
	if (r<w_next*q_current/w_current[0]/q_next) {
	  accept=true;
	}
	if (verbose>=1) {
	  std::cout << "mcmc: " << w_current[0] << " " << q_current
		    << " " << w_next << " " << q_next << " ratio: "
		    << w_next*q_current/w_current[0]/q_next
		    << " accept: " << accept << std::endl;
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
	  if (aff_inv) {
	    meas_ret=meas(next,w_next,ik,true);
	  } else {
	    meas_ret=meas(next,w_next,0,true);
	  }
	}
	  
	// Prepare for next point
	if (aff_inv) {
	  current[ik]=next;
	  w_current[ik]=w_next;
	} else {
	  current[0]=next;
	  w_current[0]=w_next;
	}
	  
      } else {
	    
	// Point was rejected
	n_reject++;

	// Repeat measurement of old point
	if (!warm_up) {
	  if (aff_inv) {
	    meas_ret=meas(current[ik],w_current[ik],ik,false);
	  } else {
	    meas_ret=meas(current[0],w_current[0],0,false);
	  }
	}

      }

      if (meas_ret!=0) {
	main_done=true;
	if (meas_ret!=-1 && err_nonconv) {
	  O2SCL_ERR("Measurement function returned error value.",
		    o2scl::exc_efailed);
	}
      }
      mcmc_iters++;
      
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
  template<class func_t=o2scl::multi_funct11,
    class measure_t=measure_funct, class vec_t=ubvector>
    class mcmc_table : public mcmc_base<func_t,measure_t,vec_t> {
    
  protected:
    
  /// Parameter names
  std::vector<std::string> param_names;
    
  /// Parameter units
  std::vector<std::string> param_units;
  //@}
    
  /// Main data table for Markov chain
  std::shared_ptr<o2scl::table_units<> > tab;
    
  /** \brief Create strings which contain column names and units
   */
  virtual void table_names_units(std::string &s, std::string &u) {
    s="mult weight ";
    u=". . ";
    for(size_t i=0;i<param_names.size();i++) {
      s+=this->param_names[i]+" ";
      if (this->param_units[i].length()>0) {
	u+=this->param_units[i]+" ";
      } else {
	u+=". ";
      }
    }
    return;
  }

  /** \brief Desc
   */
  virtual int mcmc_init() {

    // -----------------------------------------------------------
    // Init table
      
    std::string s, u;
    table_names_units(s,u);
    std::cout << "Here: " << s << std::endl;
    tab->line_of_names(s);
      
    {
      size_t ctr=0;
      std::string unit;
      std::istringstream is(u);
      while(is >> unit) {
	if (unit!=((std::string)".")) {
	  tab->set_unit(tab->get_column_name(ctr),unit);
	}
	ctr++;
      } 
      if (ctr!=tab->get_ncolumns()) {
	O2SCL_ERR("Column/unit alignment in mcmc_table::mcmc().",
		  o2scl::exc_esanity);
      }
    }
    return 0;
  }

  /** \brief Fill \c line with data for insertion into the table
   */
  virtual void fill_line(const vec_t &pars, double weight, 
			 std::vector<double> &line) {
    
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
    
  /** \brief Desc
   */
  virtual int add_line(const vec_t &pars, double weight,
			size_t ix, bool new_meas) {

    // Test to see if we need to add a new line of data or
    // increment the weight on the previous line
    if (tab->get_nlines()<=(this->nwalk-1) || new_meas==true) {
	
      std::vector<double> line;
      fill_line(pars,weight,line);
      tab->line_of_data(line.size(),line);
	
    } else if (tab->get_nlines()>0) {
	
      // Otherwise, just increment the multiplier on the previous line
      tab->set("mult",tab->get_nlines()-this->nwalk,
	       tab->get("mult",tab->get_nlines()-this->nwalk)+1.0);
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
