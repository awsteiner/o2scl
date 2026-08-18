/* 
   ───────────────────────────────────────────────────────────────────
   
   Copyright (C) 2010-2026, Edwin van Leeuwen and Andrew W. Steiner
   
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
#ifndef O2SCL_DIFF_EVO_ADAPT_H
#define O2SCL_DIFF_EVO_ADAPT_H

/** \file diff_evo_adapt.h
    \brief File defining \ref o2scl::diff_evo_adapt
*/

#include <vector>
#include <algorithm>
#include <limits>
#include <exception>

#include <o2scl/mm_funct.h>
#include <o2scl/diff_evo.h>

namespace o2scl {

  /** \brief Multidimensional minimization by the differential
      evolution method
       
      This class minimizes a function using differential evolution.
      This method is a genetic algorithm and as such works well for
      discontinuous problems, since it does not require the gradient
      of the function to be minimized.

      This is an adaptive version of \ref diff_evo as described in
      \verbatim embed:rst
      [Brest06]_ .
      \endverbatim

      By default this class uses the same steady-state ("in-place")
      update rule as \ref diff_evo, where a trial vector which
      improves on its parent replaces that parent in the population
      immediately, so that later agents in the same generation may
      already see some already-updated neighbors. Setting \ref
      generational to <tt>true</tt> switches to the classic
      generational update described in
      \verbatim embed:rst
      [Storn97]_ ,
      \endverbatim
      where every trial vector in a generation is evaluated against
      a fixed snapshot of the population taken at the start of that
      generation, and all accepted replacements are applied only
      after the whole generation has been evaluated. The two modes
      typically give very similar convergence behavior; the
      generational form is required for the safe, race-free
      parallel evaluation used by \ref o2scl::diff_evo_para, which
      always operates in generational mode.
  */
  template<class func_t=multi_funct, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class init_funct_t=mm_funct> 
    class diff_evo_adapt : public diff_evo<func_t, vec_t, init_funct_t> 
    {

    public:

    typedef boost::numeric::ublas::vector<double> ubvector;

    /// Probability of adjusting f (default 0.1)
    double tau_1;
    /// Probability of adjusting cr (default 0.1)
    double tau_2;

    /// \name Lower bound and range of F (defaults 0.1 and 0.9)
    //@{
    double fl, fr;
    //@}

    /** \brief If true, use the generational (synchronous) update
        rule instead of the default steady-state (in-place) rule
        (default false)

        See the class documentation above for the distinction
        between the two update rules. \ref o2scl::diff_evo_para
        forces this to <tt>true</tt> internally, since the
        generational rule is required for its parallel evaluation
        of a generation's trial vectors to be race-free.
    */
    bool generational;

    diff_evo_adapt() : diff_evo<func_t,vec_t,init_funct_t>() {
      tau_1=0.1;
      tau_2=0.1;
      fl=0.1;
      fr=0.9;
      generational=false;
    }

    /// Return string denoting type ("diff_evo_adapt")
    virtual const char *type() { return "diff_evo_adapt"; }
      
    /** \brief Calculate the minimum \c fmin of \c func w.r.t the
	array \c x of size \c nvar.
    */
    virtual int mmin(size_t nvar, vec_t &x0, double &fmin,
		     func_t &func) {

      // Keep track of number of generation without better solutions
      int nconverged=0;

      size_t pop_size_loc=select_pop_size(nvar);

      initialize_population(nvar,x0,pop_size_loc);
      fmins.resize(pop_size_loc);
      this->n_evals=0;

      // Set initial fmin
      for (size_t x=0;x<pop_size_loc;++x) {
	vec_t agent_x;
	agent_x.resize(nvar);
	for (size_t i=0;i<nvar;++i) {
	  agent_x[i]=this->population[x*nvar+i];
	}
	double fmin_x=0;
	fmin_x=func(nvar,agent_x);
	this->n_evals++;
	fmins[x]=fmin_x;
	if (x==0) {
	  fmin=fmin_x;
	  for (size_t i=0;i<nvar;++i) {
	    x0[i]=agent_x[i];
	  }
	} else if (fmin_x<fmin) {
	  fmin=fmin_x;
	  for (size_t i=0;i<nvar;++i) {
	    x0[i]=agent_x[i];
	  }
	}

      }

      int gen=0;
      while (gen<this->ntrial && nconverged <= ((int)this->nconv) &&
             (this->max_evals==0 || this->n_evals<this->max_evals)) {
	++nconverged;
	++gen;

        if (generational) {

          // Synchronous/generational update: evaluate the whole
          // generation against a fixed snapshot, then apply all
          // accepted replacements at once. See run_generation() and
          // the class documentation above.
          run_generation(nvar,pop_size_loc,func,fmin,x0,nconverged);

        } else {

	// For each agent x in the population do:
	for (size_t x=0;x<pop_size_loc;++x) {

	  std::vector<int> others;

	  // Create a copy agent_x and agent_y of the current agent vector
	  vec_t agent_x, agent_y;
	  agent_x.resize(nvar);
	  agent_y.resize(nvar);
	  for (size_t i=0;i<nvar;++i) {
	    agent_x[i]=this->population[x*nvar+i];
	    agent_y[i]=this->population[x*nvar+i];
	  }
	  // Value of f and cr for this agent
	  double f_x, cr_x;
	  if (this->gr.random()>=tau_1) {
	    f_x=variables[x*2];
	  } else {
	    f_x=fl+this->gr.random()*fr;
	  } if (this->gr.random()>=tau_2) {
	    cr_x=variables[x*2+1];
	  } else {
	    cr_x=this->gr.random();
	  }

	  // Pick three agents a, b, and c from the population at
	  // random, they must be distinct from each other as well
	  // as from agent x
	  others=this->pick_unique_agents(3,x,pop_size_loc);

	  // Pick a random index R in {1, ..., n}, where the highest
	  // possible value n is the dimensionality of the problem
	  // to be optimized.
	  size_t r=floor(this->gr.random()*nvar);

	  for (size_t i=0;i<nvar;++i) {
	    // Pick ri~U(0,1) uniformly from the open range (0,1)
	    double ri=this->gr.random();
	    // If (i=R) or (ri<CR) let yi=ai + F(bi - ci),
	    // otherwise let yi=xi
	    if (i==r || ri<cr_x) {
	      agent_y[i]=this->population[others[0]*nvar+i] +
		f_x*(this->population[others[1]*nvar+i]-
		     this->population[others[2]*nvar+i]);
	    }
	  }
	  // If (f(y) < f(x)) then replace the agent in the
	  // population with the improved candidate solution, that is,
	  // set x=y in the population
	  double fmin_y;

	  fmin_y=func(nvar,agent_y);
	  this->n_evals++;
	  if (fmin_y<fmins[x]) {
	    for (size_t i=0;i<nvar;++i) {
	      this->population[x*nvar+i]=agent_y[i];
	      fmins[x]=fmin_y;
	    }

	    variables[x*2]=f_x;
	    variables[x*2+1]=cr_x;

	    if (fmin_y<fmin) {
	      fmin=fmin_y;
	      for (size_t i=0;i<nvar;++i) {
		x0[i]=agent_y[i];
	      }
	      nconverged=0;
	    }
	  }

	  if (this->max_evals>0 && this->n_evals>=this->max_evals) break;

	}

        }

	if (this->verbose>0) {
	  this->print_iter(nvar,fmin,gen,x0,pop_size_loc,nconverged,
                           this->nconv);
	}
      }

      this->last_ntrial=gen;
      this->last_n_evals=this->n_evals;

      if (gen>=this->ntrial) {
	std::string str="Exceeded maximum number of iterations ("+
	itos(this->ntrial)+") in diff_evo_adapt::mmin().";
	O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      } else if (this->max_evals>0 && this->n_evals>=this->max_evals) {
	std::string str="Exceeded maximum number of function "
          "evaluations ("+itos(this->max_evals)+
          ") in diff_evo_adapt::mmin().";
	O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      }
      return 0;
    };

    /** \brief Print out iteration information
     */
    virtual void print_iter(size_t nvar, double fmin, 
			    int iter, vec_t &best_fit,
                            size_t pop_size_loc,
                            size_t nconverged_loc, size_t nconv_loc) {
      
      std::cout << type() << "::print_iter():\n  "
                << "Generation, min., n_converged:\n  "
                << iter << " of " << this->ntrial << ", " 
                << fmin << ", " << nconverged_loc << " of "
                << nconv_loc << std::endl;
      std::cout << "  Parameters: ";
      std::cout.setf(std::ios::showpos);
      for (size_t i=0;i<nvar;++i) {
	std::cout << best_fit[i] << " ";
      }
      std::cout.unsetf(std::ios::showpos);
      std::cout << std::endl;
      
      if (this->verbose>1) {
        std::cout << "  Population (index, parameters, minimum, "
                  << "weight, crossover):" << std::endl;
        for (size_t i=0;i<pop_size_loc;++i) {
          std::cout << "  " << i << ": ";
          std::cout.setf(std::ios::showpos);
          for (size_t j=0;j<nvar;++j) {
            std::cout << this->population[i*nvar+j] << " ";
          }
          std::cout << "fmin: " << fmins[i];
          std::cout.unsetf(std::ios::showpos);
          std::cout << " F: " << variables[i*2] <<
            " CR: " << variables[i*2+1] << std::endl;
        }
	char ch;
	std::cin >> ch;
      }
      return;
    }

    protected:

    /** \brief Vector containing the tunable variable F and CR
     */
    vec_t variables;
    
    /// Vector that keeps track of fmins values
    ubvector fmins;

    /** \brief Automatically select the population size based on the
        dimensionality of the problem, if the user hasn't set \ref
        o2scl::diff_evo::pop_size directly

        This is virtual so that \ref o2scl::diff_evo_para can also
        take the number of OpenMP threads into account, rounding
        the result up so that no thread is left idle during a
        generation's parallel evaluation loop.
    */
    virtual size_t select_pop_size(size_t nvar) {
      if (this->pop_size==0) {
        // Automatically select pop_size based on on dimensionality.
        return 10*nvar;
      }
      return this->pop_size;
    }

    /** \brief Select the (possibly self-adapting) values of f and
        cr to use for agent \c x in the current generation, using
        the specified random number generator

        This performs the jDE self-adaptation step: with
        probability \ref tau_1 (\ref tau_2), a new value of f (cr)
        is drawn at random; otherwise the agent's own value from
        \ref variables, inherited from the last generation in which
        it was updated, is reused unchanged.
    */
    virtual void adapt_f_cr(size_t x, double &f_x, double &cr_x,
                             rng<> &r_gen) {
      if (r_gen.random()>=tau_1) {
        f_x=variables[x*2];
      } else {
        f_x=fl+r_gen.random()*fr;
      }
      if (r_gen.random()>=tau_2) {
        cr_x=variables[x*2+1];
      } else {
        cr_x=r_gen.random();
      }
    }

    /** \brief Construct the "rand/1/bin" trial vector for agent \c
        x, reading positions from \c pop_src and using the
        specified random number generator

        The source population \c pop_src is passed explicitly
        (rather than always reading \ref o2scl::diff_evo::population
        directly) so that this can be safely called against a
        read-only snapshot of the population from a parallel
        generational update (see \ref run_generation() and \ref
        o2scl::diff_evo_para), as well as against the live
        population from the default steady-state update.
    */
    virtual void build_trial_vector(size_t nvar, size_t x,
                                     size_t pop_size_loc,
                                     const vec_t &pop_src,
                                     double f_x, double cr_x,
                                     vec_t &agent_y, rng<> &r_gen) {
      // Pick three agents a, b, and c from the population at
      // random, they must be distinct from each other as well
      // as from agent x
      std::vector<int> others=this->pick_unique_agents
        (3,x,pop_size_loc,r_gen);

      // Pick a random index R in {1, ..., n}, where the highest
      // possible value n is the dimensionality of the problem
      // to be optimized.
      size_t r=floor(r_gen.random()*nvar);

      for (size_t i=0;i<nvar;++i) {
        agent_y[i]=pop_src[x*nvar+i];
      }
      for (size_t i=0;i<nvar;++i) {
        // Pick ri~U(0,1) uniformly from the open range (0,1)
        double ri=r_gen.random();
        // If (i=R) or (ri<CR) let yi=ai + F(bi - ci),
        // otherwise let yi=xi
        if (i==r || ri<cr_x) {
          agent_y[i]=pop_src[others[0]*nvar+i] +
            f_x*(pop_src[others[1]*nvar+i]-
                 pop_src[others[2]*nvar+i]);
        }
      }
    }

    /** \brief Compute the trial vector and its function value for
        every agent in the population, given a fixed snapshot \c
        pop_src of the population taken at the start of the
        generation

        The default implementation here is a simple serial loop
        over all agents, using \ref gr as the random number
        generator and \c func as the function to evaluate. \ref
        o2scl::diff_evo_para overrides this method to evaluate the
        loop over agents in an OpenMP-parallel region, using an
        independent random number generator and (optionally) a
        separate copy of \c func_t per thread. This is the only
        method which needs to be overridden to parallelize
        generational differential evolution, since the population
        snapshot \c pop_src is read-only here and every output
        array is indexed by the agent number \c x, so no two agents
        write to the same location.

        A trial vector's construction can land \c func on an
        invalid or nonphysical point (differential evolution
        scatters trial vectors well away from any single guess, by
        design), and \c func may report that by having the error
        handler throw rather than returning some ordinary large
        value. Since generational evaluation of a whole population
        happens without any opportunity for the caller to intervene
        agent-by-agent -- and, in \ref o2scl::diff_evo_para, happens
        inside an OpenMP parallel region, where an exception that
        escapes uncaught is fatal (it cannot cross the parallel
        region boundary and calls <tt>std::terminate()</tt>
        immediately, and several threads doing so at once is exactly
        the "terminate called recursively" pattern) -- any exception
        thrown while evaluating a trial vector is caught here and
        treated as a very bad (but finite) fitness value instead of
        being allowed to propagate. This matches how black-box
        optimizers generally need to treat infeasible points, and
        means a single bad trial vector just gets discarded by \ref
        apply_trials() rather than aborting the whole minimization.
    */
    virtual void compute_trials(size_t nvar, size_t pop_size_loc,
                                 const vec_t &pop_src, func_t &func,
                                 std::vector<vec_t> &trial,
                                 ubvector &fmin_trial,
                                 vec_t &f_trial, vec_t &cr_trial) {
      for (size_t x=0;x<pop_size_loc;++x) {
        double f_x, cr_x;
        adapt_f_cr(x,f_x,cr_x,this->gr);
        f_trial[x]=f_x;
        cr_trial[x]=cr_x;
        trial[x].resize(nvar);
        build_trial_vector(nvar,x,pop_size_loc,pop_src,f_x,cr_x,
                            trial[x],this->gr);
        try {
          fmin_trial[x]=func(nvar,trial[x]);
        } catch (const std::exception &e) {
          fmin_trial[x]=std::numeric_limits<double>::max();
        }
        this->n_evals++;
      }
    }

    /** \brief Apply the trial vectors computed by \ref
        compute_trials() to the population, replacing any agent
        whose trial vector improved on it

        This always runs serially (over the OpenMP thread which
        calls it), since it mutates the shared population, \ref
        fmins, and \ref variables arrays. It is the same
        acceptance rule used by the default steady-state update,
        just applied all at once at the end of the generation
        rather than as each trial vector is computed.
    */
    virtual void apply_trials(size_t nvar, size_t pop_size_loc,
                               std::vector<vec_t> &trial,
                               ubvector &fmin_trial,
                               vec_t &f_trial, vec_t &cr_trial,
                               double &fmin, vec_t &x0,
                               int &nconverged) {
      for (size_t x=0;x<pop_size_loc;++x) {
        if (fmin_trial[x]<fmins[x]) {
          for (size_t i=0;i<nvar;++i) {
            this->population[x*nvar+i]=trial[x][i];
          }
          fmins[x]=fmin_trial[x];
          variables[x*2]=f_trial[x];
          variables[x*2+1]=cr_trial[x];
          if (fmin_trial[x]<fmin) {
            fmin=fmin_trial[x];
            for (size_t i=0;i<nvar;++i) {
              x0[i]=trial[x][i];
            }
            nconverged=0;
          }
        }
      }
    }

    /** \brief Perform one generational (synchronous) update of the
        population

        Takes a snapshot of the current population, computes trial
        vectors and function values for the whole generation
        against that snapshot (\ref compute_trials(), potentially
        in parallel), and then applies the accepted replacements
        (\ref apply_trials(), always serially).
    */
    virtual void run_generation(size_t nvar, size_t pop_size_loc,
                                 func_t &func, double &fmin, vec_t &x0,
                                 int &nconverged) {
      vec_t pop_snapshot=this->population;
      std::vector<vec_t> trial(pop_size_loc);
      ubvector fmin_trial(pop_size_loc);
      vec_t f_trial(pop_size_loc), cr_trial(pop_size_loc);
      compute_trials(nvar,pop_size_loc,pop_snapshot,func,trial,
                     fmin_trial,f_trial,cr_trial);
      apply_trials(nvar,pop_size_loc,trial,fmin_trial,f_trial,cr_trial,
                   fmin,x0,nconverged);
    }

    /** \brief Initialize a population of random agents
     */
      virtual int initialize_population(size_t nvar, vec_t &x0,
                                        size_t pop_size_loc) {
        
      this->population.resize(nvar*pop_size_loc);
      variables.resize(2*pop_size_loc);
      if (this->rand_init_funct==0) {
	for(size_t i=0;i<pop_size_loc;i++) {
	  for(size_t j=0;j<nvar;j++) {
            if (this->use_initial_point && i==0) {
              this->population[i*nvar+j]=x0[j];
            } else {
              double stepj=this->step[j%this->step.size()];
              this->population[i*nvar+j]=x0[j]-stepj/2.0+
                stepj*this->gr.random();
            }
          }
	  variables[i*2]=fl+this->gr.random()*fr;
	  variables[i*2+1]=this->gr.random();
	}
      } else {
	for (size_t i=0;i<pop_size_loc;++i) {
	  vec_t y(nvar);
	  (*this->rand_init_funct)(nvar,x0,y);
	  for (size_t j=0;j<nvar;++j) {
	    this->population[i*nvar+j]=y[j];
	  }
	  variables[i*2]=fl+this->gr.random()*fr;
	  variables[i*2+1]=this->gr.random();
	}
      }
      return 0;
    }

    private:
      
      diff_evo_adapt<func_t,vec_t,init_funct_t>
      (const diff_evo_adapt<func_t,vec_t,init_funct_t> &);
      
      diff_evo_adapt<func_t,vec_t,init_funct_t> &operator=
      (const diff_evo_adapt<func_t,vec_t,init_funct_t>&);
      
    };

}

#endif
