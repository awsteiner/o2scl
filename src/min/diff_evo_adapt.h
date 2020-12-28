 /* 
   -------------------------------------------------------------------
   
   Copyright (C) 2010-2021, Edwin van Leeuwen and Andrew W. Steiner
   
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
#ifndef O2SCL_DIFF_EVO_ADAPT_H
#define O2SCL_DIFF_EVO_ADAPT_H

/** \file diff_evo_adapt.h
    \brief File defining \ref o2scl::diff_evo_adapt
*/

#include <vector>
#include <algorithm>

#include <o2scl/rng_gsl.h>
#include <o2scl/mm_funct.h>

#include <o2scl/diff_evo.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multidimensional minimization by the differential
      evolution method
       
      This class minimizes a function using differential evolution.
      This method is a genetic algorithm and as such works well for
      non continuous problems, since it does not require the gradient
      of the function to be minimized.

      This is an adaptive version of \ref diff_evo as described in
      \verbatim embed:rst
      [Brest06]_ .
      \endverbatim
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

    diff_evo_adapt() : diff_evo<func_t,vec_t,init_funct_t>() {
      tau_1 = 0.1;
      tau_2 = 0.1;
      fl = 0.1;
      fr = 0.9;
    }

    /** \brief Calculate the minimum \c fmin of \c func w.r.t the 
	array \c x of size \c nvar.
    */
    virtual int mmin(size_t nvar, vec_t &x0, double &fmin,
		     func_t &func) {

      // Keep track of number of generation without better solutions
      int nconverged = 0;
      if (this->pop_size==0) {
	// Automatically select pop_size dependent on dimensionality.
	this->pop_size = 10*nvar;
      }

      initialize_population( nvar, x0 );
      fmins.resize(this->pop_size);

      // Set initial fmin
      for (size_t x = 0; x < this->pop_size; ++x) {
	vec_t agent_x;
	agent_x.resize(nvar);
	for (size_t i = 0; i < nvar; ++i) {
	  agent_x[i] = this->population[x*nvar+i];
	}
	double fmin_x = 0;
	fmin_x=func(nvar,agent_x);
	fmins[x]=fmin_x;
	if (x==0) {
	  fmin = fmin_x;
	  for (size_t i = 0; i<nvar; ++i) {
	    x0[i] = agent_x[i];
	  }
	} else if (fmin_x<fmin) {
	  fmin = fmin_x;
	  for (size_t i = 0; i<nvar; ++i) {
	    x0[i] = agent_x[i];
	  }
	}

      }

      int gen = 0;
      while (gen < this->ntrial && nconverged <= ((int)this->nconv)) {
	++nconverged;
	++gen;

	// For each agent x in the population do: 
	for (size_t x = 0; x < this->pop_size; ++x) {

	  std::vector<int> others;

	  // Create a copy agent_x and agent_y of the current agent vector
	  vec_t agent_x, agent_y;
	  agent_x.resize(nvar);
	  agent_y.resize(nvar);
	  for (size_t i = 0; i < nvar; ++i) {
	    agent_x[i] = this->population[x*nvar+i];
	    agent_y[i] = this->population[x*nvar+i];
	  }
	  // Value of f and cr for this agent
	  double f_x, cr_x;
	  if (this->gr.random() >= tau_1) {
	    f_x = variables[x*2];
	  } else {
	    f_x = fl+this->gr.random()*fr;
	  } if (this->gr.random() >= tau_2) {
	    cr_x = variables[x*2+1];
	  } else {
	    cr_x = this->gr.random();
	  }
                            
	  // Pick three agents a, b, and c from the population at 
	  // random, they must be distinct from each other as well 
	  // as from agent x
	  others = this->pick_unique_agents( 3, x );

	  // Pick a random index R in {1, ..., n}, where the highest 
	  // possible value n is the dimensionality of the problem 
	  // to be optimized.
	  size_t r = floor(this->gr.random()*nvar);

	  for (size_t i = 0; i < nvar; ++i) {
	    // Pick ri~U(0,1) uniformly from the open range (0,1)
	    double ri = this->gr.random();
	    // If (i=R) or (ri<CR) let yi = ai + F(bi - ci), 
	    // otherwise let yi = xi
	    if (i == r || ri < cr_x) {
	      agent_y[i] = this->population[others[0]*nvar+i] + 
		f_x*(this->population[others[1]*nvar+i]-
		     this->population[others[2]*nvar+i]);
	    }
	  }
	  // If (f(y) < f(x)) then replace the agent in the 
	  // population with the improved candidate solution, that is, 
	  // set x = y in the population
	  double fmin_y;
                            
	  fmin_y=func(nvar,agent_y);
	  if (fmin_y<fmins[x]) {
	    for (size_t i = 0; i < nvar; ++i) {
	      this->population[x*nvar+i] = agent_y[i];
	      fmins[x] = fmin_y;
	    }
                                
	    variables[x*2] = f_x;
	    variables[x*2+1] = cr_x;

	    if (fmin_y<fmin) {
	      fmin = fmin_y;
	      for (size_t i = 0; i<nvar; ++i) {  
		x0[i] = agent_y[i];
	      }
	      nconverged = 0;
	    }
	  }

	}
	if (this->verbose > 0) {
	  this->print_iter( nvar, fmin, gen, x0 );
	}
      }
      
      this->last_ntrial=gen;
      
      if (gen>=this->ntrial) {
	std::string str="Exceeded maximum number of iterations ("+
	itos(this->ntrial)+") in diff_evo_adapt::mmin().";
	O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      }
      return 0;
    };

    /** \brief Print out iteration information
     */
    virtual void print_iter(size_t nvar, double fmin, 
			    int iter, vec_t &best_fit) {
      
      std::cout << "Generation " << iter << std::endl;
      std::cout << "Fmin: " << fmin << std::endl;
      std::cout << "Parameters: ";
      for (size_t i=0; i<nvar; ++i) {
	std::cout << best_fit[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Population: " << std::endl;
      for (size_t i=0; i<this->pop_size; ++i ) {
	std::cout << i << ": ";
	for (size_t j = 0; j<nvar; ++j ) {
	  std::cout << this->population[i*nvar+j] << " ";
	}
	std::cout << "fmin: " << fmins[i] << 
	  " F: " << variables[i*2] <<
	  " CR: " << variables[i*2+1] << std::endl;
      }
      if (this->verbose>1) {
	char ch;
	std::cin >> ch;
      }
      return;
    }

 
#ifndef DOXYGEN_INTERNAL

    protected:

    /** \brief Vector containing the tunable variable F and CR
     */
    vec_t variables;
    
    /// Vector that keeps track of fmins values
    ubvector fmins;

    /** \brief Initialize a population of random agents
     */
    virtual int initialize_population( size_t nvar, vec_t &x0 ) {
      this->population.resize(nvar*this->pop_size);
      variables.resize(2*this->pop_size);
      if (this->rand_init_funct==0) {
	for(size_t i=0;i<this->pop_size;i++) {
	  for(size_t j=0;j<nvar;j++) {
	    double stepj=this->step[j%this->step.size()];
	    this->population[i*nvar+j]=x0[j]-stepj/2.0+stepj*this->gr.random();
	  }
	  variables[i*2] = fl + this->gr.random()*fr;
	  variables[i*2+1] = this->gr.random();
	}
      } else {
	for (size_t i = 0; i < this->pop_size; ++i) {
	  vec_t y(nvar);
	  (*this->rand_init_funct)(nvar,x0,y);
	  for (size_t j = 0; j < nvar; ++j) {
	    this->population[ i*nvar+j ] = y[j];
	  }
	  variables[i*2] = fl + this->gr.random()*fr;
	  variables[i*2+1] = this->gr.random();
	}
      }
      return 0;
    }

    private:

    diff_evo_adapt<func_t,vec_t,init_funct_t>
    (const diff_evo_adapt<func_t,vec_t,init_funct_t> &);
    diff_evo_adapt<func_t,vec_t,init_funct_t> &operator=
    (const diff_evo_adapt<func_t,vec_t,init_funct_t>&);

#endif

    };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
