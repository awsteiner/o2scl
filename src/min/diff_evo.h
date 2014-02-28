/* 
  -------------------------------------------------------------------
   
  Copyright (C) 2010-2014, Edwin van Leeuwen
   
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
#ifndef O2SCL_DIFF_EVO_H
#define O2SCL_DIFF_EVO_H

/** \file diff_evo.h
    \brief File defining \ref o2scl::diff_evo
*/

#include <vector>
#include <algorithm>

#include <o2scl/rng_gsl.h>
#include <o2scl/mmin.h>
#include <o2scl/mm_funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multidimensional minimization by the differential
      evolution method
      
      This class minimizes a function using differential evolution.
      This method is a genetic algorithm and as such works well for
      non continuous problems, since it does not rely on a gradient of
      the function that is being mind.
      
      The method starts by initializing a random population of
      candidate parameters. To do this the user needs to define a
      function to create these random parameters, which can be
      provided using \ref set_init_function().
      
      After the initial population is created the algorithm will
      repeat a number of standard steps until a solution is found or
      the maximum number of iterations is reached. Based on 
      \ref Storn97.

      \note The constructor sets \ref o2scl::mmin_base::ntrial to 1000 .

      If the population converges prematurely, then \ref diff_evo::f
      and \ref pop_size should be increased.
  */
    template<class func_t=multi_funct<>, 
      class vec_t=boost::numeric::ublas::vector<double> , 
      class init_funct_t=mm_funct_fptr<vec_t > > class diff_evo : 
      public mmin_base<func_t,func_t,vec_t> {
      
    public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;

    /** \brief Population size (default 0)

	Should be at least 4. Typically between \f$ 5 d \f$ and \f$ 10
	d \f$ where \f$ d \f$ is the dimensionality of the problem.

	If this is 0 (the default), then it is set by mmin to be 
	equal to \f$ 10 d \f$ .
    */
    size_t pop_size;

    /** \brief The number of generations without a better fit before we
	assume that the algorithm has converged (default 25)
    */
    size_t nconv;

    /** \brief Differential weight (default 0.75)
	
	A parameter which controls the amplification of the
	differential variation. Usually between 0 and 2.
    */
    double f;

    /** \brief Crossover probability (default 0.8)
	
	Usually between 0 and 1. 
    */
    double cr;

    diff_evo() {
      this->ntrial=1000;
      f = 0.75;
      cr = 0.8;
      rand_init_funct = 0;
      pop_size = 0;
      nconv = 25;
    }

    virtual ~diff_evo() {
    }

    /** \brief Set the function that is used to produce random 
	init variables
	
	REQUIRED
	
	The init function is called in the beginning to fill 
	the population with random individuals, so it is best 
	to make this cover the part of the parameter space you 
	are interested in. The method will find solutions outside 
	this parameter space, but choosing a good init function will 
	help finding solutions faster.
    */
    virtual void set_init_function( init_funct_t &function ) {
      rand_init_funct = &function;
    }

    /** \brief Calculate the minimum \c fmin of \c func w.r.t the 
	array \c x of size \c nvar.

	Initialize all agents x with random positions in the
	search-space. Until a termination criterion is met (e.g.
	number of iterations performed, or adequate fitness reached),
	repeat the following: For each agent x in the population do:
	Pick three agents a, b, and c from the population at random,
	they must be distinct from each other as well as from agent x
	Pick a random index {1, ..., n}, where the highest possible
	value n is the dimensionality of the problem to be optimized.
	Compute the agent's potentially new position y = [y1, ..., yn]
	by iterating over each i {1, ..., n} as follows: Pick
	ri~U(0,1) uniformly from the open range (0,1) If (i=R) or
	(ri<CR) let yi = ai + F(bi - ci), otherwise let yi = xi If
	(f(y) < f(x)) then replace the agent in the population with
	the improved candidate solution, that is, set x = y in the
	population.
                  
	Pick the agent from the population that has the lowest fmin
	and return it as the best found candidate solution.
    */
    virtual int mmin(size_t nvar, vec_t &x0, double &fmin, func_t &func) {

      // Keep track of number of generation without better solutions
      size_t nconverged = 0;

      if (pop_size==0) {
	// Automatically select pop_size dependent on dimensionality.
	pop_size = 10*nvar;
      }

      initialize_population( nvar, x0 );
      
      fmins.resize(pop_size);

      // Set initial fmin
      for (size_t x = 0; x < pop_size; ++x) {
	vec_t agent_x;
	agent_x.resize(nvar);
	for (size_t i = 0; i < nvar; ++i) {
	  agent_x[i] = population[x*nvar+i];
	}
	double fmin_x = 0;
	fmin_x=func(nvar,agent_x);
	fmins[x]=fmin_x;
	if (x==0) {
	  fmin = fmin_x;
	  for (size_t i = 0; i<nvar; ++i)  
	    x0[i] = agent_x[i];
	  //x0 = agent_x;
	} else if (fmin_x<fmin) {
	  fmin = fmin_x;
	  for (size_t i = 0; i<nvar; ++i)  
	    x0[i] = agent_x[i];
	  //x0 = agent_x;
	}
      }

      int gen = 0;
      while (gen < this->ntrial && nconverged <= nconv) {
	     
	++nconverged;
	++gen;

	// For each agent x in the population do: 
	for (size_t x = 0; x < pop_size; ++x) {

	  std::vector<int> others;

	  // Create a copy agent_x and agent_y of the current agent
	  // vector
	  vec_t agent_x, agent_y;
	  agent_x.resize(nvar);
	  agent_y.resize(nvar);
	  for (size_t i = 0; i < nvar; ++i) {
	    agent_x[i] = population[x*nvar+i];
	    agent_y[i] = population[x*nvar+i];
	  }
                            
	  // Pick three agents a, b, and c from the population at 
	  // random, they must be distinct from each other as well as
	  // from agent x
	  others = pick_unique_agents( 3, x );

	  // Pick a random index R in {1, ..., n}, where the highest 
	  // possible value n is the dimensionality of the problem 
	  // to be optimized.
	  size_t r = floor(gr.random()*nvar);

	  for (size_t i = 0; i < nvar; ++i) {
	    // Pick ri~U(0,1) uniformly from the open range (0,1)
	    double ri = gr.random();
	    // If (i=R) or (ri<CR) let yi = ai + F(bi - ci), otherwise 
	    // let yi = xi
	    if (i == r || ri < cr) {
	      agent_y[i] = population[others[0]*nvar+i] + 
		f*(population[others[1]*nvar+i]-
		   population[others[2]*nvar+i]);
	    }
	  }
	  // If (f(y) < f(x)) then replace the agent in the population 
	  // with the improved candidate solution, that is, set x = y 
	  // in the population
	  double fmin_y;
                            
	  fmin_y=func(nvar,agent_y);
	  if (fmin_y<fmins[x]) {
	    for (size_t i = 0; i < nvar; ++i) {
	      population[x*nvar+i] = agent_y[i];
	      fmins[x] = fmin_y;
	    }
	    if (fmin_y<fmin) {
	      fmin = fmin_y;
	      for (size_t i = 0; i<nvar; ++i) {  
		x0[i] = agent_y[i];
	      }
	      nconverged = 0;
	    }
	  }

	}
	if (this->verbose > 0)
	  this->print_iter( nvar, fmin, gen, x0 );
      }

      if(gen>=this->ntrial) {
	std::string str="Exceeded maximum number of iterations ("+
	itos(this->ntrial)+") in diff_evo::mmin().";
	O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      }

      return 0;
    };

    /** \brief Print out iteration information.

	\comment
	Depending on the value of the variable verbose, this prints
	out the iteration information. If verbose=0, then no
	information is printed, while if verbose>1, then after each
	iteration, the present values of x and y are output to
	std::cout along with the iteration number. If verbose>=2 then
	each iteration waits for a character.
	\endcomment
    */
    virtual void print_iter( size_t nvar, double fmin, 
			     int iter, vec_t &best_fit ) {
      std::cout << "Generation " << iter << std::endl;
      std::cout << "Fmin: " << fmin << std::endl;
      std::cout << "Parameters: ";
      for (size_t i=0; i<nvar; ++i) {
	std::cout << best_fit[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Population: " << std::endl;
      for (size_t i=0; i<pop_size; ++i ) {
	std::cout << i << ": ";
	for (size_t j = 0; j<nvar; ++j ) {
	  std::cout << population[i*nvar+j] << " ";
	}
	std::cout << "fmin: " << fmins[i] << std::endl;
      }
    }

#ifndef DOXYGEN_INTERNAL

    protected:

    /** \brief Vector containing the population.
	
	For now using one long vector with all agents after each other
    */
    vec_t population;

    /// Vector that keeps track of fmins values
    ubvector fmins;

    /** \brief Function that is used to produce random init variables
	
	This function is used to fill the population with random agents
    */
    init_funct_t *rand_init_funct;

    /// Random number generator
    rng_gsl gr;

    /** \brief Initialize a population of random agents
     */
    virtual int initialize_population( size_t nvar, vec_t &x0 ) {
      if (rand_init_funct==0) {
	O2SCL_ERR_RET("No initialization function provided.",
		      exc_ebadfunc);

      }
      population.resize(nvar*pop_size );
      for (size_t i = 0; i < pop_size; ++i) {
	vec_t y(nvar);
	(*rand_init_funct)( nvar, x0, y );
	for (size_t j = 0; j < nvar; ++j) {
	  population[ i*nvar+j ] = y[j];

	}
      }
      return 0;
    }

    /** \brief Pick number of unique agent id's
	
	Unique from x and each other
	
	Uses the Fisher-Yates algorithm.  
	
    */
    virtual std::vector<int> pick_unique_agents( int nr, size_t x ) {
      std::vector<int> ids;
      std::vector<int> agents;
      // Fill array with ids
      for (size_t i=0; i<pop_size-1; ++i){
	if (i<x) {
	  ids.push_back( i );
	} else {
	  ids.push_back( i+1 );
	}
      }
      // Shuffle according to Fisher-Yates
      for (size_t i=ids.size()-1; i>ids.size()-nr-1; --i) {
	int j = round(gr.random()*i);
	std::swap( ids[i], ids[j] );
      }
      for (size_t i=ids.size()-1; i>ids.size()-nr-1; --i) {
	agents.push_back( ids[i] );
      }
      return agents;
    }

#endif

#ifndef DOXYGEN_INTERNAL

  private:

    diff_evo<func_t,vec_t,init_funct_t>
    (const diff_evo<func_t,vec_t,init_funct_t> &);
    diff_evo<func_t,vec_t,init_funct_t> &operator=
    (const diff_evo<func_t,vec_t,init_funct_t>&);

#endif

    }; 

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
