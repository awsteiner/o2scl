/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

/* Example: ex_methast.cpp
   -------------------------------------------------------------------

   An example to demonstrate generation of an arbitrary distribution
   through the creation of a Markov chain using the
   Metropolis-Hastings algorithm.

*/
  
#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hist.h>
#include <o2scl/rng.h>
#include <o2scl/expval.h>
#include <o2scl/inte_qag_gsl.h>

using namespace std;
using namespace o2scl;

/*
  An unusual two-peaked probability distribution. 
*/
double expsin_fun(double x) {
  return exp(-x*x)*(sin(x-1.4)+1.0);
}

// Same function but now weighted with x
double expsin_fun2(double x) {
  return x*expsin_fun(x);
}

// The proposal distribution
double weight(double x) {
  return exp(-(x-1)*(x-1)/0.5)*2.0/3.0+exp(-(x+1)*(x+1)/0.5)/3.0;
}

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  // The square root of the number of iterations
  static const size_t N=1000;
  
  // The random number generator
  rng<> gr;
  if (argc>=2) {
    // We have to use the full o2scl::stoi() specification to 
    // avoid confusion with std::stoi()
    gr.set_seed(o2scl::stoi(argv[1]));
  } 

  // Run twice, first time just with simple Metropolis
  for(size_t ell=0;ell<2;ell++) {

    bool hastings;
    if (ell==0) hastings=false;
    else hastings=true;

    // Construct a histogrmam of the distribution. Initialize the grid 
    // spacing with a unform grid between -5 and 5. 
    hist h;
    h.set_bin_edges(uniform_grid_end_width<>(-5.0,5.0,0.1));
    h.clear_wgts();

    // Compute the average value of the distribution, and use
    // a expval_scalar object to estimate the uncertainty
    static const size_t M=20;
    expval_scalar sev(M,N*N/M);

    // An initial point and the initial weight
    double x_old=0.01;
    double w_old=expsin_fun(x_old);
    double Q_old=weight(x_old);

    // GSL RNG
    const gsl_rng_type *T=gsl_rng_default;
    gsl_rng *grng=gsl_rng_alloc(T);
    gsl_rng_env_setup();

    for(size_t k=0;k<N;k++) {
      for(size_t i=0;i<N;i++) {

	// Perform a step and compute the new weight
	double x_new;

	if (!hastings) {
	  x_new=gr.random()*10.0-5.0;
	} else {
	  do {
	    if (gr.random()<1.0/3.0) {
	      x_new=gsl_ran_gaussian(grng,0.5)-1.0;
	    } else {
	      x_new=gsl_ran_gaussian(grng,0.5)+1.0;
	    }
	  } while (fabs(x_new)>5.0);
	}

	double w_new=expsin_fun(x_new);
	double Q_new=weight(x_new);

	double r=gr.random();
	if (hastings) {
	  // The Metrpolis-Hastings algorithm
	  if (r<w_new*Q_old/w_old/Q_new && x_new>-5.0 && x_new<5.0) {
	    x_old=x_new;
	    w_old=w_new;
	    Q_old=Q_new;
	  }
	} else {
	  // The Metropolis algorithm
	  if (r<w_new/w_old && x_new>-5.0 && x_new<5.0) {
	    x_old=x_new;
	    w_old=w_new;
	  }
	}

	// Update the histogram
	h.update(x_old);

	// Update the average value
	sev.add(x_old);

      }
    }
  
    // Normalize by the value of the distribution at x=1.05
    double norm=h.get_wgt(1.05)/expsin_fun(1.05);

    // Output the generated distribution and compare to the real
    // distribution.

    for(size_t i=0;i<h.size();i++) {
    
      double xx=h.get_rep_i(i);

      // Limit output to every third point
      if (i%3==0) {
	cout.width(2);
	cout << i << " ";
	cout.setf(ios::showpos);
	cout << xx << " "; 
	cout.unsetf(ios::showpos);
	cout << h.get_wgt_i(i)/norm << " " 
	     << expsin_fun(xx) << " " << weight(xx)/4.0 << " " 
	     << fabs(h.get_wgt_i(i)/norm-expsin_fun(xx))/expsin_fun(xx) << endl;
      }
    
      // Only test regions where we'll have enough statistics
      if (expsin_fun(xx)>5.0e-3) {
	t.test_rel(h.get_wgt_i(i)/norm,expsin_fun(xx),3.0e-1,"dist");
      }
    
    }
    cout << endl;

    // Result from direct integration
    funct ff=expsin_fun;
    funct ff2=expsin_fun2;
    inte_qag_gsl<> qag;
    double exact=qag.integ(ff2,-5.0,5.0)/qag.integ(ff,-5.0,5.0);
    cout << "Exact: " << exact << endl;

    // Output the average value and its uncertainty
    double avg, std_dev, avg_err;
    sev.current_avg(avg,std_dev,avg_err);
    cout << "Simul: " << avg << " " << avg_err << " " 
	 << fabs(avg-exact) << endl;

    t.test_abs(avg,exact,3.0*avg_err,"Average value");
    cout << endl;

  }
  
  t.report();

  return 0;
}
// End of example
