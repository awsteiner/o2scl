/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

/* Example: ex_markov.cpp
   -------------------------------------------------------------------

   An example which demonstrates (i) the generation of an arbitrary
   distribution through the creation of a Markov chain using the
   Metropolis algorithm, and (ii) the evaluation of a integral over
   that distribution using an expval object.
*/
  
#include <iostream>
#include <cmath>
#include <o2scl/test_mgr.h>
#include <o2scl/hist.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/expval.h>
#include <o2scl/inte_qag_gsl.h>

using namespace std;
using namespace o2scl;

/*
  An unusual two-peaked probability distribution. This distribution is
  is strongly bimodal, and thus the algorithm has to step over a
  region where the probability distribution is small. 
*/
double expsin_fun(double x) {
  return exp(-x*x)*(sin(x-1.4)+1.0);
}

// Same function but now weighted with x
double expsin_fun2(double x) {
  return x*exp(-x*x)*(sin(x-1.4)+1.0);
}

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  // The square root of the number of iterations
  static const size_t N=1000;

  // The stepsize for the Metropolis algorithm. If this stepsize was
  // too small, the algorithm would fail to step between the two
  // peaks and give incorrect results.
  double step=1.0;

  // The random number generator
  rng_gsl gr;
  if (argc>=2) {
    // We have to use the full o2scl::stoi() specification to 
    // avoid confusion with std::stoi()
    gr.set_seed(o2scl::stoi(argv[1]));
  }

  // Construct histograms of the distribution. Initialize the grid 
  // spacing with a unform grid between -5 and 5. 
  hist h, h2;
  h.set_bin_edges(uniform_grid_end_width<>(-5.0,5.0,0.1));
  h2.set_bin_edges(uniform_grid_end_width<>(-5.0,5.0,0.1));

  // Compute the average value of the distribution, and use
  // a expval_scalar object to estimate the uncertainty
  static const size_t M=20;
  expval_scalar sev(M,N*N/M);

  // Store the coordinates to get the unblocked average value
  std::vector<double> store;

  // An initial point and the initial weight
  double x_curr=0.01;
  double w_curr=expsin_fun(x_curr);

  double dev=0.0, dev2=0.0;

  for(size_t k=0;k<N;k++) {
    for(size_t i=0;i<N;i++) {

      // Perform a step and compute the new weight
      double r=gr.random();
      double x_new=x_curr+r*step*2.0-step;
      double w_new=expsin_fun(x_new);

      // The Metrpolis algorithm
      double r2=gr.random();
      if (r2<w_new/w_curr && x_new>-5.0 && x_new<5.0) {
	x_curr=x_new;
	w_curr=w_new;
      }

      /*
	There are two ways to do this, one is to just count up the
	number of hits in each bin, the second is to record the
	function value (the weight) in each bin, and divide by the
	number of counts afterwards to get the average.
      */
      
      // Update the histogram
      h.update(x_curr);

      // Update the second histogram with the weight
      h2.update(x_curr,w_curr);

      // Update the average value
      sev.add(x_curr);

      // Compute the average value using an simple unblocked average
      store.push_back(x_curr);
    }
  }

  // Normalize by the value of the distribution at x=1.05
  double norm=h.get_wgt(1.05)/expsin_fun(1.05);

  // Normalize the second histogram 
  for(size_t i=0;i<h.size();i++) {
    if (h[i]>0.0) {
      h2[i]/=h[i];
    } else {
      h2[i]=0.0;
    }
  }

  // Output the generated distribution and compare to the real
  // distribution.

  for(size_t i=0;i<h.size();i++) {
    
    double xx=h.get_rep_i(i);

    // Limit output to every third point
    if (i%3==0) {
      cout.width(2);
      cout.setf(ios::left);
      cout << i << " ";
      cout.unsetf(ios::left);
      cout.setf(ios::showpos);
      cout << xx << " "; 
      cout.unsetf(ios::showpos);
      cout << h[i]/norm << " "; 
      cout << h2[i] << " " << expsin_fun(xx) << endl;
    }

    dev+=fabs(h[i]/norm-expsin_fun(xx))/expsin_fun(xx);
    dev2+=fabs(h2[i]-expsin_fun(xx))/expsin_fun(xx);

    // Only test regions where we'll have enough statistics
    if (expsin_fun(xx)>5.0e-3) {
      t.test_rel(h.get_wgt_i(i)/norm,expsin_fun(xx),0.1,"dist");
    }
    
  }
  cout << endl;

  // Compare the deviations of the first and the second
  // histogram
  cout << "Avg. devs.: " << dev/((double)h.size()) << " " 
       << dev2/((double)h.size()) << endl;
  cout << endl;

  cout.precision(4);

  // Result from direct integration
  funct_fptr ff(expsin_fun);
  funct_fptr ff2(expsin_fun2);
  inte_qag_gsl<> qag;
  double exact=qag.integ(ff2,-5.0,5.0)/qag.integ(ff,-5.0,5.0);
  cout << "Exact result: " << exact << endl;

  // Output the average value and its uncertainty
  double avg, std_dev, avg_err;
  sev.current_avg(avg,std_dev,avg_err);
  cout << "Avg, err, |avg-exact| inc. autocorr.: " 
       << avg << " " << avg_err << " " << fabs(avg-exact) << endl;

  t.test_abs(avg,exact,3.0*avg_err,"Average value");

  // Compute average and uncertainty ignoring autocorrelations
  double avg2, std_dev2, avg_err2;
  avg2=vector_mean(store.size(),store);
  std_dev2=vector_variance(store.size(),store);
  avg_err2=vector_variance(store.size(),store)/
    sqrt(((double)store.size()));
  cout << "Avg, err, |avg-exact| w/o autocorr. : " 
       << avg2 << " " << avg_err2 << " " << fabs(avg2-exact) << endl;

  // Show that expval and vector_mean() report the same average
  t.test_rel(avg,avg2,1.0e-12,"Consistency");

  cout << endl;

  t.report();
  
  return 0;
}
// End of example
