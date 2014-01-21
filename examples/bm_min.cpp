/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2014, Andrew W. Steiner

  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------
*/
/*
  Compare O2scl minimizers with Boost
*/

#include <cmath>
#include <functional>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/min_brent_gsl.h>
#include <o2scl/min_cern.h>
#include <o2scl/min_brent_boost.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/hist.h>
#include <o2scl/expval.h>

#include <boost/math/tools/minima.hpp>

using namespace std;
using namespace o2scl;

size_t cnt;

double func(double x) {
  cnt++;
  return -exp(-pow(x-0.5,4.0));
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  const size_t N=500;

  rng_gsl gr;

  // Choose a logarithmic grid
  hist h;
  h.set_bin_edges(uniform_grid_log_end_width<>(1.0e-16,1.0e-2,10.0));
  
  // Store timing data in expect_val objects with 10 blocks
  vector<expval_scalar> time_mbg(h.size());
  vector<expval_scalar> time_mcn(h.size());
  vector<expval_scalar> time_mbb(h.size());
  vector<expval_scalar> evals_mbg(h.size());
  vector<expval_scalar> evals_mcn(h.size());
  vector<expval_scalar> evals_mbb(h.size());
  for(size_t i=0;i<h.size();i++) {
    time_mbg[i].set_blocks(10,1);
    time_mcn[i].set_blocks(10,1);
    time_mbb[i].set_blocks(10,1);
    evals_mbg[i].set_blocks(10,1);
    evals_mcn[i].set_blocks(10,1);
    evals_mbb[i].set_blocks(10,1);
  }

  vector<bool> has_data1(h.size());
  vector<bool> has_data2(h.size());
  vector<bool> has_data3(h.size());
  
  funct_fptr ff(func);
  std::pair<double,double> res;
  std::function<double (double)> f2=func;
  double x0,x1,x2,min;

  cout << "MBG   ";

  min_brent_gsl<> mbg;

  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    mbg.tol_abs=pow(10.0,-4-((double)j)/18.0);
    mbg.tol_rel=pow(10.0,-4-((double)j)/18.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();

      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=1.9+gr.random()*0.01;
	x2=2.0+gr.random()*0.01;

	mbg.min_bkt(x1,x0,x2,min,ff);
      }
      
      double acc=sqrt(pow(x1-0.5,2.0)+pow(min+1.0,2.0));
      if (acc>1.0e-16 && acc<1.0e-2) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_mbg[ix].add(t1);
	evals_mbg[ix].add(((double)(cnt-1)));
	if (!has_data1[ix]) has_data1[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "MCN   ";

  min_cern<> mcn;

  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    mcn.tol_abs=pow(10.0,-4-((double)j)/18.0);
    mcn.tol_rel=pow(10.0,-4-((double)j)/18.0);

    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      
      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=1.9+gr.random()*0.01;
	x2=2.0+gr.random()*0.01;

	mcn.min_bkt(x1,x0,x2,min,ff);
      }
      
      double acc=sqrt(pow(x1-0.5,2.0)+pow(min+1.0,2.0));
      if (acc<1.0e-2 && acc>1.0e-16) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_mcn[ix].add(t1);
	evals_mcn[ix].add(((double)(cnt-1)));
	if (!has_data2[ix]) has_data2[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "MBB   ";

  min_brent_boost<> mbb;

  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    mbb.tol_abs=pow(10.0,-6-((double)j)/18.0);
    mbb.tol_rel=pow(10.0,-6-((double)j)/18.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();

      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=1.9+gr.random()*0.01;
	x2=2.0+gr.random()*0.01;

	mbb.min_bkt(x1,x0,x2,min,f2);
      }
      
      double acc=sqrt(pow(x1-0.5,2.0)+pow(min+1.0,2.0));
      if (acc>1.0e-16 && acc<1.0e-2) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_mbb[ix].add(t1);
	evals_mbb[ix].add(((double)(cnt-1)));
	if (!has_data3[ix]) has_data3[ix]=true;
      }
    }
  }
  cout << endl;

  cout << endl;

  cout << "Class min_brent_gsl" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (has_data1[i]) {
      double avg, sd, avge;
      time_mbg[i].current_avg(avg,sd,avge);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_mbg[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << endl;
    }
  }
  cout << endl;
  
  cout << "Class min_cern" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (has_data2[i]) {
      double avg, sd, avge;
      time_mcn[i].current_avg(avg,sd,avge);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_mcn[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << endl;
    }
  }
  cout << endl;

  cout << "Class min_brent_boost" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (has_data3[i]) {
      double avg, sd, avge;
      time_mbb[i].current_avg(avg,sd,avge);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_mbb[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << endl;
    }
  }
  cout << endl;
  
  t.report();
  return 0;
}
