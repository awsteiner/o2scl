/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2024, Andrew W. Steiner

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
  Compare O2scl solvers with Boost
*/

#include <cmath>
#include <functional>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/root_bkt_cern.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/root_stef.h>
#include <o2scl/root_toms748.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/hist.h>
#include <o2scl/expval.h>

#include <boost/math/tools/roots.hpp>

using namespace std;
using namespace o2scl;

size_t cnt;

static const double flatp=2.0;

/* A function which challenges a bracketing solver by including
   several near solutions inside the interval
*/
double func(double x) {
  cnt++;
  return atan((x-0.01)*4)*(1.0+sin((x-0.01)*50.0)/1.1);
}

/* A function which challenges a bracketing solver by including
   several near solutions inside the interval and has 
   a solution exactly at x=0
*/
double func3(double x) {
  cnt++;
  return atan((x-0.0)*4)*(1.0+sin((x-0.01)*50.0)/1.1);
}

/* A function designed to be difficult for derivative solvers
   because it is nearly flat except for a small interval 
   near the solution
 */
double func2(double x) {
  cnt++;
  double fermi1=1.0/(1.0+exp(flatp*(x-0.5)));
  double fermi2=1.0/(1.0+exp(flatp*(x+0.5)));
  double fermi3=1.0/(1.0+exp(flatp*(x-0.01)));
  return fermi1+fermi2+fermi3-1.5;
}

/* The derivative of func2() */
double dfunc2(double x) {
  double dfermi1=-flatp*exp(flatp*(x-0.5))/pow(1.0+exp(flatp*(x-0.5)),2.0);
  double dfermi2=-flatp*exp(flatp*(x+0.5))/pow(1.0+exp(flatp*(x+0.5)),2.0);
  double dfermi3=-flatp*exp(flatp*(x-0.01))/pow(1.0+exp(flatp*(x-0.01)),2.0);
  return dfermi1+dfermi2+dfermi3;
}

/* The form of func2() and dfunc2() for Boost */
boost::math::tuple<double,double> bfunc2(double x) {
  cnt++;
  double fermi1=1.0/(1.0+exp(flatp*(x-0.5)));
  double fermi2=1.0/(1.0+exp(flatp*(x+0.5)));
  double fermi3=1.0/(1.0+exp(flatp*(x-0.01)));
  double dfermi1=-flatp*exp(flatp*(x-0.5))/pow(1.0+exp(flatp*(x-0.5)),2.0);
  double dfermi2=-flatp*exp(flatp*(x+0.5))/pow(1.0+exp(flatp*(x+0.5)),2.0);
  double dfermi3=-flatp*exp(flatp*(x-0.01))/pow(1.0+exp(flatp*(x-0.01)),2.0);
  return boost::math::make_tuple(fermi1+fermi2+fermi3-1.5,
				 dfermi1+dfermi2+dfermi3);
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  const size_t N=300;

  rng_gsl gr;

  // Choose a logarithmic grid
  hist h;
  // We have to adjust the bin size a bit to avoid problems with
  // finite precision
  h.set_bin_edges(uniform_grid_log_end_width<>(1.0e-16,1.0e-4,10.0-1.0e-8));
  
  // Store timing data in expect_val objects with 10 blocks
  vector<expval_scalar> time_grb(h.size());
  vector<expval_scalar> time_cr(h.size());
  vector<expval_scalar> time_rt(h.size());
  vector<expval_scalar> time_grb2(h.size());
  vector<expval_scalar> time_cr2(h.size());
  vector<expval_scalar> time_rt2(h.size());
  vector<expval_scalar> time_grs(h.size());
  vector<expval_scalar> time_bnr(h.size());

  vector<expval_scalar> evals_grb(h.size());
  vector<expval_scalar> evals_cr(h.size());
  vector<expval_scalar> evals_rt(h.size());
  vector<expval_scalar> evals_grb2(h.size());
  vector<expval_scalar> evals_cr2(h.size());
  vector<expval_scalar> evals_rt2(h.size());
  vector<expval_scalar> evals_grs(h.size());
  vector<expval_scalar> evals_bnr(h.size());

  for(size_t i=0;i<h.size();i++) {

    time_grb[i].set_blocks(10,1);
    time_cr[i].set_blocks(10,1);
    time_rt[i].set_blocks(10,1);
    time_grb2[i].set_blocks(10,1);
    time_cr2[i].set_blocks(10,1);
    time_rt2[i].set_blocks(10,1);
    time_grs[i].set_blocks(10,1);
    time_bnr[i].set_blocks(10,1);
    
    evals_grb[i].set_blocks(10,1);
    evals_cr[i].set_blocks(10,1);
    evals_rt[i].set_blocks(10,1);
    evals_grb2[i].set_blocks(10,1);
    evals_cr2[i].set_blocks(10,1);
    evals_rt2[i].set_blocks(10,1);
    evals_grs[i].set_blocks(10,1);
    evals_bnr[i].set_blocks(10,1);
  }

  vector<bool> data_grb(h.size());
  vector<bool> data_cr(h.size());
  vector<bool> data_rt(h.size());
  vector<bool> data_grb2(h.size());
  vector<bool> data_cr2(h.size());
  vector<bool> data_rt2(h.size());
  vector<bool> data_grs(h.size());
  vector<bool> data_bnr(h.size());

  funct11 ff=func;
  funct11 ff2=func2;
  funct11 df2=dfunc2;
  funct11 ff3=func3;

  std::function<double (double)> fb=func;
  std::function<double (double)> fb3=func3;

  root_brent_gsl<> grb;
  root_bkt_cern<> cr;
  root_toms748<> rt;
  root_stef<> grs;

  cout << "Collecting data:" << endl;
  cout << "RBG    ";
  grb.test_form=0;
  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    grb.tol_abs=pow(10.0,-4-((double)j)/9.0);
    grb.tol_rel=pow(10.0,-4-((double)j)/9.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      double x0, x1;

      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=2.0+gr.random()*0.01;

	grb.solve_bkt(x0,x1,ff);
      }
      
      double acc=fabs(func(x0));
      if (acc>1.0e-16 && acc<1.0e-4) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_grb[ix].add(t1);
	evals_grb[ix].add(((double)(cnt-1)));
	if (!data_grb[ix]) data_grb[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "RBC    ";
  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    cr.tol_abs=pow(10.0,-4-((double)j)/7.0);
    cr.tol_rel=pow(10.0,-4-((double)j)/7.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      double x0, x1;
      
      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=2.0+gr.random()*0.01;

	cr.solve_bkt(x0,x1,ff);
      }
      
      double acc=fabs(func(x0));

      if (acc<1.0e-4 && acc>1.0e-16) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_cr[ix].add(t1);
	evals_cr[ix].add(((double)(cnt-1)));
	if (!data_cr[ix]) data_cr[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "RT748  ";
  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    rt.tol_abs=pow(10.0,-4-((double)j)/7.0);
    rt.tol_rel=pow(10.0,-4-((double)j)/7.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      double x0, x1;
      
      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=2.0+gr.random()*0.01;

	rt.solve_bkt(x0,x1,fb);
      }
      
      double acc=fabs(func(x0));

      if (acc<1.0e-4 && acc>1.0e-16) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_rt[ix].add(t1);
	evals_rt[ix].add(((double)(cnt-1)));
	if (!data_rt[ix]) data_rt[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "RBG2   ";
  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    grb.tol_abs=pow(10.0,-4-((double)j)/9.0);
    grb.tol_rel=pow(10.0,-4-((double)j)/9.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      double x0, x1;

      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=2.0+gr.random()*0.01;

	grb.solve_bkt(x0,x1,ff3);
      }
      
      double acc=fabs(func3(x0));
      if (acc>1.0e-16 && acc<1.0e-4) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_grb2[ix].add(t1);
	evals_grb2[ix].add(((double)(cnt-1)));
	if (!data_grb2[ix]) data_grb2[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "RBC2   ";
  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    cr.tol_abs=pow(10.0,-4-((double)j)/7.0);
    cr.tol_rel=pow(10.0,-4-((double)j)/7.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      double x0, x1;
      
      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=2.0+gr.random()*0.01;

	cr.solve_bkt(x0,x1,ff3);
      }
      
      double acc=fabs(func3(x0));

      if (acc<1.0e-4 && acc>1.0e-16) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_cr2[ix].add(t1);
	evals_cr2[ix].add(((double)(cnt-1)));
	if (!data_cr2[ix]) data_cr2[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "RT7482 ";
  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    rt.tol_abs=pow(10.0,-4-((double)j)/7.0);
    rt.tol_rel=pow(10.0,-4-((double)j)/7.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      double x0, x1;
      
      for(size_t ell=0;ell<N;ell++) {
	cnt=0;

	// Random numbers near -1 and 2
	x0=-1.0-gr.random()*0.01;
	x1=2.0+gr.random()*0.01;

	rt.solve_bkt(x0,x1,fb3);
      }
      
      double acc=fabs(func3(x0));

      if (acc<1.0e-4 && acc>1.0e-16) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_rt2[ix].add(t1);
	evals_rt2[ix].add(((double)(cnt-1)));
	if (!data_rt2[ix]) data_rt2[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "RS     ";

  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    grs.tol_abs=pow(10.0,-4-((double)j)/9.0);
    grs.tol_rel=pow(10.0,-4-((double)j)/9.0);
    grs.ntrial=1000;
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      
      double x;
      for(size_t ell=0;ell<N;ell++) {
	cnt=0;
	x=0.01+gr.random();
	grs.solve_de(x,ff2,df2);
      }
      
      double acc=fabs(func2(x));

      if (acc>1.0e-16 && acc<1.0e-4) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_grs[ix].add(t1);
	evals_grs[ix].add(((double)(cnt-1)));
	if (!data_grs[ix]) data_grs[ix]=true;
      }
    }
  }
  cout << endl;

  cout << "BNR    ";
  std::function<boost::math::tuple<double,double>(double)> f3=bfunc2;

  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      
      double x;
      for(size_t ell=0;ell<N;ell++) {
	cnt=0;
	x=0.01+gr.random();
	double min=-10.0;
	double max=10.0;
	x=boost::math::tools::newton_raphson_iterate
	  (f3,x,min,max,((double)(j))/2.0);
      }
      
      double acc=fabs(func2(x));

      if (acc>1.0e-16 && acc<1.0e-4) {
	double t1=(clock()-tmp)/10000.0;
	size_t ix=h.get_bin_index(acc);
	time_bnr[ix].add(t1);
	evals_bnr[ix].add(((double)(cnt-1)));
	if (!data_bnr[ix]) data_bnr[ix]=true;
      }
    }
  }
  cout << endl;
  cout << endl;

  cout << "Class root_brent_gsl, function func()" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (data_grb[i]) {
      double avg, sd, avge;
      size_t b1, b2;
      time_grb[i].current_avg_stats(avg,sd,avge,b1,b2);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_grb[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << " ";
      cout << b1 << " " << b2 << endl;
    }
  }
  cout << endl;

  cout << "Class root_bkt_cern, function func()" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (data_cr[i]) {
      double avg, sd, avge;
      size_t b1, b2;
      time_cr[i].current_avg_stats(avg,sd,avge,b1,b2);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_cr[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << " ";
      cout << b1 << " " << b2 << endl;
    }
  }
  cout << endl;
  
  cout << "Class root_toms748, function func()" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (data_rt[i]) {
      double avg, sd, avge;
      size_t b1, b2;
      time_rt[i].current_avg_stats(avg,sd,avge,b1,b2);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_rt[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << " ";
      cout << b1 << " " << b2 << endl;
    }
  }
  cout << endl;

  cout << "Class root_brent_gsl, function func3()" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (data_grb2[i]) {
      double avg, sd, avge;
      size_t b1, b2;
      time_grb2[i].current_avg_stats(avg,sd,avge,b1,b2);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_grb2[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << " ";
      cout << b1 << " " << b2 << endl;
    }
  }
  cout << endl;

  cout << "Class root_bkt_cern, function func3()" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (data_cr2[i]) {
      double avg, sd, avge;
      size_t b1, b2;
      time_cr2[i].current_avg_stats(avg,sd,avge,b1,b2);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_cr2[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << " ";
      cout << b1 << " " << b2 << endl;
    }
  }
  cout << endl;
  
  cout << "Class root_toms748, function func3()" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (data_rt2[i]) {
      double avg, sd, avge;
      size_t b1, b2;
      time_rt2[i].current_avg_stats(avg,sd,avge,b1,b2);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_rt2[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << " ";
      cout << b1 << " " << b2 << endl;
    }
  }
  cout << endl;

  cout << "Class root_stef, function func2()" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (data_grs[i]) {
      double avg, sd, avge;
      time_grs[i].current_avg(avg,sd,avge);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_grs[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << endl;
    }
  }
  cout << endl;

  cout << "Boost Newton-Raphson, function func2()" << endl;
  cout << "Accuracy     Avg. time    Error        Evaluations  Error" << endl;
  for(size_t i=0;i<h.size();i++) {
    if (data_bnr[i]) {
      double avg, sd, avge;
      time_bnr[i].current_avg(avg,sd,avge);
      cout << h.get_rep_i(i) << " " << avg << " " << avge << " ";
      evals_bnr[i].current_avg(avg,sd,avge);
      cout << avg << " " << avge << endl;
    }
  }
  cout << endl;

  t.report();
  return 0;
}
