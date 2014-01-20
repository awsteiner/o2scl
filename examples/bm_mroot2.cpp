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
  This program compares the accuracy and speed of gsl_mroot_hybrids
  with cern_mroot. The GSL solver is about 10 to 20 percent faster for
  this example and has the advantage of being able to avoid
  singularities. However, the CERN solver doesn't require all
  equations to be computed for every function evaulation, and may thus
  be much faster in other selected problems.
*/

#include <cmath>
#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/ovector_tlate.h>
#include <o2scl/gsl_mroot_hybrids.h>
#include <o2scl/cern_mroot.h>
#include <o2scl/gsl_rnga.h>
#include <o2scl/hist_ev.h>

using namespace std;
using namespace o2scl;

int gfn(size_t nv, const ovector_base &x, ovector_base &y, double &pa) {
  y[0]=sin(x[1]-pa);
  y[1]=sin(x[0]-pa*2.0)+exp(x[0]-pa);
  return 0;
}

int main(void) {

  const size_t N=100;

  ovector x(2), y(2);
  gsl_rnga gr;

  // Store timing data in a hist_ev object
  hist_ev hev_gmh, hev_cmr;
  // Choose a logarithmic grid with 10 blocks
  hist_new hist;
  hev_gmh.set_grid_blocks(hist.grid_end_size(1.0e-16,1.0e-4,10.0,true),10);
  hev_cmr.set_grid_blocks(hist.grid_end_size(1.0e-16,1.0e-4,10.0,true),10);

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // Time gsl_mroot_hybrids
  gsl_mroot_hybrids<> mroot_gmh;
  cout << "GSL  ";
  for(size_t j=0;j<60;j++) {
    cout << "." << flush;
    
    mroot_gmh.tol_rel=pow(10.0,-4-((double)j)/9.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      double pa=0.0;

      for(size_t ell=0;ell<N;ell++) {
	pa=gr.random()/2.0;
	mm_funct_fptr_param<double> mff(gfn,pa);
	x[0]=0.5;
	x[1]=0.5;
	
	mroot_gmh.msolve(2,x,mff);
      }
      
      gfn(2,x,y,pa);
      double acc=fabs(y[0])+fabs(y[1]);
      if (acc<fabs(1.0e-16)) acc=1.0001e-16;
      
      double t1=(clock()-tmp)/10000.0;
      hev_gmh.add(acc,t1);
    }
  }
  cout << endl;
  
  // Time cern_mroot
  cern_mroot<> mroot_cmr;
  cout << "CERN ";
  for(size_t j=0;j<60;j++) {

    cout << "." << flush;
    
    mroot_cmr.tol_abs=pow(10.0,-2-((double)j)/18.0);
    mroot_cmr.tol_rel=pow(10.0,-2-((double)j)/18.0);
    
    for(size_t k=0;k<N;k++) {
      size_t tmp=clock();
      double pa=0.0;

      for(size_t ell=0;ell<N;ell++) {
	pa=gr.random()/2.0;
	mm_funct_fptr_param<double> mff(gfn,pa);
	x[0]=0.5;
	x[1]=0.5;
	
	mroot_cmr.msolve(2,x,mff);
      }
      gfn(2,x,y,pa);
      double acc=fabs(y[0])+fabs(y[1]);
      if (acc<fabs(1.0e-16)) acc=1.0001e-16;
      
      double t1=(clock()-tmp)/10000.0;
      hev_cmr.add(acc,t1);
    }
  }
  cout << endl;
  cout << endl;

  cout << "GSL, gsl_mroot_hybrids" << endl;
  uvector rep, avg, sd, avge;
  hev_gmh.current_avg(rep,avg,sd,avge);
  cout << "Accuracy     Avg. time    Error in avg. time" << endl;
  for(size_t i=0;i<rep.size();i++) {
    cout << rep[i] << " " << avg[i] << " " << avge[i] << endl;
  }
  cout << endl;
  
  cout << "CERNLIB, cern_mroot" << endl;
  hev_cmr.current_avg(rep,avg,sd,avge);
  cout << "Accuracy     Avg. time    Error in avg. time" << endl;
  for(size_t i=0;i<rep.size();i++) {
    cout << rep[i] << " " << avg[i] << " " << avge[i] << endl;
  }
  cout << endl;

  t.report();
  return 0;
}
