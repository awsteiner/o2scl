/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2022, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>

#include <o2scl/fit_min.h>
#include <o2scl/fit_bayes.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

double func2(size_t np, const ubvector &p, double x) {
  return p[0]+x*p[1];
}

double func3(size_t np, const ubvector &p, double x) {
  return p[0]+x*p[1]+x*x*p[2];
}

int main(void) {
  cout.setf(ios::scientific);
  cout.precision(8);

  test_mgr t;
  t.set_output_level(2);

  fit_min<> mf;
  fit_bayes<> bf;
  
  fit_funct_fptr<> ff2(func2);
  fit_funct_fptr<> ff3(func3);

  ubvector xdat(4), ydat(4), sigma(4);
  xdat[0]=0.0;
  xdat[1]=1.0;
  xdat[2]=2.0;
  xdat[3]=3.0;
  ydat[0]=1.1;
  ydat[1]=3.07;
  ydat[2]=5.1;
  ydat[3]=8.0;
  sigma[0]=0.01;
  sigma[1]=0.08;
  sigma[2]=0.15;
  sigma[3]=0.5;

  // ------------------------------------------------------
  // Two parameter frequentist fit

  {
    ubvector x_init2(2), xerr2(2);
    ubmatrix mycovar2(2,2);
    double chi2;

    x_init2[0]=0.0;
    x_init2[1]=1.0;
    chi_fit_funct<> cff2(4,xdat,ydat,sigma,ff2);
    
    mf.fit(2,x_init2,mycovar2,chi2,cff2);
    
    cout << "Parameters: " << x_init2[0] << " " << x_init2[1] << endl;
    cout << "Errors: " 
	 << sqrt(mycovar2(0,0)) << " " << sqrt(mycovar2(1,1)) << endl;
    cout << "chi^2/dof: " << chi2/2.0 << endl;
    cout << "Covariance matrix: " << endl;
    matrix_out(cout,mycovar2);
    cout << endl;
  }

  // ------------------------------------------------------
  // Three parameter frequentist fit

  {
    ubvector x_init3(3), xerr3(3);
    ubmatrix mycovar3(3,3);
    double chi2;

    x_init3[0]=1.0;
    x_init3[1]=1.0;
    x_init3[2]=1.0;
    chi_fit_funct<> cff3(4,xdat,ydat,sigma,ff3);
    
    mf.fit(3,x_init3,mycovar3,chi2,cff3);
    
    cout << "Parameters: " << x_init3[0] << " " << x_init3[1] << " " 
	 << x_init3[2] << endl;
    cout << "Errors: " << sqrt(mycovar3(0,0)) << " " 
	 << sqrt(mycovar3(1,1)) << " " << sqrt(mycovar3(2,2)) << endl;
    cout << "Chi^2/dof: " << chi2 << endl;
    cout << "Covariance matrix: " << endl;
    matrix_out(cout,mycovar3);
    cout << endl;
  }

  // ------------------------------------------------------
  // Two parameter Bayesian fit

  {
    ubvector par_lo(2), par_hi(2), err_lo(2), err_hi(2);
    ubvector par_max(2), err_max(2);
    uniform_prior<> pri;

    par_lo[0]=1.0;
    par_hi[0]=1.2;
    par_lo[1]=1.7;
    par_hi[1]=2.1;
    
    bf.fit(4,xdat,ydat,sigma,2,par_lo,par_max,par_hi,
	   err_lo,err_max,err_hi,ff2,pri);
    
    cout << "Parameters: ";
    vector_out(cout,par_lo,true);
    cout << "            "; 
    vector_out(cout,par_max,true);
    cout << "            "; 
    vector_out(cout,par_hi,true);
    cout << "Errors:     "; 
    vector_out(cout,err_lo,true);
    cout << "            "; 
    vector_out(cout,err_max,true);
    cout << "            "; 
    vector_out(cout,err_hi,true);
    
    double evi, err;
    bf.evidence(4,xdat,ydat,sigma,2,par_lo,par_hi,pri,evi,err);
    cout << "Evidence: " << evi << " " << err << endl;
    cout << endl;
  }
  
  // ------------------------------------------------------
  // Three parameter Bayesian fit

  {
    ubvector par_lo(3), par_hi(3), err_lo(3), err_hi(3);
    ubvector par_max(3), err_max(3);
    uniform_prior<> pri;

    par_lo[0]=1.04;
    par_hi[0]=1.16;
    par_lo[1]=1.0;
    par_hi[1]=2.5;
    par_lo[2]=-0.3;
    par_hi[2]=0.45;
    
    bf.fit(4,xdat,ydat,sigma,3,par_lo,par_max,par_hi,
	   err_lo,err_max,err_hi,ff3,pri);
    
    cout << "Parameters: ";
    vector_out(cout,par_lo,true);
    cout << "            "; 
    vector_out(cout,par_max,true);
    cout << "            "; 
    vector_out(cout,par_hi,true);
    cout << "Errors:     "; 
    vector_out(cout,err_lo,true);
    cout << "            "; 
    vector_out(cout,err_max,true);
    cout << "            "; 
    vector_out(cout,err_hi,true);
    
    double evi, err;
    bf.evidence(4,xdat,ydat,sigma,3,par_lo,par_hi,pri,evi,err);
    cout << "Evidence: " << evi << " " << err << endl;
    cout << endl;
  }
  
  // ------------------------------------------------------
  // Three parameter Bayesian fit with histograms

  {
    vector<hist> par_hist;
    ubvector par_lo(3), par_hi(3);
    uniform_prior<> pri;

    par_lo[0]=1.04;
    par_hi[0]=1.16;
    par_lo[1]=1.0;
    par_hi[1]=2.5;
    par_lo[2]=-0.3;
    par_hi[2]=0.45;
    
    bf.fit_hist(4,xdat,ydat,sigma,3,par_lo,par_hi,par_hist,ff3,pri);

    cout.precision(5);
    cout << "Parameter histograms: " << endl;
    for(size_t j=0;j<par_hist[0].size();j++) {
      for(size_t i=0;i<3;i++) {
	if (i==2) {
	  cout.setf(ios::showpos);
	  cout << par_hist[i].get_rep_i(j) << " ";
	  cout.unsetf(ios::showpos);
	} else {
	  cout << par_hist[i].get_rep_i(j) << " ";
	}
	cout << par_hist[i].get_wgt_i(j) << " ";
      }
      cout << endl;
    }
    cout << endl;
  }

  t.report();
  return 0;

}
