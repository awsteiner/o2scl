/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
// sphinx-example-start
/* Example: ex_interp.cpp
   -------------------------------------------------------------------
   A simple example for interpolation
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/interp.h>
#include <o2scl/interp_krige.h>
#include <o2scl/test_mgr.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

// A function for filling the data and comparing results
double f(double x, const double &mean, const double &sd) {
  return (sin(1.0/(0.3+x))-mean)/sd;
}

double covar(double x1, double x2) {
  //return 1.0*exp(-(x1-x2)*(x1-x2)/0.15/0.15);
  //return 1.0*exp(-(x1-x2)*(x1-x2)/0.13455555/0.1345555);
  return exp(-(x1-x2)*(x1-x2)/0.131226/0.131226/2.0);
}

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  // Create the sample data

  static const size_t N=20;
  ubvector x(N), y(N);
  x[0]=0.0;
  y[0]=f(x[0],0.0,1.0);
  for(size_t i=1;i<N;i++) {
    x[i]=x[i-1]+pow(((double)i)/40.0,2.0);
    y[i]=f(x[i],0.0,1.0);
  }
  
  table<> tdata;
  tdata.line_of_names("x y");

  double y_mean=vector_mean(y);
  for(size_t i=0;i<N;i++) {
    y[i]-=y_mean;
  }
  double y_sd=vector_stddev(y);
  cout << y_mean << " " << y_sd << endl;
  for(size_t i=0;i<N;i++) {
    y[i]/=y_sd;
    double line[2]={x[i],y[i]};
    tdata.line_of_data(2,line);
  }
  double y_mean2=vector_mean(y);
  double y_sd2=vector_stddev(y);
  cout << y_mean2 << " " << y_sd2 << endl;

  interp_vec<ubvector> iv_lin(N,x,y,itp_linear);
  interp_vec<ubvector> iv_csp(N,x,y,itp_cspline);
  interp_vec<ubvector> iv_aki(N,x,y,itp_akima);
  interp_vec<ubvector> iv_mon(N,x,y,itp_monotonic);
  interp_vec<ubvector> iv_stef(N,x,y,itp_steffen);

  interp_krige_optim<ubvector> iko;
  iko.verbose=2;
  iko.nlen=100;
  iko.mode=interp_krige_optim<ubvector>::mode_loo_cv;
  iko.set_noise(N,x,y,1.0e-10);
  cout << endl;

  interp_krige_optim<ubvector> iko2;
  iko2.verbose=2;
  iko2.nlen=100;
  iko2.mode=interp_krige_optim<ubvector>::mode_max_lml;
  iko2.set_noise(N,x,y,1.0e-10);
  cout << endl;

  interp_krige<ubvector> ik;
  std::function<double(double,double)> fc=covar;
  ik.set_covar_noise(N,x,y,fc,1.0e-10);

  double max=x[x.size()-1];

  cout << "\nx              y             iko           iko2: " << endl;
  for(size_t i=0;i<N;i++) {
    cout.setf(ios::showpos);
    cout << x[i] << " " << f(x[i],y_mean,y_sd) << " "
	 << iko.eval(x[i]) << " "
	 << iko2.eval(x[i]) << endl;
    cout.unsetf(ios::showpos);
  }
  cout << endl;

  cout << "\nx              y             iko           iko2: " << endl;
  size_t N2=N*100;
  table<> tresult;
  tresult.line_of_names("x y ylin ycsp yaki ymon ystef yiko yiko_lml"); 
  for(size_t i=0;i<=N2;i++) {
    double x=((double)i)/((double)N2)*max;
    double line[9]={x,f(x,y_mean,y_sd),iv_lin.eval(x),iv_csp.eval(x),
		    iv_aki.eval(x),iv_mon.eval(x),iv_stef.eval(x),
		    iko.eval(x),iko2.eval(x)};
    tresult.line_of_data(9,line);
    if (i%50==0) {
      cout.setf(ios::showpos);
      cout << x << " " << f(x,y_mean,y_sd) << " " << iko.eval(x) << " "
	   << iko2.eval(x) << endl;
      cout.unsetf(ios::showpos);
    }
  }

  hdf_file hf;
  hf.open_or_create("ex_interp.o2");
  hdf_output(hf,tdata,"tdata");
  hdf_output(hf,tresult,"tresult");
  hf.close();

  
  t.report();
  return 0;
}
// End of example
