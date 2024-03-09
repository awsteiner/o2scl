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
// sphinx-example-start
/* Example: ex_ode.cpp
   -------------------------------------------------------------------
   An example to demonstrate solving differential equations. See "License 
   Information" section of the documentation for license information.
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_gamma.h>
#include <o2scl/test_mgr.h>
#include <o2scl/ode_funct.h>
#include <o2scl/ode_rkck_gsl.h>
#include <o2scl/ode_rk8pd_gsl.h>
#include <o2scl/astep_gsl.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;

// Differential equation defining the Bessel function. This assumes
// the second derivative at x=0 is 0 and thus only works for odd alpha.
int derivs(double x, size_t nv, const ubvector &y, 
	   ubvector &dydx, double &alpha) {
  dydx[0]=y[1];
  if (x==0.0) dydx[1]=0.0;
  else dydx[1]=(-x*y[1]+(-x*x+alpha*alpha)*y[0])/x/x;
  return 0;
}

// Differential equation defining the Airy function, Ai(x)
int derivs2(double x, size_t nv, const ubvector &y, 
	    ubvector &dydx) {
  dydx[0]=y[1];
  dydx[1]=y[0]*x;
  return 0;
}

// Differential equation defining the Bessel function using
// a vector object for o2scl::ode_iv_solve_grid
int derivs3(double x, size_t nv, const o2scl::solve_grid_mat_row &y, 
	    o2scl::solve_grid_mat_row &dydx, double &alpha) {
  dydx[0]=y[1];
  if (x==0.0) dydx[1]=0.0;
  else dydx[1]=(-x*y[1]+(-x*x+alpha*alpha)*y[0])/x/x;
  return 0;
}

int main(void) {

  cout.setf(ios::scientific);
  cout.setf(ios::showpos);

  // The independent variable and stepsize
  double x, dx=1.0e-1;

  // The function and derivative values and the estimated errors
  ubvector y(2), dydx(2), yout(2), yerr(2), dydx_out(2);

  test_mgr t;
  t.set_output_level(1);

  // The parameter for the Bessel function
  double alpha=1.0;
  
  // Specify the differential equations to solve
  ode_funct od=
    std::bind(derivs,std::placeholders::_1,std::placeholders::_2,
	      std::placeholders::_3,std::placeholders::_4,alpha);
  ode_funct od2=derivs2;
  ode_funct_solve_grid od3=
    std::bind(derivs3,std::placeholders::_1,std::placeholders::_2,
	      std::placeholders::_3,std::placeholders::_4,alpha);
  
  // The basic ODE steppers
  ode_rkck_gsl<> ode;
  ode_rk8pd_gsl<> ode2;

  // ------------------------------------------------------------
  // Store the results in tables

  table<> tab[8];
  tab[0].line_of_names("x calc exact diff err");
  tab[1].line_of_names("x calc exact diff err");
  tab[2].line_of_names("x calc exact diff err");
  tab[3].line_of_names("x calc exact diff err");
  tab[4].line_of_names("x calc exact diff err0 err1");
  tab[5].line_of_names("x calc exact diff err0 err1");
  tab[6].line_of_names("x calc exact diff err0 err1");
  tab[7].line_of_names("x calc exact diff err0");

  // ------------------------------------------------------------
  // Solution 1: Solve using the non-adaptive Cash-Karp stepper.

  cout << "Bessel function, Cash-Karp: " << endl;

  // Initial values at x=0
  x=0.0;
  y[0]=0.0;
  y[1]=0.5;

  // The non-adaptive ODE steppers require the derivatives as
  // input
  derivs(x,2,y,dydx,alpha);

  cout << " x             J1(calc)      J1(exact)     rel. diff.    "
       << "err" << endl;

  while (x<1.0) {

    // Perform a step. Since the fourth and sixth arguments are 
    // the same, the values in 'y' are updated with the new values
    // at x+dx. 
    ode.step(x,dx,2,y,dydx,y,yerr,dydx,od);

    // Update the x value
    x+=dx;

    // Print and test
    cout << x << " " << y[0] << " " 
	 << gsl_sf_bessel_J1(x) << " ";
    cout << fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)) << " ";
    cout << yerr[0] << endl;
    t.test_rel(y[0],gsl_sf_bessel_J1(x),5.0e-5,"rkck");

    // Also output the results to a table
    vector<double> line={x,y[0],gsl_sf_bessel_J1(x),
      fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)),
      yerr[0]};
    tab[0].line_of_data(line.size(),line);
  }

  // Compare with the exact result at the last point
  cout << "Accuracy at end: " 
       << fabs(y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x) << endl;
  cout << endl;

  // End of Solution 1
  // ------------------------------------------------------------
  // Solution 2: Solve using the non-adaptive Prince-Dormand stepper.
  // Note that for the Bessel function, the 8th order stepper performs
  // worse than the 4th order. The error returned by the stepper is
  // larger near x=0, as expected.
  
  cout << "Bessel function, Prince-Dormand: " << endl;
  x=0.0;
  y[0]=0.0;
  y[1]=0.5;
  derivs(x,2,y,dydx,alpha);
  cout << " x             J1(calc)      J1(exact)     rel. diff.    "
       << "err" << endl;
  while (x<1.0) {
    ode2.step(x,dx,2,y,dydx,y,yerr,dydx,od);
    x+=dx;
    cout << x << " " << y[0] << " " 
	 << gsl_sf_bessel_J1(x) << " ";
    cout << fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)) << " ";
    cout << yerr[0] << endl;
    t.test_rel(y[0],gsl_sf_bessel_J1(x),5.0e-4,"rk8pd");

    // Also output the results to a table
    vector<double> line={x,y[0],gsl_sf_bessel_J1(x),
      fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)),
      yerr[0]};
    tab[1].line_of_data(line.size(),line);
  }
  cout << "Accuracy at end: " 
       << fabs(y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x) << endl;
  cout << endl;

  // End of Solution 2
  // ------------------------------------------------------------
  // Solution 3: Solve using the non-adaptive Cash-Karp stepper.

  cout << "Airy function, Cash-Karp: " << endl;
  x=0.0;
  y[0]=1.0/pow(3.0,2.0/3.0)/gsl_sf_gamma(2.0/3.0);
  y[1]=-1.0/pow(3.0,1.0/3.0)/gsl_sf_gamma(1.0/3.0);
  derivs2(x,2,y,dydx);
  cout << " x             Ai(calc)      Ai(exact)     rel. diff.    "
       << "err" << endl;
  while (x<1.0) {
    ode.step(x,dx,2,y,dydx,y,yerr,dydx,od2);
    x+=dx;
    cout << x << " " << y[0] << " " 
	 << gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE) << " ";
    cout << fabs((y[0]-gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE))/
		 gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE)) << " ";
    cout << yerr[0] << endl;
    t.test_rel(y[0],gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE),1.0e-8,"rkck");

    // Also output the results to a table
    vector<double> line={x,y[0],gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE),
      fabs((y[0]-gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE))/
           gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE)),
      yerr[0]};
    tab[2].line_of_data(line.size(),line);
  }
  cout << "Accuracy at end: " 
       << fabs(y[0]-gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE))/
    gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE) << endl;
  cout << endl;
  
  // End of Solution 3
  // ------------------------------------------------------------
  // Solution 4: Solve using the non-adaptive Prince-Dormand stepper.
  // On this function, the higher-order routine performs significantly
  // better.
  
  cout << "Airy function, Prince-Dormand: " << endl;
  x=0.0;
  y[0]=1.0/pow(3.0,2.0/3.0)/gsl_sf_gamma(2.0/3.0);
  y[1]=-1.0/pow(3.0,1.0/3.0)/gsl_sf_gamma(1.0/3.0);
  derivs2(x,2,y,dydx);
  cout << " x             Ai(calc)      Ai(exact)     rel. diff.    "
       << "err" << endl;
  while (x<1.0) {
    ode2.step(x,dx,2,y,dydx,y,yerr,dydx,od2);
    x+=dx;
    cout << x << " " << y[0] << " " 
	 << gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE) << " ";
    cout << fabs((y[0]-gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE))/
		 gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE)) << " ";
    cout << yerr[0] << endl;
    t.test_rel(y[0],gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE),1.0e-14,"rk8pd");

    // Also output the results to a table
    vector<double> line={x,y[0],gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE),
      fabs((y[0]-gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE))/
           gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE)),
      yerr[0]};
    tab[3].line_of_data(line.size(),line);
  }
  cout << "Accuracy at end: " 
       << fabs(y[0]-gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE))/
    gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE) << endl;
  cout << endl;

  // End of Solution 4
  // ------------------------------------------------------------
  // Solution 5: Solve using the GSL adaptive stepper

  // Lower the output precision to fit in 80 columns
  cout.precision(5);

  cout << "Adaptive stepper: " << endl;
  astep_gsl<> ode3;
  x=0.0;
  y[0]=0.0;
  y[1]=0.5;
  cout << "   x            J1(calc)     J1(exact)    rel. diff.";
  cout << "   err_0        err_1" << endl;
  int k=0;
  while (x<10.0) {
    int retx=ode3.astep(x,10.0,dx,2,y,dydx,yerr,od);
    if (k%3==0) {
      cout << retx << " " << x << " " << y[0] << " " 
	   << gsl_sf_bessel_J1(x) << " ";
      cout << fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)) << " ";
      cout << yerr[0] << " " << yerr[1] << endl;
    }
    t.test_rel(y[0],gsl_sf_bessel_J1(x),5.0e-3,"astep");
    t.test_rel(y[1],0.5*(gsl_sf_bessel_J0(x)-gsl_sf_bessel_Jn(2,x)),
               5.0e-3,"astep 2");
    t.test_rel(yerr[0],0.0,4.0e-6,"astep 3");
    t.test_rel(yerr[1],0.0,4.0e-6,"astep 4");
    t.test_gen(retx==0,"astep 5");

    // Also output the results to a table
    vector<double> line={x,y[0],gsl_sf_bessel_J1(x),
      fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)),
      yerr[0],yerr[1]};
    tab[4].line_of_data(line.size(),line);

    k++;
  }
  cout << "Accuracy at end: "
       << fabs(y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x) << endl;
  cout << endl;

  // End of Solution 5
  // ------------------------------------------------------------
  // Solution 6: Solve using the GSL adaptive stepper.
  // Decrease the tolerances, and the adaptive stepper takes
  // smaller step sizes.

  cout << "Adaptive stepper with smaller tolerances: " << endl;
  ode3.con.eps_abs=1.0e-8;
  ode3.con.a_dydt=1.0;
  x=0.0;
  y[0]=0.0;
  y[1]=0.5;
  cout << "   x            J1(calc)     J1(exact)    rel. diff.";
  cout << "   err_0        err_1" << endl;
  k=0;
  while (x<10.0) {
    int retx=ode3.astep(x,10.0,dx,2,y,dydx,yerr,od);
    if (k%3==0) {
      cout << retx << " " << x << " " << y[0] << " " 
	   << gsl_sf_bessel_J1(x) << " ";
      cout << fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)) << " ";
      cout << yerr[0] << " " << yerr[1] << endl;
    }
    t.test_rel(y[0],gsl_sf_bessel_J1(x),5.0e-3,"astep");
    t.test_rel(y[1],0.5*(gsl_sf_bessel_J0(x)-gsl_sf_bessel_Jn(2,x)),
               5.0e-3,"astep 2");
    t.test_rel(yerr[0],0.0,4.0e-8,"astep 3");
    t.test_rel(yerr[1],0.0,4.0e-8,"astep 4");
    t.test_gen(retx==0,"astep 5");

    // Also output the results to a table
    vector<double> line={x,y[0],gsl_sf_bessel_J1(x),
      fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)),
      yerr[0],yerr[1]};
    tab[5].line_of_data(line.size(),line);

    k++;
  }
  cout << "Accuracy at end: "
       << fabs(y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x) << endl;
  cout << endl;

  // End of Solution 6
  // ------------------------------------------------------------
  // Solution 7: Solve using the GSL adaptive stepper.
  // Use the higher-order stepper, and less steps are required. The
  // stepper automatically takes more steps near x=0 in order since
  // the higher-order routine has more trouble there.

  cout << "Adaptive stepper, Prince-Dormand: " << endl;
  ode3.set_step(ode2);
  x=0.0;
  y[0]=0.0;
  y[1]=0.5;
  cout << "   x            J1(calc)     J1(exact)    rel. diff.";
  cout << "   err_0        err_1" << endl;
  k=0;
  while (x<10.0) {
    int retx=ode3.astep(x,10.0,dx,2,y,dydx,yerr,od);
    if (k%3==0) {
      cout << retx << " " << x << " " << y[0] << " " 
	   << gsl_sf_bessel_J1(x) << " ";
      cout << fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)) << " ";
      cout << yerr[0] << " " << yerr[1] << endl;
    }
    t.test_rel(y[0],gsl_sf_bessel_J1(x),5.0e-3,"astep");
    t.test_rel(y[1],0.5*(gsl_sf_bessel_J0(x)-gsl_sf_bessel_Jn(2,x)),
               5.0e-3,"astep");
    t.test_rel(yerr[0],0.0,4.0e-8,"astep 3");
    t.test_rel(yerr[1],0.0,4.0e-8,"astep 4");
    t.test_gen(retx==0,"astep 5");

    // Also output the results to a table
    vector<double> line={x,y[0],gsl_sf_bessel_J1(x),
      fabs((y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x)),
      yerr[0],yerr[1]};
    tab[6].line_of_data(line.size(),line);

    k++;
  }
  cout << "Accuracy at end: "
       << fabs(y[0]-gsl_sf_bessel_J1(x))/gsl_sf_bessel_J1(x) << endl;
  cout << endl;
  
  // End of Solution 7
  // ------------------------------------------------------------
  // Solution 8: Solve using the O2scl initial value solver

  // Return the output precision to the default
  cout.precision(6);
  
  cout << "Initial value solver: " << endl;
  ode_iv_solve_grid<> ode4;
  ode4.ntrial=10000;

  // Define the grid and the storage for the solution
  const size_t ngrid=101;
  ubvector xg(ngrid), yinit(2);
  ubmatrix yg(ngrid,2), ypg(ngrid,2), yerrg(ngrid,2);
  for(size_t i=0;i<ngrid;i++) xg[i]=((double)i)/10.0;

  // Set the initial value 
  xg[0]=0.0;
  xg[ngrid-1]=10.0;
  yg(0,0)=0.0;
  yg(0,1)=0.5;

  // Perform the full solution
  ode4.solve_grid<ubvector,ubmatrix>(0.1,2,ngrid,xg,yg,yerrg,ypg,od3);
  
  // Output and test the results
  cout << " x             J1(calc)      J1(exact)     rel. diff." << endl;
  for(size_t i=1;i<ngrid;i++) {

    if (i%10==0) {
      cout << xg[i] << " " << yg(i,0) << " "
	   << gsl_sf_bessel_J1(xg[i]) << " ";
      cout << fabs(yg(i,0)-gsl_sf_bessel_J1(xg[i])) << " "
	   << fabs(yerrg(i,0)) << endl;
    }
    t.test_rel(yg(i,0),gsl_sf_bessel_J1(xg[i]),5.0e-7,"astep");
    t.test_rel(yg(i,1),0.5*(gsl_sf_bessel_J0(xg[i])-
                             gsl_sf_bessel_Jn(2,xg[i])),5.0e-7,"astep 2");

    // Also output the results to a table
    vector<double> line={xg[i],yg(i,0),gsl_sf_bessel_J1(xg[i]),
      fabs(yg(i,0)-gsl_sf_bessel_J1(xg[i])),
      fabs(yerrg(i,0))};
    tab[7].line_of_data(line.size(),line);
  }

  cout << "Accuracy at end: "
       << fabs(yg(ngrid-1,0)-gsl_sf_bessel_J1(xg[ngrid-1]))/
    gsl_sf_bessel_J1(xg[ngrid-1]) << endl;
  cout << endl;

  // End of Solution 8
  // ------------------------------------------------------------
  // Output results to a file

  hdf_file hf;
  hf.open_or_create("ex_ode.o2");
  for(size_t i=0;i<8;i++) {
    hdf_output(hf,tab[i],((string)"table_")+itos(i));
  }
  hf.close();

  // ------------------------------------------------------------

  cout.unsetf(ios::showpos);
  t.report();

  return 0;
}
