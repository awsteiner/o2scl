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
/* Example: ex_stiff.cpp
   -------------------------------------------------------------------
   An example to demonstrate solving stiff differential equations
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/ode_funct.h>
#include <o2scl/astep_gsl.h>
#include <o2scl/ode_bsimp_gsl.h>
#include <o2scl/table.h>

#ifdef O2SCL_HDF
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#endif

using namespace std;
using namespace o2scl;
#ifdef O2SCL_HDF
using namespace o2scl_hdf;
#endif

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int derivs(double x, size_t nv, const ubvector &y, 
	   ubvector &dydx) {
  dydx[0]=480.0*y[0]+980.0*y[1];
  dydx[1]=-490.0*y[0]-990.0*y[1];
  return 0;
}

int jac(double x, size_t nv, const ubvector &y, 
	ubmatrix &dfdy, ubvector &dfdx) {
  dfdy(0,0)=480.0;
  dfdy(0,1)=980.0;
  dfdy(1,0)=-490.0;
  dfdy(1,1)=-990.0;
  dfdx[0]=0.0;
  dfdx[1]=0.0;
  return 0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);
  cout.precision(3);

  // Specification of the differential equations and the Jacobian
  ode_funct od11=derivs;
  ode_jac_funct oj=jac;
  
  table<> tab[2];
  tab[0].line_of_names("x calc exact rel_err rel_diff");
  tab[1].line_of_names("x calc exact rel_err rel_diff");
  
  // ------------------------------------------------------------
  // First solve with ode_bsimp_gsl, designed to handle stiff ODEs

  ode_bsimp_gsl<> gb;

  double x1, dx=1.0e-1;
  ubvector y1(2), dydx1(2), yout1(2), yerr1(2), dydx_out1(2);

  x1=0.0;
  y1[0]=1.0;
  y1[1]=0.0;

  derivs(x1,2,y1,dydx1);

  for(size_t i=1;i<=40;i++) {
    
    gb.step(x1,dx,2,y1,dydx1,y1,yerr1,dydx1,od11,oj);
    x1+=dx;

    double exact[2]={-exp(-500.0*x1)+2.0*exp(-10.0*x1),
		     exp(-500.0*x1)-exp(-10.0*x1)};
    
    cout.setf(ios::showpos);
    cout << x1 << " " << y1[0] << " " << y1[1] << " " 
	 << yerr1[0] << " " << yerr1[1] << " " << exact[0] << " "
	 << exact[1] << endl;
    cout.unsetf(ios::showpos);
    
    double line[5]={x1,y1[0],exact[0],fabs(yerr1[0]/exact[0]),
		    fabs((y1[0]-exact[0])/exact[0])};
    tab[0].line_of_data(5,line);

    t.test_rel(y1[0],exact[0],3.0e-4,"y0");
    t.test_rel(y1[1],exact[1],3.0e-4,"y1");

  }

  cout << endl;

  // ------------------------------------------------------------
  // Now compare to the traditional adaptive stepper

  astep_gsl<> ga;

  double x2;
  ubvector y2(2), dydx2(2), yout2(2), yerr2(2), dydx_out2(2);

  x2=0.0;
  y2[0]=1.0;
  y2[1]=0.0;

  derivs(x2,2,y2,dydx2);

  size_t j=0;
  while (x2<4.0) {

    ga.astep(x2,4.0,dx,2,y2,dydx2,yerr2,od11);

    double exact[2]={-exp(-500.0*x2)+2.0*exp(-10.0*x2),
		     exp(-500.0*x2)-exp(-10.0*x2)};
    
    if (j%25==0) {
      cout.setf(ios::showpos);
      cout << x2 << " " << y2[0] << " " << y2[1] << " " 
	   << yerr2[0] << " " << yerr2[1] << " " << exact[0] << " "
	   << exact[1] << endl;
      cout.unsetf(ios::showpos);
    }
    j++;
    
    double line[5]={x2,y2[0],exact[0],fabs(yerr2[0]/exact[0]),
		    fabs((y2[0]-exact[0])/exact[0])};
    tab[1].line_of_data(5,line);
    
  }

  cout << endl;

  // ------------------------------------------------------------
  // Output results to a file

#ifdef O2SCL_HDF
  hdf_file hf;
  hf.open_or_create("ex_stiff.o2");
  for(size_t i=0;i<2;i++) {
    hdf_output(hf,tab[i],((string)"table_")+itos(i));
  }
  hf.close();
#endif

  // ------------------------------------------------------------

  t.report();

  return 0;

}
// End of example
