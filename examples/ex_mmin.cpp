/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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

/* Example: ex_mmin.cpp
   -------------------------------------------------------------------
   Example usage of the multidimensional minimizers with and without
   gradients. 
*/
#include <fstream>
#include <string>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/multi_funct.h>
#include <o2scl/constants.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/mmin_conf.h>
#include <o2scl/mmin_conp.h>
#include <o2scl/mmin_bfgs2.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

class cl {

public:

  cl() {
    param=30.0;
  }

  // To output function evaluations to a file
  ofstream fout;

  // Parameter of the quadratic
  double param;
  
  // Updated spring function
  double spring_two(size_t nv, const ubvector &x) {
    double theta=atan2(x[1],x[0]);
    double r=sqrt(x[0]*x[0]+x[1]*x[1]);
    double z=x[2];
    double tmz=theta-z;
    double fact=8.0-pow(sin(tmz+o2scl_const::pi/2.0)+1.0,3.0);
    double rm1=r-1.0;
    double ret=fact+exp(rm1*rm1)+z*z/param;
    fout << x[0] << " " << x[1] << " " << x[2] << " " << ret << endl;
    return ret;
  }

  // Gradient of the spring function
  int sgrad(size_t nv, ubvector &x, ubvector &g) {

    double theta=atan2(x[1],x[0]);
    double r=sqrt(x[0]*x[0]+x[1]*x[1]);
    double z=x[2];
    double tmz=theta-z;
    double rm1=r-1.0;
    double fact=8.0-pow(sin(tmz+o2scl_const::pi/2.0)+1.0,3.0);

    double dtdx=-x[1]/r/r;
    double dtdy=x[0]/r/r;
    double drdx=x[0]/r;
    double drdy=x[1]/r;
    double dfdt=-3.0*pow(sin(tmz+o2scl_const::pi/2.0)+1.0,2.0)*
      cos(tmz+o2scl_const::pi/2.0);
    double dfdz=2.0*z/param+3.0*pow(sin(tmz+o2scl_const::pi/2.0)+1.0,2.0)*
      cos(tmz+o2scl_const::pi/2.0);
    double dfdr=2.0*rm1*exp(rm1*rm1);

    g[0]=dfdr*drdx+dfdt*dtdx;
    g[1]=dfdr*drdy+dfdt*dtdy;
    g[2]=dfdz;

    return 0;
  }
  
};

int main(void) {
  cl acl;
  ubvector x(3);
  double fmin;
  test_mgr t;

  t.set_output_level(1);
  cout.setf(ios::scientific);

  // Using a member function
#ifdef O2SCL_NO_CPP11
  multi_funct_mfptr<cl> f1(&acl,&cl::spring_two);
  grad_funct_mfptr<cl> f1g(&acl,&cl::sgrad);
#else
  multi_funct11 f1=std::bind
    (std::mem_fn<double(size_t,const ubvector &)>(&cl::spring_two),
     &acl,std::placeholders::_1,std::placeholders::_2);
  grad_funct11 f1g=std::bind
    (std::mem_fn<int(size_t,ubvector &,ubvector &)>(&cl::sgrad),
     &acl,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
#endif

  mmin_simp2<> gm1;
  mmin_conf<> gm2;
  mmin_conp<> gm3;
  mmin_bfgs2<> gm4;

  // This function is difficult to minimize, so more trials
  // are required.
  gm1.ntrial*=10;
  gm2.ntrial*=10;
  gm3.ntrial*=10;
  gm4.ntrial*=10;

  // Simplex minimization
  acl.fout.open("ex_mmin1.dat");
  x[0]=1.0;
  x[1]=1.0;
  x[2]=7.0*o2scl_const::pi;
  gm1.mmin(3,x,fmin,f1);
  acl.fout.close();
  cout << gm1.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,1.0e-4,"1a");
  t.test_rel(x[1],0.0,1.0e-4,"1b");
  t.test_rel(x[2],0.0,1.0e-4,"1c");

  // Fletcher-Reeves conjugate 
  acl.fout.open("ex_mmin2.dat");
  x[0]=1.0;
  x[1]=0.0;
  x[2]=7.0*o2scl_const::pi;
  gm2.mmin(3,x,fmin,f1);
  acl.fout.close();
  cout << gm2.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"2a");
  t.test_rel(x[1],0.0,4.0e-3,"2b");
  t.test_rel(x[2],0.0,4.0e-3,"2c");

  // Fletcher-Reeves conjugate with gradients
  acl.fout.open("ex_mmin2g.dat");
  x[0]=1.0;
  x[1]=0.0;
  x[2]=7.0*o2scl_const::pi;
  gm2.mmin_de(3,x,fmin,f1,f1g);
  acl.fout.close();
  cout << gm2.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"2a");
  t.test_rel(x[1],0.0,4.0e-3,"2b");
  t.test_rel(x[2],0.0,4.0e-3,"2c");

  // Polak-Ribere conjugate
  acl.fout.open("ex_mmin3.dat");
  x[0]=1.0;
  x[1]=0.0;
  x[2]=7.0*o2scl_const::pi;
  gm3.mmin(3,x,fmin,f1);
  acl.fout.close();
  cout << gm3.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"3a");
  t.test_rel(x[1],0.0,4.0e-3,"3b");
  t.test_rel(x[2],0.0,4.0e-3,"3c");

  // Polak-Ribere conjugate with gradients
  acl.fout.open("ex_mmin3g.dat");
  x[0]=1.0;
  x[1]=0.0;
  x[2]=7.0*o2scl_const::pi;
  gm3.mmin_de(3,x,fmin,f1,f1g);
  acl.fout.close();
  cout << gm3.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"3a");
  t.test_rel(x[1],0.0,4.0e-3,"3b");
  t.test_rel(x[2],0.0,4.0e-3,"3c");

  // BFGS method

  // BFGS has trouble converging (especially to zero, since the
  // minimimum of x[0] is exactly at zero) if the derivative is not
  // very accurate.
  gm4.def_grad.epsrel=1.0e-8;

  gm4.err_nonconv=false;
  acl.fout.open("ex_mmin4.dat");
  x[0]=1.0;
  x[1]=0.0;
  x[2]=7.0*o2scl_const::pi;
  gm4.mmin(3,x,fmin,f1);
  acl.fout.close();
  cout << gm4.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,1.0e-4,"4a");
  t.test_rel(x[1],0.0,1.0e-4,"4b");
  t.test_rel(x[2],0.0,1.0e-4,"4c");

  t.report();
  return 0;
}
// End of example

