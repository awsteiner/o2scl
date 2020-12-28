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
#include <o2scl/diff_evo.h>
#include <o2scl/diff_evo_adapt.h>
#include <o2scl/rng_gsl.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

class cl {

public:

  cl() {
    param=30.0;
    rg.clock_seed();
  }

  // To output function evaluations to a file
  ofstream fout;

  // Parameter of the quadratic
  double param;
  
  rng_gsl rg;

  // Updated spring function
  double spring_two(size_t nv, const ubvector &x) {
    double theta=atan2(x[1],x[0]);
    double r=hypot(x[0],x[1]);
    double z=x[2];
    double tmz=theta-z;
    double fact=8.0-pow(sin(tmz+o2scl_const::pi/2.0)+1.0,3.0);
    double rm1=r-1.0;
    double ret=fact+exp(rm1*rm1)+z*z/param;
    fout << x[0] << " " << x[1] << " " << x[2] << " " << ret << " "
	 << fact << " " << rm1 << " " << endl;
    return ret;
  }

  // Gradient of the spring function
  int sgrad(size_t nv, ubvector &x, ubvector &g) {

    double theta=atan2(x[1],x[0]);
    double r=hypot(x[0],x[1]);
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

  // Set up the function objects
  multi_funct f1=std::bind
    (std::mem_fn<double(size_t,const ubvector &)>(&cl::spring_two),
     &acl,std::placeholders::_1,std::placeholders::_2);
  grad_funct f1g=std::bind
    (std::mem_fn<int(size_t,ubvector &,ubvector &)>(&cl::sgrad),
     &acl,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  mmin_simp2<> gm1;
  mmin_conf<> gm2;
  mmin_conp<> gm3;
  mmin_bfgs2<> gm4;
  diff_evo<> gm5;
  diff_evo_adapt<> gm6;

  vector<double> guess={2.0,1.0,7.0*o2scl_const::pi};
  
  // This function is difficult to minimize, so more trials
  // are required.
  gm1.ntrial*=10;
  gm2.ntrial*=10;
  gm3.ntrial*=10;
  gm4.ntrial*=10;

  // Simplex minimization
  acl.fout.open("ex_mmin1.dat");
  vector_copy(3,guess,x);
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
  vector_copy(3,guess,x);
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
  vector_copy(3,guess,x);
  gm2.mmin_de(3,x,fmin,f1,f1g);
  acl.fout.close();
  cout << gm2.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"2ga");
  t.test_rel(x[1],0.0,4.0e-3,"2gb");
  t.test_rel(x[2],0.0,4.0e-3,"2gc");

  // Polak-Ribere conjugate
  acl.fout.open("ex_mmin3.dat");
  vector_copy(3,guess,x);
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
  vector_copy(3,guess,x);
  gm3.mmin_de(3,x,fmin,f1,f1g);
  acl.fout.close();
  cout << gm3.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"3ga");
  t.test_rel(x[1],0.0,4.0e-3,"3gb");
  t.test_rel(x[2],0.0,4.0e-3,"3gc");

  // de
  acl.fout.open("ex_mmin5.dat");
  vector_copy(3,guess,x);
  gm5.mmin(3,x,fmin,f1);
  acl.fout.close();
  cout << gm5.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"5a");
  t.test_rel(x[1],0.0,4.0e-3,"5b");
  t.test_rel(x[2],0.0,4.0e-3,"5c");

  // dea
  acl.fout.open("ex_mmin6.dat");
  vector_copy(3,guess,x);
  gm6.mmin(3,x,fmin,f1);
  acl.fout.close();
  cout << gm6.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"6a");
  t.test_rel(x[1],0.0,4.0e-3,"6b");
  t.test_rel(x[2],0.0,4.0e-3,"6c");

  t.report();
  return 0;
}

/*
  The BFGS minimizer doesn't appear to work for this particular
  example. This may be a result of finite precision in the object
  function rather than a failure of the BFGS.

  // BFGS with gradients
  acl.fout.open("ex_mmin4g.dat");
  vector_copy(3,guess,x);
  gm4.mmin_de(3,x,fmin,f1,f1g);
  acl.fout.close();
  cout << gm4.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"4ga");
  t.test_rel(x[1],0.0,4.0e-3,"4gb");
  t.test_rel(x[2],0.0,4.0e-3,"4gc");

  gm4.def_grad.epsrel=1.0e-8;
  
  acl.fout.open("ex_mmin4.dat");
  vector_copy(3,guess,x);
  gm4.mmin(3,x,fmin,f1);
  acl.fout.close();
  cout << gm4.last_ntrial << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],1.0,4.0e-3,"4a");
  t.test_rel(x[1],0.0,4.0e-3,"4b");
  t.test_rel(x[2],0.0,4.0e-3,"4c");

*/
