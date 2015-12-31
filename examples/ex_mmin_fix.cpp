/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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

/* Example: ex_mmin_fix.cpp
   -------------------------------------------------------------------
   Example usage of the mmin_fix class, which fixes some of the
   paramters for a multidimensional minimization.
*/

#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/multi_funct.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/mmin_fix.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

class cl {

public:

  double mfn(size_t nv, const ubvector &x) {
    return (x[0]-2.0)*(x[0]-2.0)+(x[1]-1.0)*(x[1]-1.0)+x[2]*x[2];
  }
  
};

int main(void) {
  cl acl;
  ubvector x(3);
  double fmin;
  test_mgr t;

  t.set_output_level(1);
  cout.setf(ios::scientific);

  /*
    Perform the minimization the standard way, with the 
    simplex2 minimizer
  */
  multi_funct11 f1c11=
    std::bind(std::mem_fn<double(size_t,const ubvector &)>(&cl::mfn),
              acl,std::placeholders::_1,std::placeholders::_2);
  mmin_simp2<> gm1;
    
  x[0]=0.5;
  x[1]=0.5;
  x[2]=0.5;
  gm1.mmin(3,x,fmin,f1c11);
  cout << gm1.last_ntrial << " iterations." << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],2.0,1.0e-4,"1a");
  t.test_rel(x[1],1.0,1.0e-4,"1b");
  t.test_rel(x[2],0.0,1.0e-4,"1c");
  
  // Create a new mmin_fix_params object
  mmin_fix_params<> gmf;

  // Create a base minimizer which can be used by the mmin_fix_params
  // object. Note that we can't use 'gm1' here, because it has a
  // different type than 'gm2', even though its functionality is
  // effectively the same.
  mmin_simp2<mmin_fix_params<> > gm2;
  
  // Set the base minimizer
  gmf.set_mmin(gm2);

  /*
    First perform the minimization as above.
  */
  x[0]=0.5;
  x[1]=0.5;
  x[2]=0.5;
  gmf.mmin(3,x,fmin,f1c11);
  cout << gmf.last_ntrial << " iterations." << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],2.0,1.0e-4,"2a");
  t.test_rel(x[1],1.0,1.0e-4,"2b");
  t.test_rel(x[2],0.0,1.0e-4,"2c");

  /*
    Now fix the 2nd variable, and re-minimize.
  */
  bool fix[3]={false,true,false};
  x[0]=0.5;
  x[1]=0.5;
  x[2]=0.5;
  gmf.mmin_fix(3,x,fmin,fix,f1c11);
  cout << gmf.last_ntrial << " iterations." << endl;
  cout << "Found minimum at: " 
       << x[0] << " " << x[1] << " " << x[2] << endl;
  t.test_rel(x[0],2.0,1.0e-4,"3a");
  t.test_rel(x[1],0.5,1.0e-4,"3b");
  t.test_rel(x[2],0.0,1.0e-4,"3c");

  t.report();
  return 0;
}
// End of example
