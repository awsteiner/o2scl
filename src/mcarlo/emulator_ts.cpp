/*
  -------------------------------------------------------------------
  
  Copyright (C) 2022, Andrew W. Steiner
  
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
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/emulator.h>

#ifdef O2SCL_ARMA
#include <armadillo>
#endif
#ifdef O2SCL_EIGEN
#include <eigen3/Eigen/Dense>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_linalg;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

double ft(double x, double y) {
  return 3.0-2.0*x*x+7.0*y;
}

void generate_table(table<> &tab, size_t N=100) {

  tab.clear();
  
  tab.line_of_names("x y z");
  
  for(size_t i=0;i<N;i++) {
    double x=3.0*sin(i*i);
    double y=5.0*cos(pow(i,4.0));
    vector<double> line={x,y,ft(x,y)};
    tab.line_of_data(3,line);
  }
  
  return;
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);

#ifdef O2SCL_NEVER_DEFINED
  
  vector<string> col_list={"x","y","z"};
    
  if (true) {
    
    table<> tab;
    generate_table(tab);

    emulator_interpm_idw_table<> em1;
    em1.set(2,1,2,tab,col_list);
    
    if (false) {
      for(size_t j=0;j<20;j++) {
        
        ubvector p(2);
        p[0]=1.0;
        p[1]=2.0;
        double z;
        double dz;
        ubvector dat(1), datu(1);
        em1.eval_unc(2,p,z,dz,dat,datu);
        cout << z << " " << dz << endl;
        exit(-1);
        //em1.eval_unc(
        //virtual int eval_unc(size_t n, const vec_t &p, double &log_wgt,
        //double &log_wgt_unc, vec2_t &dat, vec2_t &dat_unc) {
        
      }
    }
    
  }
  
  if (true) {
    
    table<> tab;
    generate_table(tab);

    emulator_interpm_krige_table<> em2;
    em2.set(2,1,2,tab,col_list);
    
    if (false) {
      for(size_t j=0;j<20;j++) {
        
        ubvector p(2);
        p[0]=1.0;
        p[1]=2.0;
        double z;
        double dz;
        ubvector dat(1), datu(1);
        em2.eval_unc(2,p,z,dz,dat,datu);
        cout << z << " " << dz << endl;
        exit(-1);
        //em1.eval_unc(
        //virtual int eval_unc(size_t n, const vec_t &p, double &log_wgt,
        //double &log_wgt_unc, vec2_t &dat, vec2_t &dat_unc) {
        
      }
    }
    
  }

#endif
  
  t.report();
  return 0;
}

