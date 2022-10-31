/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/exp_max.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

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

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout.setf(ios::scientific);
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  ubvector mean1(2), mean2(2);
  ubmatrix covar1(2,2), covar2(2,2);

  mean1[0]=0.75;
  mean1[1]=0.6;
  mean2[0]=0.25;
  mean2[1]=0.25;
  covar1(0,0)=0.2;
  covar1(1,1)=0.2;
  covar1(0,1)=0.039;
  covar1(1,0)=0.039;
  covar2(0,0)=0.2;
  covar2(1,1)=0.2;
  covar2(0,1)=0.0;
  covar2(1,0)=0.0;

  prob_dens_mdim_gaussian<> pdmg1, pdmg2;
  pdmg1.set_covar(2,mean1,covar1);
  pdmg2.set_covar(2,mean2,covar2);
  
  table<> tab;
  tab.line_of_names("x y");

  ubvector tvec(2);
  for(size_t i=0;i<40;i++) {
    pdmg1(tvec);
    tab.line_of_data(2,tvec);
  }
  for(size_t i=0;i<20;i++) {
    pdmg2(tvec);
    tab.line_of_data(2,tvec);
  }

  exp_max_gmm<> emg;
  const_matrix_view_table<> cmvt(tab,{"x","y"});
  emg.set_data(2,tab.get_nlines(),cmvt);

  emg.compute_auto(2);

  t.report();
  return 0;
}

