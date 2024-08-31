/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/test_mgr.h>
#include <o2scl/exp_max.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

#ifdef O2SCL_SET_ARMA
#include <armadillo>
#endif
#ifdef O2SCL_SET_EIGEN
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

  if (false) {
    
    ubvector mean1(2), mean2(2);
    ubmatrix covar1(2,2), covar2(2,2);

    mean1[0]=0.75;
    mean1[1]=0.6;
    mean2[0]=0.25;
    mean2[1]=0.25;
    double covar=1.6e-3;
    covar1(0,0)=covar;
    covar1(1,1)=covar;
    covar1(0,1)=covar*0.9;
    covar1(1,0)=covar*0.9;
    covar2(0,0)=covar;
    covar2(1,1)=covar;
    covar2(0,1)=0.0;
    covar2(1,0)=0.0;

    prob_dens_mdim_gaussian<> pdmg1, pdmg2;
    pdmg1.set_covar(2,mean1,covar1);
    pdmg2.set_covar(2,mean2,covar2);

    table<> tab;
    tab.line_of_names("x y");

    ubvector tvec(2);
    for(size_t i=0;i<80;i++) {
      pdmg1(tvec);
      tab.line_of_data(2,tvec);
    }
    for(size_t i=0;i<40;i++) {
      pdmg2(tvec);
      tab.line_of_data(2,tvec);
    }

    hdf_file hf;
    hf.open_or_create("em.o2");
    hdf_output(hf,tab,"em");
    hf.close();
  
    exp_max_gmm<> emg;
    const_matrix_view_table<> cmvt(tab,{"x","y"});
    emg.set_data(2,tab.get_nlines(),cmvt);
    
    emg.calc_auto(2);
    const ubvector &test_mean1=emg.get_gmm().pdmg[0].get_peak();
    const ubvector &test_mean2=emg.get_gmm().pdmg[1].get_peak();
    if (test_mean1[0]>test_mean2[0]) {
      t.test_rel(test_mean1[0],0.75,1.0e-2,"mean 1");
      t.test_rel(test_mean1[1],0.6,1.0e-2,"mean 2");
      //t.test_rel(test_mean2[0],0.25,2.0e-2,"mean 3");
      //t.test_rel(test_mean2[1],0.25,2.0e-2,"mean 4");
    } else {
      t.test_rel(test_mean2[0],0.75,1.0e-2,"mean 1b");
      t.test_rel(test_mean2[1],0.6,1.0e-2,"mean 2b");
      //t.test_rel(test_mean1[0],0.25,2.0e-2,"mean 3b");
      //t.test_rel(test_mean1[1],0.25,2.0e-2,"mean 4b");
    }

  }

#ifdef O2SCL_SET_EIGEN

  if (false) {
  
    ubvector mean1(2), mean2(2);
    ubmatrix covar1(2,2), covar2(2,2);

    mean1[0]=0.75;
    mean1[1]=0.6;
    mean2[0]=0.25;
    mean2[1]=0.25;
    double covar=1.6e-3;
    covar1(0,0)=covar;
    covar1(1,1)=covar;
    covar1(0,1)=covar*0.9;
    covar1(1,0)=covar*0.9;
    covar2(0,0)=covar;
    covar2(1,1)=covar;
    covar2(0,1)=0.0;
    covar2(1,0)=0.0;

    prob_dens_mdim_gaussian<> pdmg1, pdmg2;
    pdmg1.set_covar(2,mean1,covar1);
    pdmg2.set_covar(2,mean2,covar2);

    table<> tab;
    tab.line_of_names("x y");

    ubvector tvec(2);
    for(size_t i=0;i<80;i++) {
      pdmg1(tvec);
      tab.line_of_data(2,tvec);
    }
    for(size_t i=0;i<40;i++) {
      pdmg2(tvec);
      tab.line_of_data(2,tvec);
    }

    hdf_file hf;
    hf.open_or_create("em.o2");
    hdf_output(hf,tab,"em");
    hf.close();

    /*

      AWS, 2/27/23: this doesn't work. We need to fix exp_max_gmm
      to allow different vector types for the underlying 
      prob_dens_mdim_gmm object.

      exp_max_gmm<const_matrix_view_table<>,Eigen::VectorXd,
      Eigen::MatrixXd> emg;
      
      const_matrix_view_table<> cmvt(tab,{"x","y"});
      emg.set_data(2,tab.get_nlines(),cmvt);
      
      emg.calc_auto(2);
      const Eigen::VectorXd &test_mean1=emg.get_gmm().pdmg[0].get_peak();
      const Eigen::VectorXd &test_mean2=emg.get_gmm().pdmg[1].get_peak();
      if (test_mean1[0]>test_mean2[0]) {
      t.test_rel(test_mean1[0],0.75,1.0e-2,"mean 1");
      t.test_rel(test_mean1[1],0.6,1.0e-2,"mean 2");
      t.test_rel(test_mean2[0],0.25,2.0e-2,"mean 3");
      t.test_rel(test_mean2[1],0.25,2.0e-2,"mean 4");
      } else {
      t.test_rel(test_mean2[0],0.75,1.0e-2,"mean 1b");
      t.test_rel(test_mean2[1],0.6,1.0e-2,"mean 2b");
      t.test_rel(test_mean1[0],0.25,2.0e-2,"mean 3b");
      t.test_rel(test_mean1[1],0.25,2.0e-2,"mean 4b");
      }
    */

  }
  
#endif
  
  t.report();
  return 0;
}

