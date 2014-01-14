/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <climits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/vector.h>
#include <o2scl/columnify.h>
#include <o2scl/test_mgr.h>
#include <o2scl/permutation.h>
#include <o2scl/tensor_grid.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  {
    ubmatrix ub1(3,3);
    ub1(0,0)=0.0;
    ub1(0,1)=1.0;
    ub1(0,2)=2.0;
    ub1(1,0)=1.0;
    ub1(1,1)=2.0;
    ub1(1,2)=3.0;
    ub1(2,0)=2.0;
    ub1(2,1)=3.0;
    ub1(2,2)=4.0;
    matrix_row_gen<ubmatrix> r1=
      o2scl::matrix_row<ubmatrix,matrix_row_gen<ubmatrix> >(ub1,2);
    t.test_rel(r1[0],2.0,1.0e-12,"matrix row 1");
    t.test_rel(r1[1],3.0,1.0e-12,"matrix row 2");
    t.test_rel(r1[2],4.0,1.0e-12,"matrix row 3");
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    r1[0]=-1.0;
    r1[1]=-2.0;
    r1[2]=-3.0;
    cout << ub1(2,0) << " " << ub1(2,1) << " " << ub1(2,2) << endl;
    t.test_rel(ub1(2,0),-1.0,1.0e-12,"matrix row 4");
    t.test_rel(ub1(2,1),-2.0,1.0e-12,"matrix row 5");
    t.test_rel(ub1(2,2),-3.0,1.0e-12,"matrix row 6");
    ub1(2,0)=4.0;
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    t.test_rel(r1[0],4.0,1.0e-12,"matrix row 7");
    matrix_column_gen<ubmatrix> r2=
      o2scl::matrix_column<ubmatrix,matrix_column_gen<ubmatrix> >(ub1,2);
    t.test_rel(r2[0],2.0,1.0e-12,"matrix col 1");
    t.test_rel(r2[1],3.0,1.0e-12,"matrix col 2");
    t.test_rel(r2[2],-3.0,1.0e-12,"matrix col 3");
    r2[0]=1.0;
    t.test_rel(ub1(0,2),1.0,1.0e-12,"matrix col 4");

  }
  
  {
    ubmatrix ub1(3,3);
    ub1(0,0)=0.0;
    ub1(0,1)=1.0;
    ub1(0,2)=2.0;
    ub1(1,0)=1.0;
    ub1(1,1)=2.0;
    ub1(1,2)=3.0;
    ub1(2,0)=2.0;
    ub1(2,1)=3.0;
    ub1(2,2)=4.0;
    ubmatrix_row r1=o2scl::matrix_row<ubmatrix,ubmatrix_row>(ub1,2);
    t.test_rel(r1[0],2.0,1.0e-12,"matrix row 1");
    t.test_rel(r1[1],3.0,1.0e-12,"matrix row 2");
    t.test_rel(r1[2],4.0,1.0e-12,"matrix row 3");
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    r1[0]=-1.0;
    r1[1]=-2.0;
    r1[2]=-3.0;
    cout << ub1(2,0) << " " << ub1(2,1) << " " << ub1(2,2) << endl;
    t.test_rel(ub1(2,0),-1.0,1.0e-12,"matrix row 4");
    t.test_rel(ub1(2,1),-2.0,1.0e-12,"matrix row 5");
    t.test_rel(ub1(2,2),-3.0,1.0e-12,"matrix row 6");
    ub1(2,0)=4.0;
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    t.test_rel(r1[0],4.0,1.0e-12,"matrix row 7");
    ubmatrix_column r2=o2scl::matrix_column<ubmatrix,ubmatrix_column>(ub1,2);
    t.test_rel(r2[0],2.0,1.0e-12,"matrix col 1");
    t.test_rel(r2[1],3.0,1.0e-12,"matrix col 2");
    t.test_rel(r2[2],-3.0,1.0e-12,"matrix col 3");
    r2[0]=1.0;
    t.test_rel(ub1(0,2),1.0,1.0e-12,"matrix col 4");
  }

#ifdef O2SCL_ARMA
  {
    arma::mat ub1(3,3);
    ub1(0,0)=0.0;
    ub1(0,1)=1.0;
    ub1(0,2)=2.0;
    ub1(1,0)=1.0;
    ub1(1,1)=2.0;
    ub1(1,2)=3.0;
    ub1(2,0)=2.0;
    ub1(2,1)=3.0;
    ub1(2,2)=4.0;
    arma::subview_row<double> r1=
      matrix_row<arma::mat,arma::subview_row<double> >(ub1,2);
    t.test_rel(r1[0],2.0,1.0e-12,"matrix row 1");
    t.test_rel(r1[1],3.0,1.0e-12,"matrix row 2");
    t.test_rel(r1[2],4.0,1.0e-12,"matrix row 3");
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    r1[0]=-1.0;
    r1[1]=-2.0;
    r1[2]=-3.0;
    cout << ub1(2,0) << " " << ub1(2,1) << " " << ub1(2,2) << endl;
    t.test_rel(ub1(2,0),-1.0,1.0e-12,"matrix row 4");
    t.test_rel(ub1(2,1),-2.0,1.0e-12,"matrix row 5");
    t.test_rel(ub1(2,2),-3.0,1.0e-12,"matrix row 6");
    ub1(2,0)=4.0;
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    t.test_rel(r1[0],4.0,1.0e-12,"matrix row 7");
    arma::subview_col<double> r2=
      matrix_column<arma::mat,arma::subview_col<double> >(ub1,2);
    t.test_rel(r2[0],2.0,1.0e-12,"matrix col 1");
    t.test_rel(r2[1],3.0,1.0e-12,"matrix col 2");
    t.test_rel(r2[2],-3.0,1.0e-12,"matrix col 3");
    r2[0]=1.0;
    t.test_rel(ub1(0,2),1.0,1.0e-12,"matrix col 4");
  }
#endif

#ifdef O2SCL_EIGEN
  {
    Eigen::MatrixXd ub1(3,3);
    ub1(0,0)=0.0;
    ub1(0,1)=1.0;
    ub1(0,2)=2.0;
    ub1(1,0)=1.0;
    ub1(1,1)=2.0;
    ub1(1,2)=3.0;
    ub1(2,0)=2.0;
    ub1(2,1)=3.0;
    ub1(2,2)=4.0;
    Eigen::MatrixXd::RowXpr r1=
    matrix_row<Eigen::MatrixXd,Eigen::MatrixXd::RowXpr>(ub1,2);
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    t.test_rel(r1[0],2.0,1.0e-12,"matrix row 1");
    t.test_rel(r1[1],3.0,1.0e-12,"matrix row 2");
    t.test_rel(r1[2],4.0,1.0e-12,"matrix row 3");
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    r1[0]=-1.0;
    r1[1]=-2.0;
    r1[2]=-3.0;
    cout << ub1(2,0) << " " << ub1(2,1) << " " << ub1(2,2) << endl;
    t.test_rel(ub1(2,0),-1.0,1.0e-12,"matrix row 4");
    t.test_rel(ub1(2,1),-2.0,1.0e-12,"matrix row 5");
    t.test_rel(ub1(2,2),-3.0,1.0e-12,"matrix row 6");
    ub1(2,0)=4.0;
    cout << r1[0] << " " << r1[1] << " " << r1[2] << endl;
    t.test_rel(r1[0],4.0,1.0e-12,"matrix row 7");
    Eigen::MatrixXd::ColXpr r2=
      matrix_column<Eigen::MatrixXd,Eigen::MatrixXd::ColXpr>(ub1,2);
    t.test_rel(r2[0],2.0,1.0e-12,"matrix col 1");
    t.test_rel(r2[1],3.0,1.0e-12,"matrix col 2");
    t.test_rel(r2[2],-3.0,1.0e-12,"matrix col 3");
    r2[0]=1.0;
    t.test_rel(ub1(0,2),1.0,1.0e-12,"matrix col 4");
  }
#endif

  // Test vector_sort
  double st1[5]={3.1,4.1,5.9,2.6,5.3};
  vector_sort<double[5],double>(5,st1);
  vector_out(cout,5,st1,true);
  t.test_rel(st1[0],2.6,1.0e-10,"sort 1");
  t.test_rel(st1[1],3.1,1.0e-10,"sort 1");
  t.test_rel(st1[2],4.1,1.0e-10,"sort 1");
  t.test_rel(st1[3],5.3,1.0e-10,"sort 1");
  t.test_rel(st1[4],5.9,1.0e-10,"sort 1");

  // Test vector_sort_double
  double *st2=new double[5];
  st2[0]=3.1;
  st2[1]=4.1;
  st2[2]=5.9;
  st2[3]=2.6;
  st2[4]=5.3;
  vector_sort_double(5,st2);
  vector_out(cout,5,st2,true);
  t.test_rel(st2[0],2.6,1.0e-10,"sort 2");
  t.test_rel(st2[1],3.1,1.0e-10,"sort 2");
  t.test_rel(st2[2],4.1,1.0e-10,"sort 2");
  t.test_rel(st2[3],5.3,1.0e-10,"sort 2");
  t.test_rel(st2[4],5.9,1.0e-10,"sort 2");
  delete[] st2;

  // Test vector_sort_index
  double *st3=new double[5];
  size_t ind[5];
  st3[0]=3.1;
  st3[1]=4.1;
  st3[2]=5.9;
  st3[3]=2.6;
  st3[4]=5.3;
  vector_sort_index(5,st3,ind);
  vector_out(cout,5,st3,true);
  vector_out(cout,5,ind,true);
  t.test_gen(ind[0]==3,"sort 3");
  t.test_gen(ind[1]==0,"sort 3");
  t.test_gen(ind[2]==1,"sort 3");
  t.test_gen(ind[3]==4,"sort 3");
  t.test_gen(ind[4]==2,"sort 3");

  // Show that it works with a permutation
  double *st4=new double[5];
  permutation p(5);
  st4[0]=3.1;
  st4[1]=4.1;
  st4[2]=5.9;
  st4[3]=2.6;
  st4[4]=5.3;
  vector_sort_index(5,st4,p);
  vector_out(cout,5,st4,true);
  vector_out(cout,5,p,true);
  t.test_gen(p[0]==3,"sort 3");
  t.test_gen(p[1]==0,"sort 3");
  t.test_gen(p[2]==1,"sort 3");
  t.test_gen(p[3]==4,"sort 3");
  t.test_gen(p[4]==2,"sort 3");

  // Compute the inverse
  permutation p2=p.inverse();
  vector_out(cout,5,p2,true);

  // Test vector_grid
  ubvector ovg(11);
  vector_grid(uniform_grid_end<>(0,1,10),ovg);
  //cout << ovg << endl;

  t.report();
  return 0;
}

