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
#include <o2scl/interpm_krige.h>
#include <o2scl/interp_krige.h>
#include <o2scl/table.h>
#include <o2scl/interp2_neigh.h>
#include <o2scl/interp2_planar.h>
#include <o2scl/rng.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

double f(double x, double mean, double sd) {
  return (sin(1.0/(0.3+x))-mean)/sd;
}

template<class vec_t, class vec2_t>
double covar(const vec_t &x, const vec2_t &y, double len) {
  double ret=exp(-(pow(x[0]-y[0],2.0)+pow(x[1]-y[1],2.0))/len/len/2.0);
  return ret;
}

template<class vec_t, class vec2_t>
double covar_deriv(const vec_t &x, const vec2_t &y, size_t i,
                   double len) {
  double ret;
  if (i==0) {
    ret=-exp(-(pow(x[0]-y[0],2.0)+pow(x[1]-y[1],2.0))/len/len/2.0)/
      len/len*(x[0]-y[0]);
  } else {
    ret=-exp(-(pow(x[0]-y[0],2.0)+pow(x[1]-y[1],2.0))/len/len/2.0)/
      len/len*(x[1]-y[1]);
  }
  return ret;
}

double ft(double x, double y) {
  return 3.0-2.0*x*x+7.0*y;
}

double ftdx(double x, double y) {
  return -4.0*x;
}

double ftdy(double x, double y) {
  return 7.0;
}

void generate_table(table<> &tab) {

  tab.clear();
  
  tab.line_of_names("x y z");
  
  for(size_t i=0;i<100;i++) {
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
  
  vector<string> col_list_x={"x","y"};
  vector<string> col_list_y={"z"};
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef o2scl::matrix_view_table<> mat_x_t;
  typedef const matrix_row_gen<mat_x_t> mat_x_row_t;
  typedef const matrix_column_gen<mat_x_t> mat_x_col_t;
  typedef o2scl::matrix_view_table_transpose<> mat_y_t;
  typedef const matrix_row_gen<mat_y_t> mat_y_row_t;
  typedef vector<function<double(mat_x_row_t &, mat_x_row_t &) > > f1_t;
  typedef vector<function<double(mat_x_row_t &, const ubvector &) > > f2_t;
    
  {
    cout << "interpm_krige, not rescaled" << endl;
    
    table<> tab;
    generate_table(tab);
    
    hdf_file hf;
    hf.open_or_create("interpm_krige_ts_data.o2");
    hdf_output(hf,tab,"tab");
    hf.close();
    
    interpm_krige<ubvector,mat_x_t,mat_x_row_t,mat_x_col_t,
                  mat_y_t,mat_y_row_t,ubmatrix,f1_t,
                  o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > ik;
    
    f1_t fa1={std::bind(&covar<mat_x_row_t,mat_x_row_t>,
                        std::placeholders::_1,std::placeholders::_2,1.1)};
    f2_t fa2={std::bind(&covar<mat_x_row_t,ubvector>,
                        std::placeholders::_1,std::placeholders::_2,1.1)};
    
    matrix_view_table<> mvt_x(tab,col_list_x);
    matrix_view_table_transpose<> mvt_y(tab,col_list_y);

    ik.verbose=2;
    ik.set_data(2,1,tab.get_nlines(),mvt_x,mvt_y,fa1);

    gen_test_number<> gtn_x, gtn_y;
    gtn_x.set_radix(1.9);
    gtn_y.set_radix(1.9);

    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x.gen();
      point[1]=gtn_y.gen();
        
      ik.eval(point,out,fa2);
      cout.setf(ios::showpos);
      cout << point[0] << " " << point[1] << " "
           << out[0] << " " << ft(point[0],point[1]) << endl;
      cout.unsetf(ios::showpos);
      t.test_rel(out[0],ft(point[0],point[1]),1.0,"unscaled 1");

    }
    cout << endl;

  }

  {
    
    // Now with rescaled=true
    cout << "interpm_krige, rescaled" << endl;
    table<> tab2;
    
    generate_table(tab2);
    
    matrix_view_table<> mvt_x2(tab2,col_list_x);
    matrix_view_table_transpose<> mvt_y2(tab2,col_list_y);

    gen_test_number<> gtn_x2, gtn_y2;
    gtn_x2.set_radix(1.9);
    gtn_y2.set_radix(1.9);
    
    interpm_krige<ubvector,mat_x_t,mat_x_row_t,mat_x_col_t,
                  mat_y_t,mat_y_row_t,ubmatrix,f1_t,
                  o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > ik;
    
    f1_t fa1={std::bind(&covar<mat_x_row_t,mat_x_row_t>,
                        std::placeholders::_1,std::placeholders::_2,1.1)};
    f2_t fa2={std::bind(&covar<mat_x_row_t,ubvector>,
                        std::placeholders::_1,std::placeholders::_2,1.1)};
    
    ik.verbose=2;
    ik.set_data(2,1,tab2.get_nlines(),mvt_x2,mvt_y2,fa1,true);

    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x2.gen();
      point[1]=gtn_y2.gen();
      
      ik.eval(point,out,fa2);
      cout.setf(ios::showpos);
      cout << point[0] << " " << point[1] << " "
           << out[0] << " " << ft(point[0],point[1]) << endl;
      cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),1.0,"rescaled 1");
    }
    cout << endl;

  }
  {
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,mat_x_col_t,
       mat_y_t,mat_y_row_t,ubmatrix,
       o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > iko;

    table<> tab3;
    generate_table(tab3);
    
    matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table_transpose<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3, gtn_y3;
    gtn_x3.set_radix(1.9);
    gtn_y3.set_radix(1.9);
    
    iko.verbose=2;
    iko.nlen=200;
    ubvector len_precompute;
    iko.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3,
                len_precompute);
    
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_y3.gen();
      
      iko.eval2(point,out);
      cout.setf(ios::showpos);
      cout << point[0] << " " << point[1] << " "
           << out[0] << " " << ft(point[0],point[1]) << endl;
      cout.unsetf(ios::showpos);
      t.test_rel(out[0],ft(point[0],point[1]),4.0e-2,"unscaled 2");

    }
    cout << endl;
    
  }

  /*
    AWS: 5/11/22: Old interpm_krige_nn test which is not working

    if (false) {
    // Construct the data
    vector<ubvector> x;
    ubvector tmp(2);
    tmp[0]=1.04; tmp[1]=0.02;
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.01; 
    x.push_back(tmp);
    tmp[0]=0.81; tmp[1]=0.23; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.83; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.99; 
    x.push_back(tmp);
    tmp[0]=0.82; tmp[1]=0.84; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=0.24; 
    x.push_back(tmp);
    tmp[0]=0.03; tmp[1]=1.02; 
    x.push_back(tmp);
    mat_x_t x2(x);

    vector<ubvector> y;
    tmp.resize(8);
    for(size_t i=0;i<8;i++) {
    tmp[i]=ft(x[i][0],x[i][1]);
    }
    y.push_back(tmp);
    mat_x_t y2(y);

    interpm_krige_nn<ubvector,mat_x_t,
    matrix_row_gen<mat_x_t> > ik;
    vector<function<double(const ubvector &,const ubvector &)> > fa={covar};
    ik.set_data(2,1,8,x2,y2,fa,4);
  
    cout << "Interpolate at a point and compare the three methods:" << endl;
    ubvector point(2);
    ubvector out(1);
    point[0]=0.4;
    point[1]=0.5;
    ik.eval(point,out,fa);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    ik.eval(point,out,fa);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;

    }
  */

  t.report();
  return 0;
}

