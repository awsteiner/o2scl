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
#include <o2scl/table.h>
#include <o2scl/interp2_neigh.h>
#include <o2scl/interp2_planar.h>
#include <o2scl/rng.h>

using namespace std;
using namespace o2scl;

/*
typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef o2scl::matrix_view_table<> mat_x_t;
typedef const matrix_row_gen<mat_x_t> mat_x_row_t;
typedef const matrix_column_gen<mat_x_t> mat_x_col_t;
typedef o2scl::matrix_view_table_transpose<> mat_y_t;
typedef const matrix_row_gen<mat_y_t> mat_y_row_t;
typedef vector<function<double(mat_x_row_t &, mat_x_row_t &) > > f1_t;
typedef vector<function<double(mat_x_row_t &, const ubvector &) > > f2_t;
typedef vector<function<double(mat_y_row_t &, mat_y_row_t &) > > f3_t;
typedef vector<function<double(mat_y_row_t &, const ubvector &) > > f4_t;
*/

template<class vec_t, class vec2_t>
double covar(const vec_t &x, const vec2_t &y, double len) {
  double ret=exp(-(pow(x[0]-y[0],2.0)+pow(x[1]-y[1],2.0))/len/len/2.0);
  //cout << len << " " << x[0] << " " << y[0] << " "
  //<< x[1] << " " << y[1] << " " << ret << endl;
  return ret;
}

double ft(double x, double y) {
  return 3.0-2.0*x*x+7.0*y;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  {
    cout << "interpm_krige, unscaled" << endl;

    table<> tab;
    tab.line_of_names("x y z");
    tab.line_of_data(2,vector<double>({1.04,0.02}));
    tab.line_of_data(2,vector<double>({0.03,1.01}));
    tab.line_of_data(2,vector<double>({0.81,0.23}));
    tab.line_of_data(2,vector<double>({0.03,0.83}));
    tab.line_of_data(2,vector<double>({0.03,0.99}));
    tab.line_of_data(2,vector<double>({0.82,0.84}));
    tab.line_of_data(2,vector<double>({0.03,0.24}));
    tab.line_of_data(2,vector<double>({0.03,1.02}));

    for(size_t i=0;i<8;i++) {
      tab.set("z",i,ft(tab.get("x",i),tab.get("y",i)));
    }
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef o2scl::matrix_view_table<> mat_x_t;
    typedef const matrix_row_gen<mat_x_t> mat_x_row_t;
    typedef const matrix_column_gen<mat_x_t> mat_x_col_t;
    typedef o2scl::matrix_view_table_transpose<> mat_y_t;
    typedef const matrix_row_gen<mat_y_t> mat_y_row_t;
    typedef vector<function<double(mat_x_row_t &, mat_x_row_t &) > > f1_t;
    typedef vector<function<double(mat_x_row_t &, const ubvector &) > > f2_t;
    
    interpm_krige<ubvector,mat_x_t,mat_x_row_t,mat_x_col_t,
                  mat_y_t,mat_y_row_t,ubmatrix,
                  o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > ik;
    f1_t fa1={std::bind(&covar<mat_x_row_t,mat_x_row_t>,
                        std::placeholders::_1,std::placeholders::_2,0.70)};
    f2_t fa2={std::bind(&covar<mat_x_row_t,ubvector>,
                        std::placeholders::_1,std::placeholders::_2,0.70)};
    
    vector<string> col_list_x={"x","y"};
    vector<string> col_list_y={"z"};
    matrix_view_table<> mvt_x(tab,col_list_x);
    matrix_view_table_transpose<> mvt_y(tab,col_list_y);

    ik.set_data<>(2,1,8,mvt_x,mvt_y,fa1);
      
    ubvector point(2);
    ubvector out(1);

    point[0]=0.4;
    point[1]=0.5;
    ik.eval(point,out,fa2);
    t.test_rel(out[0],ft(point[0],point[1]),2.0e-1,"ik 0");
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    ik.eval(point,out,fa2);
    t.test_rel(out[0],ft(point[0],point[1]),1.0e-2,"ik 1");
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    cout << endl;

    ik.set_data<>(2,1,8,mvt_x,mvt_y,fa1,true);

    point[0]=0.4;
    point[1]=0.5;
    ik.eval(point,out,fa2);
    t.test_rel(out[0],ft(point[0],point[1]),5.0e-2,"ikr 0");
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    ik.eval(point,out,fa2);
    t.test_rel(out[0],ft(point[0],point[1]),1.0e-2,"ikr 1");
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    cout << endl;

    ik.unscale();
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,mat_x_col_t,
       mat_y_t,mat_y_row_t,ubmatrix,
       o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > iko;

    if (false) {
      ubvector len_precompute;
      iko.set_data<>(2,1,8,mvt_x,mvt_y,len_precompute);
      
      point[0]=0.4;
      point[1]=0.5;
      iko.eval(point,out,fa2);
      t.test_rel(out[0],ft(point[0],point[1]),2.0e-1,"iko 0");
      cout << out[0] << " " << ft(point[0],point[1]) << endl;
      point[0]=0.0301;
      point[1]=0.9901;
      iko.eval(point,out,fa2);
      t.test_rel(out[0],ft(point[0],point[1]),1.0e-2,"iko 1");
      cout << out[0] << " " << ft(point[0],point[1]) << endl;
      cout << endl;
    }

    
  }

#ifdef O2SCL_NEVER_DEFINED
  
  {

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef vector<ubvector> mat_x_t;
    typedef const matrix_row_gen<mat_x_t> mat_x_row_t;
    typedef const matrix_column_gen<mat_x_t> mat_x_col_t;
    typedef vector<ubvector> mat_y_t;
    typedef const matrix_row_gen<mat_y_t> mat_y_row_t;
    typedef vector<function<double(mat_x_row_t &, mat_x_row_t &) > > f1_t;
    typedef vector<function<double(mat_x_row_t &, const ubvector &) > > f2_t;
    
    interpm_krige<ubvector,mat_x_t,mat_x_row_t,mat_x_col_t,
                  mat_y_t,mat_y_row_t,ubmatrix,
                  o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > ik;
    
    cout << "interpm_krige_optim, rescaled" << endl;
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

    interpm_krige_optim<ubvector,mat_x_t,
			matrix_row_gen<mat_x_t> > iko;
    iko.verbose=2;

    /*
    iko.set_data<matrix_row_gen<mat_x_t> >(2,1,8,x2,y2,true);
    
    ubvector point(2);
    ubvector out(1);
    point[0]=0.4;
    point[1]=0.5;
    iko.eval(point,out,iko.ff2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    iko.eval(point,out,iko.ff2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    cout << endl;
    */
  }

  {
    cout << "interpm_krige_optim, rescaled, with data in a table" << endl;
    // Try a table representation
    table<> tab;
    tab.line_of_names("x y z");
    tab.line_of_data(2,vector<double>({1.04,0.02}));
    tab.line_of_data(2,vector<double>({0.03,1.01}));
    tab.line_of_data(2,vector<double>({0.81,0.23}));
    tab.line_of_data(2,vector<double>({0.03,0.83}));
    tab.line_of_data(2,vector<double>({0.03,0.99}));
    tab.line_of_data(2,vector<double>({0.82,0.84}));
    tab.line_of_data(2,vector<double>({0.03,0.24}));
    tab.line_of_data(2,vector<double>({0.03,1.02}));
    for(size_t i=0;i<8;i++) {
      tab.set("z",i,ft(tab.get("x",i),tab.get("y",i)));
    }
    matrix_view_table<> cmvtx(tab,{"x","y"});
    matrix_view_table_transpose<> cmvty(tab,{"z"});
  
    interpm_krige_optim<ubvector,mat_y_t,
			matrix_row_gen<mat_y_t> > iko;
    iko.verbose=2;

    /*
    iko.set_data<matrix_row_gen<matrix_view_table_transpose<> > >
      (2,1,8,cmvtx,cmvty,true);
    
    ubvector point(2);
    ubvector out(1);
    point[0]=0.4;
    point[1]=0.5;
    iko.eval(point,out,iko.ff2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    point[0]=0.0301;
    point[1]=0.9901;
    iko.eval(point,out,iko.ff2);
    cout << out[0] << " " << ft(point[0],point[1]) << endl;
    */
  }

  /*
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

#endif
  
  t.report();
  return 0;
}

