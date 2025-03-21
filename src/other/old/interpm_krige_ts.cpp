/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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

double f(double x, double mean, double sd) {
  return (sin(1.0/(0.3+x))-mean)/sd;
}

template<class vec_t, class vec2_t>
double covar(size_t j, const vec_t &x, const vec2_t &y, double len) {
  double ret=exp(-(pow(x[0]-y[0],2.0)+pow(x[1]-y[1],2.0))/len/len/2.0);
  return ret;
}

template<class vec_t, class vec2_t>
double covar_deriv(size_t j, const vec_t &x, const vec2_t &y, size_t i,
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
  
  vector<string> col_list_x={"x","y"};
  vector<string> col_list_y={"z"};
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef o2scl::matrix_view_table<> mat_x_t;
  typedef const matrix_row_gen<mat_x_t> mat_x_row_t;
  typedef o2scl::matrix_view_table_transpose<> mat_y_t;
  typedef const matrix_row_gen<mat_y_t> mat_y_row_t;
  
  typedef interpm_krige<ubvector,mat_x_t,mat_x_row_t,
                        mat_y_t,mat_y_row_t,ubmatrix,
                        o2scl_linalg::matrix_invert_det_cholesky
                        <ubmatrix> > ik_t;

  {

    mcovar_funct_rbf_noise mfrn;
    
    interpm_krige_optim_new
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_row_t,ubmatrix,
       o2scl_linalg::matrix_invert_det_cholesky<ubmatrix>,
       std::vector<std::vector<double>> > ikon;

    table<> tab3;
    generate_table(tab3);
    
    matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table_transpose<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    ikon.verbose=2;
    vector<double> len_list={0.01,0.03,0.1,0.3,1.0};
    vector<double> l10_list={-15,-13,-11,-9};
    vector<vector<double> > param_lists;
    param_lists.push_back(len_list);
    param_lists.push_back(len_list);
    param_lists.push_back(l10_list);
    
    ikon.set_param_lists(param_lists);
    ikon.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3,
                  [mfrn](auto &&i, auto &&t, auto &&t2) mutable {
                    return mfrn(i,t,t2); });
    
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        ikon.eval(mfrn,point,out);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),4.0e-1,"unscaled eigen 2");
      }

    }
    cout << endl;
    
  }
  
  {
    cout << "interpm_krige, not rescaled" << endl;
    
    table<> tab;
    generate_table(tab);
    
    hdf_file hf;
    hf.open_or_create("interpm_krige_ts_data.o2");
    hdf_output(hf,tab,"tab");
    hf.close();

    ik_t ik;
    
    matrix_view_table<> mvt_x(tab,col_list_x);
    matrix_view_table_transpose<> mvt_y(tab,col_list_y);
    
    ik.verbose=2;
    ik.set_data(2,1,tab.get_nlines(),mvt_x,mvt_y,
                [](auto &&i, auto &&t, auto &&t2) mutable {
                  return covar(i,t,t2,1.1); });

    gen_test_number<> gtn_x;
    gtn_x.set_radix(1.9);

    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1), sig(1);
      point[0]=gtn_x.gen();
      point[1]=gtn_x.gen();

      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        ik.eval_covar([](auto &&i, auto &&t, auto &&t2) mutable {
          return covar(i,t,t2,1.1); },point,out);
        ik.sigma_covar([](auto &&i, auto &&t, auto &&t2) mutable {
          return covar(i,t,t2,1.1); },point,sig);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),1.0,"unscaled 1");
      }

    }
    cout << endl;

  }

#ifdef O2SCL_NEVER_DEFINED  

  {
    
    // Now with rescaled=true
    cout << "interpm_krige, rescaled" << endl;
    table<> tab2;
    
    generate_table(tab2);
    
    matrix_view_table<> mvt_x2(tab2,col_list_x);
    matrix_view_table_transpose<> mvt_y2(tab2,col_list_y);

    gen_test_number<> gtn_x2;
    gtn_x2.set_radix(1.9);

    ik_t ik;
    
    ik_t::f1_t fa1={std::bind(&covar<mat_x_row_t,mat_x_row_t>,
                        std::placeholders::_1,std::placeholders::_2,1.1)};
    ik_t::f2_t fa2={std::bind(&covar<mat_x_row_t,ubvector>,
                              std::placeholders::_1,
                              std::placeholders::_2,1.1)};
    
    ik.verbose=2;
    ik.set_data(2,1,tab2.get_nlines(),mvt_x2,mvt_y2,fa1,true);

    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x2.gen();
      point[1]=gtn_x2.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        ik.eval_covar(point,out,fa2);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),1.0,"rescaled 1");
      }
    }
    cout << endl;

  }
  
  for(size_t k=0;k<4;k++) {
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_row_t,ubmatrix,
       o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > iko;

    if (k==0) {
      cout << "lml unscaled" << endl;
    } else if (k==1) {
      cout << "loo_cv unscaled" << endl;
    } else if (k==2) {
      cout << "lml rescaled" << endl;
    } else {
      cout << "loo_cv rescaled" << endl;
    }

    table<> tab3;
    generate_table(tab3);
    
    matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table_transpose<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    iko.verbose=1;
    iko.nlen=50;
    if (k==1 || k==3) iko.mode=iko.mode_loo_cv;
    if (k==2 || k==3) {
      iko.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3,true);
    } else {
      iko.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3);
    }
    
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1), sig(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        iko.eval(point,out);
        iko.sigma(point,sig);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        if (k==0) {
          t.test_rel(out[0],ft(point[0],point[1]),4.0e-1,
                     "unscaled lml 2");
        } else if (k==1) {
          t.test_rel(out[0],ft(point[0],point[1]),1.0e-1,
                     "unscaled loo_cv 2");
        } else if (k==2) {
          t.test_rel(out[0],ft(point[0],point[1]),1.0e-2,
                     "rescaled lml 2");
        } else {
          t.test_rel(out[0],ft(point[0],point[1]),1.0e-2,
                     "rescaled loo_cv 2");
        }
      }

    }
    cout << endl;
    
  }

  // More refined data set to test derivatives
  
  if (true) {
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_row_t,ubmatrix,
       o2scl_linalg::matrix_invert_det_cholesky<ubmatrix> > iko;

    table<> tab3;
    generate_table(tab3);
    
    matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table_transpose<> mvt_y3(tab3,col_list_y);
    
    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    iko.verbose=1;
    iko.nlen=50;
    iko.timing=true;
    iko.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3,true);
    
    for(size_t j=0;j<20;j++) {
      
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      cout.setf(ios::showpos);
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        cout << point[0] << " " << point[1] << endl;
        iko.eval(point,out);
        cout << "  " << out[0] << " " << ft(point[0],point[1]) << endl;
        t.test_rel(out[0],ft(point[0],point[1]),4.0e-1,"deriv");
        
        iko.deriv(point,out,0);
        cout << "  " << out[0] << " " << ftdx(point[0],point[1]) << endl;
        
        iko.deriv(point,out,1);
        cout << "  " << out[0] << " " << ftdy(point[0],point[1]) << endl;
        
        cout.unsetf(ios::showpos);
      }
      
    }
    cout << endl;
    
  }

#ifdef O2SCL_ARMA  
  
  {
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_row_t,arma::mat,
       matrix_invert_det_sympd_arma<> > iko_arma;

    table<> tab3;
    generate_table(tab3);
    
    matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table_transpose<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    iko_arma.verbose=2;
    iko_arma.nlen=50;
    iko_arma.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3);
    
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        iko_arma.eval(point,out);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),4.0e-1,"unscaled arma 2");
      }

    }
    cout << endl;
    
  }

#endif  

#ifdef O2SCL_EIGEN  
  
  {
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_row_t,Eigen::MatrixXd,
       matrix_invert_det_eigen<> > iko_eigen;

    table<> tab3;
    generate_table(tab3);
    
    matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table_transpose<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    iko_eigen.verbose=2;
    iko_eigen.nlen=50;
    iko_eigen.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3);
        
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        iko_eigen.eval(point,out);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),4.0e-1,"unscaled eigen 2");
      }

    }
    cout << endl;
    
  }

#endif    

#endif
  
  t.report();
  return 0;
}

