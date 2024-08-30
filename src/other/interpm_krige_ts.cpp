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
#include <o2scl/interpm_krige.h>
#include <o2scl/interp_krige.h>
#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

#ifdef O2SCL_ARMA
#include <armadillo>
#endif
#ifdef O2SCL_SET_EIGEN
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
  return 3.0-2.0*x*x+7.0*y+5.0*x*x*y*y;
}

double ftdx(double x, double y) {
  return -4.0*x+10.0*x*y*y;
}

double ftdy(double x, double y) {
  return 7.0+10.0*x*x*y;
}

double ftdx2(double x, double y) {
  return -4.0+10.0*y*y;
}

double ftdy2(double x, double y) {
  return 10.0*x*x;
}

double ftdxdy(double x, double y) {
  return 20.0*x*y;
}

void generate_table(table<> &tab, size_t N=100) {

  tab.clear();
  
  tab.line_of_names("x y z");
  
  for(size_t i=0;i<N;i++) {
    double ix=((double)i);
    double x=3.0*sin(ix*ix);
    double y=5.0*cos(pow(ix,4.0));
    vector<double> line={x,y,ft(x,y)};
    tab.line_of_data(3,line);
  }
  
  return;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);
  
  vector<string> col_list_x={"x","y"};
  vector<string> col_list_y={"z"};
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef o2scl::const_matrix_view_table<> mat_x_t;
  typedef const const_matrix_row_gen<mat_x_t> mat_x_row_t;
  typedef o2scl::matrix_view_table<> mat_y_t;
  typedef const matrix_column_gen<mat_y_t> mat_y_col_t;

  if (true) {

    // Test the quad_correl covariance function
    mcovar_funct_quad_correl<ubvector,mat_x_row_t> m;
    
    m.len.resize(3);
    m.pos.resize(3);
    m.slope.resize(3);
    m.len[0]=1.0;
    m.len[1]=1.1;
    m.len[2]=1.3;
    m.pos[0]=1.7;
    m.pos[1]=1.9;
    m.pos[2]=2.3;
    m.slope[0]=2.9;
    m.slope[1]=3.1;
    m.slope[2]=3.7;
    m.log10_noise=-15;

    ubvector x(3), y(3);
    x[0]=4.1;
    x[1]=4.3;
    x[2]=5.3;
    y[0]=5.9;
    y[1]=6.1;
    y[2]=6.7;

    double m0=m.covar(x,y);
    x[0]+=1.0e-4;
    double m1=m.covar(x,y);
    x[0]-=1.0e-4;
    double d01=m.deriv_tl(x,y,0);
    cout << d01 << " " << (m1-m0)/1.0e-4 << endl;
    t.test_rel(d01,(m1-m0)/1.0e-4,1.0e-4,"quad_correl deriv 1");
    
    x[1]+=1.0e-4;
    double m2=m.covar(x,y);
    x[1]-=1.0e-4;
    double d02=m.deriv_tl(x,y,1);
    cout << d02 << " " << (m2-m0)/1.0e-4 << endl;
    t.test_rel(d02,(m2-m0)/1.0e-4,1.0e-4,"quad_correl deriv 2");
    
    x[2]+=1.0e-4;
    double m3=m.covar(x,y);
    x[2]-=1.0e-4;
    double d03=m.deriv_tl(x,y,2);
    cout << d03 << " " << (m3-m0)/1.0e-4 << endl;
    t.test_rel(d03,(m3-m0)/1.0e-4,1.0e-4,"quad_correl deriv 3");

    double td0=m.deriv2_tl(x,y,0,1);
    cout << td0 << " " << 1.1759e-5 << endl;
    // Computed using mathematica
    t.test_rel(td0,1.1759e-5,1.0e-4,
               "quad_correl mixed second deriv");
    
    double td1=m.deriv2_tl(x,y,2,2);
    // Computed using mathematica
    cout << td1 << " " << -1.26167e-3 << endl;
    t.test_rel(td1,-1.26167e-3,1.0e-4,
               "quad_correl second deriv");
    
  }
  
  if (true) {

    cout << "--------------------------------------------" << endl;
    cout << "interpm_krige_optim, unscaled, loo_cv\n" << endl;
    
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(2);
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,mat_y_t,mat_y_col_t,ubmatrix> iko;
    iko.mode=iko.mode_loo_cv;

    table<> tab3;
    generate_table(tab3);
    
    const_matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    iko.verbose=1;
    vector<double> len_list={0.3,0.7,0.8,0.9,0.95,
      1.0,1.25,1.5,2.0,3.0,7.0,10.0};
    vector<double> l10_list={-15,-13,-11,-9};
    vector<vector<double>> ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);
    
    iko.set_covar(vmfrn,param_lists);

    iko.rescale=false;
    iko.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3);
    cout << endl;
        
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        iko.eval(point,out);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),1.0e+1,
                   "optim, unscaled, loo_cv");
      }

    }
    cout << endl;
  }

  if (true) {

    cout << "--------------------------------------------" << endl;
    cout << "interpm_krige_optim, rescaled, max_lml\n" << endl;
  
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(2);
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,mat_y_t,mat_y_col_t,ubmatrix> iko;
    iko.mode=iko.mode_max_lml;

    table<> tab3;
    generate_table(tab3);
    
    hdf_file hf;
    hf.open_or_create("interpm_krige_ts_data2.o2");
    hdf_output(hf,tab3,"tab");
    hf.close();
    
    const_matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    vector<double> len_list={0.3,0.7,0.8,0.9,0.95,
      1.0,1.25,1.5,2.0,3.0,7.0,10.0};
    vector<double> l10_list={-15,-13,-11,-9};
    vector<vector<double> > ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);
    
    iko.set_covar(vmfrn,param_lists);
    iko.rescale=true;
    iko.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3);
    cout << endl;
        
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        iko.eval(point,out);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),6.0,
                   "optim, rescaled, max_lml");
        
        iko.deriv(point,out,0);
        cout << "  " << out[0] << " " << ftdx(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdx(point[0],point[1]),10.0,
                   "optim, rescaled, max_lml, dfdx");
        
        iko.deriv(point,out,1);
        cout << out[0] << " " << ftdy(point[0],point[1]) << endl;
        t.test_rel(out[0],ftdy(point[0],point[1]),1.0e1,
                   "optim, rescaled, max_lml, dfdy");
        
        iko.deriv2(point,out,0,0);
        cout << "  " << out[0] << " " << ftdx2(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdx2(point[0],point[1]),10.0,
                   "optim, rescaled, max_lml, dfdx2");
        
        iko.deriv2(point,out,0,1);
        cout << out[0] << " " << ftdxdy(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdxdy(point[0],point[1]),100.0,
                   "optim, rescaled, max_lml, dfdxdy");
        
        iko.deriv2(point,out,1,1);
        cout << out[0] << " " << ftdy2(point[0],point[1]) << endl;
        t.test_rel(out[0],ftdy2(point[0],point[1]),1.0e3,
                   "optim, rescaled, max_lml, dfdy2");
      }

    }
    cout << endl;

  }

  if (true) {

    cout << "--------------------------------------------" << endl;
    cout << "interpm_krige_optim, rescaled, max_lml, full min.\n" << endl;
  
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(2);
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,mat_y_t,mat_y_col_t,ubmatrix> iko;
    iko.mode=iko.mode_max_lml;
    iko.full_min=true;

    table<> tab3;
    generate_table(tab3);
    
    hdf_file hf;
    hf.open_or_create("interpm_krige_ts_data2.o2");
    hdf_output(hf,tab3,"tab");
    hf.close();
    
    const_matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    vector<double> len_list={0.3,0.7,0.8,0.9,0.95,
      1.0,1.25,1.5,2.0,3.0,7.0,10.0};
    vector<double> l10_list={-15,-13,-11,-9};
    vector<vector<double> > ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);
    
    iko.set_covar(vmfrn,param_lists);
    iko.rescale=true;
    iko.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3);
        
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        iko.eval(point,out);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),6.0,
                   "optim, rescaled, max_lml");
        
        iko.deriv(point,out,0);
        cout << "  " << out[0] << " " << ftdx(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdx(point[0],point[1]),10.0,
                   "optim, rescaled, max_lml, dfdx");
        
        iko.deriv(point,out,1);
        cout << out[0] << " " << ftdy(point[0],point[1]) << endl;
        t.test_rel(out[0],ftdy(point[0],point[1]),1.0e1,
                   "optim, rescaled, max_lml, dfdy");
        
        iko.deriv2(point,out,0,0);
        cout << "  " << out[0] << " " << ftdx2(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdx2(point[0],point[1]),10.0,
                   "optim, rescaled, max_lml, dfdx2");
        
        iko.deriv2(point,out,0,1);
        cout << out[0] << " " << ftdxdy(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdxdy(point[0],point[1]),100.0,
                   "optim, rescaled, max_lml, dfdxdy");
        
        iko.deriv2(point,out,1,1);
        cout << out[0] << " " << ftdy2(point[0],point[1]) << endl;
        t.test_rel(out[0],ftdy2(point[0],point[1]),1.0e3,
                   "optim, rescaled, max_lml, dfdy2");
      }

    }
    cout << endl;
    
  }

  if (true) {

    cout << "--------------------------------------------" << endl;
    cout << "interpm_krige_optim, rescaled, loo_cv_bf\n" << endl;
  
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(2);
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,mat_y_t,mat_y_col_t,ubmatrix> iko;
    iko.mode=iko.mode_loo_cv_bf;

    table<> tab3;
    generate_table(tab3);
    
    hdf_file hf;
    hf.open_or_create("interpm_krige_ts_data2.o2");
    hdf_output(hf,tab3,"tab");
    hf.close();
    
    const_matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    vector<double> len_list={0.9,0.95,1.0,1.25,1.5,2.0,3.0,7.0};
    vector<double> l10_list={-15,-13,-11};
    vector<vector<double> > ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);
    
    iko.set_covar(vmfrn,param_lists);
    iko.rescale=true;
    iko.verbose=2;
    iko.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3);
    cout << endl;
        
    for(size_t j=0;j<20;j++) {
      ubvector point(2), out(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();

      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        
        iko.eval(point,out);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << out[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(out[0],ft(point[0],point[1]),6.0,
                   "optim, rescaled, max_lml");
        
        iko.deriv(point,out,0);
        cout << "  " << out[0] << " " << ftdx(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdx(point[0],point[1]),10.0,
                   "optim, rescaled, max_lml, dfdx");
        
        iko.deriv(point,out,1);
        cout << out[0] << " " << ftdy(point[0],point[1]) << endl;
        t.test_rel(out[0],ftdy(point[0],point[1]),1.0e1,
                   "optim, rescaled, max_lml, dfdy");
        
        iko.deriv2(point,out,0,0);
        cout << "  " << out[0] << " " << ftdx2(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdx2(point[0],point[1]),10.0,
                   "optim, rescaled, max_lml, dfdx2");
        
        iko.deriv2(point,out,0,1);
        cout << out[0] << " " << ftdxdy(point[0],point[1]) << " ";
        t.test_rel(out[0],ftdxdy(point[0],point[1]),100.0,
                   "optim, rescaled, max_lml, dfdxdy");
        
        iko.deriv2(point,out,1,1);
        cout << out[0] << " " << ftdy2(point[0],point[1]) << endl;
        t.test_rel(out[0],ftdy2(point[0],point[1]),1.0e3,
                   "optim, rescaled, max_lml, dfdy2");
      }

    }
    cout << endl;

  }

  if (true) {

    cout << "--------------------------------------------" << endl;
    cout << "Compare qual_fun() for interpm_krige_optim" << endl;
    cout << " and interp_krige_optim (unscaled):\n" << endl;
    
    static const size_t N=20;
    ubvector x(N), y(N);
    x[0]=0.0;
    y[0]=f(x[0],0.0,1.0);
    for(size_t i=1;i<N;i++) {
      x[i]=x[i-1]+pow(((double)i)/40.0,2.0);
      y[i]=f(x[i],0.0,1.0);
    }

    table<> tab4;
    tab4.line_of_names("x y");
    for(size_t i=0;i<N;i++) {
      vector<double> line={x[i],y[i]};
      tab4.line_of_data(2,line);
    }

    vector<double> len_list={0.01,0.03,0.1,0.3,1.0};
    vector<double> l10_list={-15,-13,-11,-9};
    vector<vector<double> > ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);
    
    const_matrix_view_table<> mvt_x4(tab4,{"x"});
    matrix_view_table<> mvt_y4(tab4,{"y"});
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_col_t,ubmatrix> iko;
    
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(1);
    
    iko.set_covar(vmfrn,param_lists);
    iko.verbose=2;
    iko.rescale=false;
    iko.set_data(1,1,tab4.get_nlines(),mvt_x4,mvt_y4);
    
    interp_krige_optim<ubvector,ubvector,covar_funct_rbf_noise> iko2;
    
    covar_funct_rbf_noise cfrn;
    
    iko2.set_covar_optim(cfrn,ptemp);
    iko2.set(N,x,y);

    ubvector p;
    p.resize(3);
    p[0]=0.1;
    p[1]=0.1;
    p[2]=1.0e-8;
    mfrn->set_params(p);
    cfrn.set_params(p);

    int success;

    t.set_output_level(2);
    iko2.mode=iko2.mode_loo_cv;
    iko.mode=iko.mode_loo_cv;
    iko.verbose=3;
    iko2.verbose=3;
    t.test_rel(iko.qual_fun(0,success),
               iko2.qual_fun(success),1.0e-10,
               "optim, compare 1d and multid, unscaled, loo_cv.");

    iko2.mode=iko2.mode_max_lml;
    iko.mode=iko.mode_max_lml;
    t.test_rel(iko.qual_fun(0,success),
               iko2.qual_fun(success),1.0e-10,
               "optim, compare multid and 1d, unscaled, max_lml.");
    t.set_output_level(1);
    cout << endl;
    
  }

  if (true) {

    cout << "--------------------------------------------" << endl;
    cout << "Compare qual_fun() for interpm_krige_optim" << endl;
    cout << " and interp_krige_optim (rescaled):\n" << endl;
    
    static const size_t N=20;
    ubvector x(N), y(N);
    x[0]=0.0;
    y[0]=f(x[0],0.0,1.0);
    for(size_t i=1;i<N;i++) {
      x[i]=x[i-1]+pow(((double)i)/40.0,2.0);
      y[i]=f(x[i],0.0,1.0);
    }

    table<> tab4;
    tab4.line_of_names("x y");
    for(size_t i=0;i<N;i++) {
      vector<double> line={x[i],y[i]};
      tab4.line_of_data(2,line);
    }

    vector<double> len_list={0.01,0.03,0.1,0.3,1.0};
    vector<double> l10_list={-15,-13,-11,-9};
    vector<vector<double> > ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);

    const_matrix_view_table<> mvt_x4(tab4,{"x"});
    matrix_view_table<> mvt_y4(tab4,{"y"});
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_col_t,ubmatrix> iko;
    
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(1);

    iko.set_covar(vmfrn,param_lists);
    iko.rescale=true;
    iko.set_data(1,1,tab4.get_nlines(),mvt_x4,mvt_y4);
    
    interp_krige_optim<ubvector,ubvector,covar_funct_rbf_noise> iko2;
    
    covar_funct_rbf_noise cfrn;
    
    iko2.set(N,x,y,cfrn,ptemp,true);

    ubvector p;
    p.resize(3);
    p[0]=0.1;
    p[1]=0.1;
    p[2]=1.0e-8;
    mfrn->set_params(p);
    cfrn.set_params(p);

    int success;
    
    t.set_output_level(2);
    iko2.mode=iko2.mode_loo_cv;
    iko.mode=iko.mode_loo_cv;
    t.test_rel(iko.qual_fun(0,success),
               iko2.qual_fun(success),1.0e-10,
               "compare multid and 1d, rescaled, loo_cv");

    iko2.mode=iko2.mode_max_lml;
    iko.mode=iko.mode_max_lml;
    t.test_rel(iko.qual_fun(0,success),
               iko2.qual_fun(success),1.0e-10,
               "compare multid and 1d, rescaled, max_lml.");
    t.set_output_level(1);
    cout << endl;
  }
  
#ifdef O2SCL_SET_EIGEN  

  {

    cout << "--------------------------------------------" << endl;
    cout << "interpm_krige_optim, eigen, rescaled, max_lml\n" << endl;
  
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(2);
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_col_t,Eigen::MatrixXd,
       matrix_invert_det_eigen<> > iko_eigen;
    
    iko_eigen.mode=iko_eigen.mode_max_lml;

    table<> tab3;
    generate_table(tab3);
    
    const_matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    iko_eigen.verbose=1;
    vector<double> len_list={0.3,0.7,0.8,0.9,0.95,
      1.0,1.25,1.5,2.0,3.0,7.0,10.0};
    vector<double> l10_list={-15,-13,-11,-9};
    vector<vector<double> > ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);
    
    iko_eigen.set_covar(vmfrn,param_lists);
    iko_eigen.rescale=true;
    iko_eigen.set_data(2,1,tab3.get_nlines(),mvt_x3,mvt_y3);
    cout << endl;
        
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
        t.test_rel(out[0],ft(point[0],point[1]),8.0,
                   "optim, rescaled, eigen, max_lml");
      }

    }
    cout << endl;

  }
  
#endif

#ifdef O2SCL_NEVER_DEFINED

#ifdef O2SCL_ARMA  
  
  {
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,
       mat_y_t,mat_y_col_t,arma::mat,
       matrix_invert_det_sympd_arma<> > iko_arma;

    table<> tab3;
    generate_table(tab3);
    
    matrix_view_table<> mvt_x3(tab3,col_list_x);
    matrix_view_table<> mvt_y3(tab3,col_list_y);

    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    iko_arma.verbose=1;
    iko_arma.nlen=50;
    iko.rescale=false;
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

#endif
  
  t.report();
  return 0;
}

