/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2025, Andrew W. Steiner

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
#include <iostream>
#include <ctime>

#include <o2scl/invert_auto.h>
#include <o2scl/interpm_krige.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_cblas;
using namespace o2scl_linalg;

double ft(double x, double y) {
  return 3.0-2.0*x*x+7.0*y+5.0*x*x*y*y;
}

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

#ifdef O2SCL_SET_CUDA
#ifdef O2SCL_SET_ARMA

  // I don't want to enable this section of code by default in case
  // O2scl was included with CUDA support but the GPU was temporarily
  // unavailable.
  if (argc>=2 && ((string)argv[1])==((string)"benchmark")) {
    
    tensor2<> t1(10,10);
    std::vector<double> t2(100);
        
    for(size_t i=0;i<10;i++) {
      for(size_t j=0;j<10;j++) {
        if (i==j) {
          t1(i,j)=((double)(i+2));
          t2[i*10+j]=((double)(i+2));
        } else if (i>j) {
          t1(i,j)=1.0e-2*exp(-2.0*
                             pow(((double)i)+((double)j),2.0));
          t2[j*10+i]=1.0e-2*exp(-2.0*
                                pow(((double)i)+((double)j),2.0));
        } else {
          t1(i,j)=0.0;
          // Store the matrix in column-major order, as expected in
          // O2scl, but store it in the upper- rather than the
          // lower-triangular part
          t2[j*10+i]=1.0e-2*exp(-2.0*
                                pow(((double)i)+((double)j),2.0));
        }
      }
    }
    
    tensor2<> tx=t1;

    cout << "Original matrix:" << endl;
    matrix_out(cout,t1);
    cout << endl;

    cholesky_decomp(10,t1);
    cout << "Native decomposition" << endl;
    matrix_out(cout,t1);
    cout << endl;

    cholesky_decomp_two(10,tx);
    cout << "Native decomposition (v2)" << endl;
    matrix_out(cout,tx);
    cout << endl;

    cholesky_decomp_cuda(10,t2);
    vector<double> t2x(100);
    vector_copy(t2,t2x);
    tensor2<> t3(10,10);
    t3.swap_data(t2x);
    cout << "CUDA decomposition" << endl;
    matrix_out(cout,t3);
    cout << endl;

    matrix_invert_cholesky_auto micf;
    tensor2<> t1b(10,10), t3b(10,10);
    vector<double> t2b(100);
    
    micf.mode=micf.force_o2;
    cout << "Native inverse" << endl;
    micf.invert(10,t1,t1b);
    matrix_out(cout,t1b);
    cout << endl;

    matrix_invert_det_cholesky_cuda midcc;
    midcc.invert(10,t2,t2b);
    tensor2<> t2c(10,10);
    t2c.swap_data(t2b);
    cout << "CUDA inverse from midcc" << endl;
    matrix_out(cout,t2c);
    cout << endl;

    t.test_rel_mat(10,10,t1b,t2c,1.0e-6,"inverse native vs. cuda");
    
    micf.mode=micf.force_cuda;
    micf.invert(10,t3,t3b);
    cout << "CUDA inverse from micf" << endl;
    matrix_out(cout,t3b);
    cout << endl;

    t.test_rel_mat(10,10,t1b,t3b,1.0e-6,"inverse native vs. cuda v2");
    
  }

  if (argc>=2 && ((string)argv[1])==((string)"benchmark")) {
    
    // We choose a nearly diagonal positive symmetric matrix which
    // is easy to invert
    for(size_t n=10;n<10000;n*=2) {

      size_t mult=10;
    
      matrix_invert_cholesky_auto micf;
    
      cout << ((double)n) << " " << std::flush;
    
      for(size_t ell=0;ell<3;ell++) {
      
        if (ell==0) {
          micf.mode=micf.force_o2;
        } else if (ell==1) {
          micf.mode=micf.force_arma;
        } else {
          micf.mode=micf.force_cuda;
        }
      
        if (ell!=0 || n<=640) {
          tensor2<> t1[mult], t2[mult];
          vector<size_t> size={n,n};
        
          struct timespec ts0;
          timespec_get(&ts0,TIME_UTC);
        
          for(size_t k=0;k<mult;k++) {
            t1[k].resize(n,n);
            t2[k].resize(n,n);
            for(size_t i=0;i<n;i++) {
              for(size_t j=0;j<n;j++) {
                if (i==j) {
                  t1[k](i,j)=((double)(i+2));
                } else {
                  t1[k](i,j)=1.0e-2*
                    exp(-2.0*pow(((double)i)+((double)j),2.0));
                }
              }
            }
          }
        
          struct timespec ts1;
          timespec_get(&ts1,TIME_UTC);
        
          double diff0=(double)(ts1.tv_sec-ts0.tv_sec);
          double ndiff0=(double)(ts1.tv_nsec-ts0.tv_nsec);
          //cout << (diff0+ndiff0*1.0e-9)/((double)mult) << " " << std::flush;

          double det;
          for(size_t k=0;k<mult;k++) {
            if (k==mult-1) {
              micf.invert_det(n,t1[k],t2[k],det);
            } else {
              micf.invert(n,t1[k],t2[k]);
            }
          }
        
          struct timespec ts2;
          timespec_get(&ts2,TIME_UTC);
        
          double diff=(double)(ts2.tv_sec-ts1.tv_sec);
          double ndiff=(double)(ts2.tv_nsec-ts1.tv_nsec);
          cout << (diff+ndiff*1.0e-9)/((double)mult) << " " << std::flush;
        
          if (mult>1 && n==10) {
            std::cout << std::endl;
            tensor2<> t3;
            t3.resize(n,n);
            cout << "det: " << det << endl;
            dgemm(o2cblas_RowMajor,o2cblas_NoTrans,o2cblas_NoTrans,
                  n,n,n,1.0,t1[mult-1],t2[mult-1],0.0,t3);
            matrix_out(cout,t3);

            tensor2<> t4;
            t4.resize(n,n);
            matrix_set_identity(t4);

            t.test_abs_mat(10,10,t3,t4,1.0e-6,"inverse");
          }
        
        } else {
          cout << 0.0 << " " << std::flush;
        }
      }
      cout << endl;

    }

  }

  if (true) {
    typedef tensor2<> mat_x_t;
    typedef const const_matrix_row_gen<tensor2<>> mat_x_row_t;
    typedef tensor2<> mat_y_t;
    typedef const matrix_column_gen<tensor2<>> mat_y_col_t;
    
    vector<std::shared_ptr<mcovar_base<ubvector,mat_x_row_t>>> vmfrn;
    vmfrn.resize(1);
    std::shared_ptr<mcovar_funct_rbf_noise<
      ubvector,mat_x_row_t>> mfrn(new mcovar_funct_rbf_noise<ubvector,
                                  mat_x_row_t>);
    vmfrn[0]=mfrn;
    mfrn->len.resize(2);
    
    interpm_krige_optim
      <ubvector,mat_x_t,mat_x_row_t,mat_y_t,mat_y_col_t,
       tensor2<>,matrix_invert_cholesky_auto> iko_auto;

    static const size_t N=100;
    tensor2<> in, out;
    in.resize(N,2);
    out.resize(N,1);
    
    for(size_t i=0;i<N;i++) {
      double ix=((double)i);
      double x=3.0*sin(ix*ix);
      double y=5.0*cos(pow(ix,4.0));
      in.set(i,0,x);
      in.set(i,1,y);
      out.set(i,0,ft(x,y));
    }
    
    iko_auto.verbose=1;
    vector<double> len_list={0.3,0.7,0.8,0.9,0.95,
      1.0,1.25,1.5,2.0,3.0,7.0,10.0};
    vector<double> l10_list={-15,-13,-11,-9};
    vector<vector<double> > ptemp;
    ptemp.push_back(len_list);
    ptemp.push_back(len_list);
    ptemp.push_back(l10_list);
    vector<vector<vector<double>>> param_lists;
    param_lists.push_back(ptemp);
    
    iko_auto.set_covar(vmfrn,param_lists);
    iko_auto.rescale=true;
    iko_auto.set_data(2,1,N,in,out);
    cout << endl;
        
    gen_test_number<> gtn_x3;
    gtn_x3.set_radix(1.9);
    
    for(size_t j=0;j<20;j++) {
      ubvector point(2), pout(1);
      point[0]=gtn_x3.gen();
      point[1]=gtn_x3.gen();
      
      if (fabs(point[0])<3.0 && fabs(point[1])<5.0) {
        iko_auto.eval(point,pout);
        cout.setf(ios::showpos);
        cout << point[0] << " " << point[1] << " "
             << pout[0] << " " << ft(point[0],point[1]) << endl;
        cout.unsetf(ios::showpos);
        t.test_rel(pout[0],ft(point[0],point[1]),8.0,
                   "optim, rescaled, eigen, max_lml");
      }

    }
    cout << endl;
    
  }

#endif
#endif
  
  t.report();
  return 0;
}

