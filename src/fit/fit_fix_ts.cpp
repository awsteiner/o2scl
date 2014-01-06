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
#include <o2scl/fit_base.h>
#include <o2scl/fit_fix.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double fixv[3]={1.0,0.4,-0.4};

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

double fitfun(size_t n, const ubvector &p, double x) {
  return p[0]*p[0]*x+p[1]*cos(x)+p[2];
}

double fitfun0(size_t n, const ubvector &p, double x) {
  return fixv[0]*fixv[0]*x+p[0]*cos(x)+p[1];
}

double fitfun1(size_t n, const ubvector &p, double x) {
  return p[0]*p[0]*x+fixv[1]*cos(x)+p[1];
}

double fitfun2(size_t n, const ubvector &p, double x) {
  return p[0]*p[0]*x+p[1]*cos(x)+fixv[2];
}

int main(void) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(1);

  // Data
  size_t ndat=40;
  ubvector xdat(ndat), ydat(ndat), yerr(ndat);
  for(size_t i=0;i<ndat;i++) {
    xdat[i]=((double)i)/50.0;
    ydat[i]=sin(xdat[i]);
    yerr[i]=0.01;
  }

  fit_funct_fptr<> mf(fitfun);
  chi_fit_funct<> cff(40,xdat,ydat,yerr,mf);
  fit_funct_fptr<> mf0(fitfun0);
  chi_fit_funct<> cff0(40,xdat,ydat,yerr,mf0);
  fit_funct_fptr<> mf1(fitfun1);
  chi_fit_funct<> cff1(40,xdat,ydat,yerr,mf1);
  fit_funct_fptr<> mf2(fitfun2);
  chi_fit_funct<> cff2(40,xdat,ydat,yerr,mf2);
  int ret;

  fit_nonlin<> gf;
  fit_fix_pars<bool *> g;

  bool fixa[3]={false,false,false};
  bool *fixp=&fixa[0];

  // Normal fit_nonlin using normal function 
  cout << "Normal fit_nonlin using normal function:" << endl;
  ubvector par_full1(3);
  ubmatrix covar_full1(3,3);
  double chi2_full1;
  par_full1[0]=1.0;
  par_full1[1]=1.0;
  par_full1[2]=1.0;
  ret=gf.fit(3,par_full1,covar_full1,chi2_full1,cff);
  //cout << par_full1 << " " << chi2_full1 << endl;
  //cout << covar_full1 << endl;
  cout << endl;

  // Normal fit_nonlin fixing first parameter
  cout << "Normal fit_nonlin fixing first parameter:" << endl;
  ubvector par_11(2);
  ubmatrix covar_11(2,2);
  double chi2_11;
  matrix_set_identity(2,2,covar_11);
  par_11[0]=1.0;
  par_11[1]=1.0;
  ret=gf.fit(2,par_11,covar_11,chi2_11,cff0);
  //cout << par_11 << " " << chi2_11 << endl;
  //cout << covar_11 << endl;
  cout << endl;

  // Normal fit_nonlin fixing second parameter
  cout << "Normal fit_nonlin fixing second parameter:" << endl;
  ubvector par_21(2);
  ubmatrix covar_21(2,2);
  double chi2_21;
  matrix_set_identity(2,2,covar_21);
  par_21[0]=1.0;
  par_21[1]=1.0;
  ret=gf.fit(2,par_21,covar_21,chi2_21,cff1);
  //cout << par_21 << " " << chi2_21 << endl;
  //cout << covar_21 << endl;
  cout << endl;

  // Normal fit_nonlin fixing third parameter
  cout << "Normal fit_nonlin fixing third parameter:" << endl;
  ubvector par_31(2);
  ubmatrix covar_31(2,2);
  double chi2_31;
  matrix_set_identity(2,2,covar_31);
  par_31[0]=1.0;
  par_31[1]=1.0;
  ret=gf.fit(2,par_31,covar_31,chi2_31,cff2);
  //cout << par_31 << " " << chi2_31 << endl;
  //cout << covar_31 << endl;
  cout << endl;

  // Using fit_fix::fit() fitting all parameters
  cout << "Using fit_fix::fit() fitting all parameters:" << endl;
  ubvector par_full2(3);
  ubmatrix covar_full2(3,3);
  double chi2_full2;
  matrix_set_identity(3,3,covar_full2);
  par_full2[0]=1.0;
  par_full2[1]=1.0;
  par_full2[2]=1.0;
  ret=g.fit(3,par_full2,covar_full2,chi2_full2,cff);
  //cout << par_full2 << " " << chi2_full2 << endl;
  //cout << covar_full2 << endl;
  cout << endl;
  
  t.test_rel_arr(3,par_full1,par_full2,1.0e-9,"par_full1/2");
  t.test_rel_mat(3,3,covar_full1,covar_full2,1.0e-9,"covar_full1/2");
  t.test_rel(chi2_full1,chi2_full2,1.0e-9,"chi2_full1/2");

  // Using fit_fix::fit_fix() fitting all parameters
  cout << "Using fit_fix::fit_fix() fitting all parameters:" << endl;
  ubvector par_full3(3);
  ubmatrix covar_full3(3,3);
  double chi2_full3;
  matrix_set_identity(3,3,covar_full3);
  par_full3[0]=1.0;
  par_full3[1]=1.0;
  par_full3[2]=1.0;
  ret=g.fit_fix(3,par_full3,covar_full3,chi2_full3,cff,fixp);
  //cout << par_full3 << " " << chi2_full3 << endl;
  //cout << covar_full3 << endl;
  cout << endl;

  t.test_rel_arr(3,par_full1,par_full3,1.0e-9,"par_full1/3");
  t.test_rel_mat(3,3,covar_full1,covar_full3,1.0e-9,"covar_full1/3");
  t.test_rel(chi2_full1,chi2_full3,1.0e-9,"chi2_full1/3");

  // Using fit_fix and fixing first parameter
  cout << "Using fit_fix and fixing first parameter: " << endl;
  ubvector par_12(3);
  ubmatrix covar_12(2,2);
  double chi2_12;

  for(size_t k=0;k<2;k++) {

    g.expand_covar=(k!=0);
    if (k==1) covar_12.resize(3,3);

    par_12[0]=fixv[0];
    par_12[1]=1.0;
    par_12[2]=1.0;
    fixa[0]=true;
    ret=g.fit_fix(3,par_12,covar_12,chi2_12,cff,fixp);
    fixa[0]=false;
    //cout << par_12 << " " << chi2_12 << endl;
    //cout << covar_12 << endl;
    cout << endl;
    
    t.test_rel(par_11[0],par_12[1],1.0e-9,"par_11/12 1");
    t.test_rel(par_11[1],par_12[2],1.0e-9,"par_11/12 2");
    if (k==0) t.test_rel_mat(2,2,covar_11,covar_12,1.0e-9,"covar_11/12");
    t.test_rel(chi2_11,chi2_12,1.0e-9,"chi2_11/12");
  }

  // Using fit_fix and fixing second parameter
  cout << "Using fit_fix and fixing second parameter: " << endl;
  ubvector par_22(3);
  ubmatrix covar_22(2,2);
  double chi2_22;

  for(size_t k=0;k<2;k++) {

    g.expand_covar=(k!=0);
    if (k==1) covar_22.resize(3,3);

    par_22[0]=1.0;
    par_22[1]=fixv[1];
    par_22[2]=1.0;
    fixa[1]=true;
    ret=g.fit_fix(3,par_22,covar_22,chi2_22,cff,fixp);
    fixa[1]=false;
    //cout << par_22 << " " << chi2_22 << endl;
    //cout << covar_22 << endl;
    cout << endl;
    
    t.test_rel(par_21[0],par_22[0],1.0e-9,"par_21/22 1");
    t.test_rel(par_21[1],par_22[2],1.0e-9,"par_21/22 2");
    if (k==0) t.test_rel_mat(2,2,covar_21,covar_22,1.0e-9,"covar_21/22");
    t.test_rel(chi2_21,chi2_22,1.0e-9,"chi2_21/22");

  }

  // Using fit_fix and fixing third parameter
  cout << "Using fit_fix and fixing third parameter: " << endl;
  ubvector par_32(3);
  ubmatrix covar_32(2,2);
  double chi2_32;

  for(size_t k=0;k<2;k++) {

    g.expand_covar=(k!=0);
    if (k==1) covar_32.resize(3,3);
    
    par_32[0]=1.0;
    par_32[1]=1.0;
    par_32[2]=fixv[2];
    fixa[2]=true;
    ret=g.fit_fix(3,par_32,covar_32,chi2_32,cff,fixp);
    fixa[2]=false;
    //cout << par_32 << " " << chi2_32 << endl;
    //cout << covar_32 << endl;
    cout << endl;
    
    t.test_rel(par_31[0],par_32[0],1.0e-9,"par_31/32 1");
    t.test_rel(par_31[1],par_32[1],1.0e-9,"par_31/32 2");
    if (k==0) t.test_rel_mat(2,2,covar_31,covar_32,1.0e-9,"covar_31/32");
    t.test_rel(chi2_31,chi2_32,1.0e-9,"chi2_31/32");

  }

  t.report();
  return 0;
}
 
