/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2022, Andrew W. Steiner
  
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
#include <o2scl/interp_krige.h>
#include <o2scl/test_mgr.h>
#include <o2scl/rng.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/deriv_gsl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef boost::numeric::ublas::vector<double> ubvector;

double f(double x, double mean, double sd) {
  return (sin(1.0/(0.3+x))-mean)/sd;
}

double df(double x, double sd) {
  return -cos(1.0/(0.3+x))/sd/pow(0.3+x,2.0);
}

double d2f(double x, double sd) {
  return 2.0*cos(1.0/(0.3+x))/sd/pow(0.3+x,3.0)-
    sin(1.0/(0.3+x))/sd/pow(0.3+x,4.0);
}

double covar(double x, double y) {
  return exp(-20.0*(x-y)*(x-y));
}

double covard(double x, double y) {
  return -exp(-20.0*(x-y)*(x-y))*(x-y)*40.0;
}

double covard2(double x, double y) {
  return -40.0*exp(-20.0*(x-y)*(x-y))+
    1600.0*exp(-20.0*(x-y)*(x-y))*(x-y)*(x-y);
}

double covari(double x, double a, double b) {
  double alpha=20.0;
  double ret=sqrt(o2scl_const::pi/alpha)/2.0*
    (gsl_sf_erf(sqrt(alpha)*(b-x))-
     gsl_sf_erf(sqrt(alpha)*(a-x)));
  return ret;
}

int main(void) {

  cout.setf(ios::scientific);

  deriv_gsl<> dg;
  inte_qag_gsl<> iqg;
  
  test_mgr t;
  t.set_output_level(2);

  // ---------------------------------------------------------------
  // Create test data

  static const size_t N=20;
  ubvector x(N), y(N);
  x[0]=0.0;
  y[0]=f(x[0],0.0,1.0);
  for(size_t i=1;i<N;i++) {
    x[i]=x[i-1]+pow(((double)i)/40.0,2.0);
    y[i]=f(x[i],0.0,1.0);
  }
  
  cout << "Data: " << endl;
  for(size_t i=0;i<N;i++) {
    cout.width(2);
    cout << x[i] << " ";
    cout.setf(ios::showpos);
    cout << y[i] << endl;
    cout.unsetf(ios::showpos);
  }
  cout << endl;

  double y_mean=vector_mean(y);
  double y_sd=vector_stddev(y,y_mean);
  for(size_t i=0;i<N;i++) {
    y[i]-=y_mean;
    y[i]/=y_sd;
  }

  cout << "Rescaled data: " << endl;
  for(size_t i=0;i<N;i++) {
    cout.width(2);
    cout << x[i] << " ";
    cout.setf(ios::showpos);
    cout << y[i] << endl;
    cout.unsetf(ios::showpos);
  }
  cout << endl;

  // ---------------------------------------------------------------
  // Test the exact derivatives

  if (true) {
    funct fptrx=std::bind(f,std::placeholders::_1,std::ref(y_mean),
                          std::ref(y_sd));
    for(size_t i=0;i<N;i++) {
      t.test_rel(dg.deriv(x[i],fptrx),
                 df(x[i],y_sd),1.0e-8,"deriv");
      t.test_rel(dg.deriv2(x[i],fptrx),
                 d2f(x[i],y_sd),4.0e-3,"deriv2");
    }
  }
  
  
  // ---------------------------------------------------------------
  // Test the interp_krige class

  interp_krige<ubvector> ik;
  std::function<double(double,double)> fp=covar;
  std::function<double(double,double)> fpd=covard;
  std::function<double(double,double)> fpd2=covard2;
  std::function<double(double,double,double)> fpi=covari;

  // ---------------------------------------------------------------
  // Test normal interpolation

  // Note that the interpolation performs poorly when the noise
  // term is zero because of the matrix inversion
  
  cout << "Normal interpolation:" << endl;
  ik.set_covar(N,x,y,fp);
  double xi0=ik.eval(x[0]);
  double xi1=ik.eval(x[N-1]);
  double xi2=ik.eval((x[0]+x[0])/2.0);
  t.test_rel(ik.eval(x[0]),y[0],0.2,"ik 1");
  t.test_rel(ik.eval(x[N-1]),y[N-1],0.2,"ik 2");
  t.test_rel(ik.eval((x[0]+x[1])/2.0),
             (y[0]+y[1])/2.0,0.2,"ik 3");
  cout << endl;

  // Show how extrapolation works
  cout << "Test extrapolation:" << endl;
  t.test_rel(ik.eval(10.0),0.0,1.0e-10,"ik 4");
  t.test_rel(ik.sigma(10.0),1.0,1.0e-10,"ik 5");
  cout << endl;
  
  // But now when the noise is non-zero, the interpolation is
  // of higher quality
  
  cout << "Normal interpolation with noise:" << endl;
  ik.set_covar_noise(N,x,y,fp,1.0e-9);
  double xn0=ik.eval(x[0]);
  double xn1=ik.eval(x[N-1]);
  double xn2=ik.eval((x[0]+x[0])/2.0);
  t.test_rel(ik.eval(x[0]),y[0],1.0e-4,"ik 6");
  t.test_rel(ik.eval(x[N-1]),y[N-1],1.0e-8,"ik 7");
  t.test_rel(ik.eval((x[0]+x[1])/2.0),
             (y[0]+y[1])/2.0,1.0e-4,"ik 8");
  t.test_abs(ik.sigma(x[0]),0.0,1.0e-6,"ik 9");
  t.test_abs(ik.sigma(x[N-1]),0.0,1.0e-7,"ik 10");
  cout << endl;

  // Test derivative -- this is extremely inaccurate
  cout << "Derivatives:" << endl;
  ik.set_covar_di_noise(N,x,y,fp,fpd,fpd2,fpi,1.0e-9);
  t.test_rel(ik.deriv(x[0]),df(x[0],y_sd),1.0,"ik 11");
  t.test_rel(ik.deriv(x[N-1]),df(x[N-1],y_sd),10.0,"ik 12");
  t.test_rel(ik.deriv((x[0]+x[1])/2.0),
             df((x[0]+x[1])/2.0,y_sd),1.0,"ik 13");
  cout << endl;
  
  // ---------------------------------------------------------------
  // Test normal interpolation with rescaling
  
  cout << "Normal interpolation with rescaling:" << endl;
  ik.set_covar_noise(N,x,y,fp,0.0,true);
  t.test_rel(ik.eval(x[0]),y[0],0.2,"ikr 1");
  t.test_rel(ik.eval(x[N-1]),y[N-1],0.2,"ikr 2");
  t.test_rel(ik.eval((x[0]+x[1])/2.0),(y[0]+y[1])/2.0,0.2,"ikr 3");
  t.test_rel(ik.eval(x[0]),xi0,5.0e-2,"ikr vs. ik 1");
  t.test_rel(ik.eval(x[N-1]),xi1,1.0e-12,"ikr vs. ik 2");
  t.test_rel(ik.eval((x[0]+x[1])/2.0),xi2,2.0e-3,"ikr vs. ik 3");
  cout << endl;

  // Just make sure this compiles
  prob_dens_gaussian pdg1=ik.gen_dist(1.0);
  prob_dens_gaussian pdg2=ik.gen_dist(1.5);
  
  // ---------------------------------------------------------------
  // Test interpolation with noise and rescaling
  
  cout << "Noisy interpolation with rescaling:" << endl;
  ik.set_covar_noise(N,x,y,fp,1.0e-9,true);
  t.test_rel(ik.eval(x[0]),y[0],0.2,"ikr 4");
  t.test_rel(ik.eval(x[N-1]),y[N-1],0.2,"ikr 5");
  t.test_rel(ik.eval((x[0]+x[1])/2.0),(y[0]+y[1])/2.0,0.2,"ikr 6");
  t.test_rel(ik.eval(x[0]),xn0,5.0e-7,"ikr vs. ik 4");
  t.test_rel(ik.eval(x[N-1]),xn1,5.0e-12,"ikr vs. ik 5");
  t.test_rel(ik.eval((x[0]+x[1])/2.0),xn2,1.0e-2,"ikr vs. ik 6");
  cout << endl;

  // ---------------------------------------------------------------
  // More precise data to test derivatives and integrals

  static const size_t N2=400;
  ubvector x2(N2), y2(N2);
  for(size_t i=0;i<N2;i++) {
    x2[i]=((double)i)*1.55/399.0;
    y2[i]=f(x2[i],0.0,1.0);
  }
  double y2_mean=vector_mean(y2);
  double y2_sd=vector_stddev(y2,y2_mean);
  for(size_t i=0;i<N2;i++) {
    y2[i]-=y2_mean;
    y2[i]/=y2_sd;
  }

  funct fptr=std::bind(f,std::placeholders::_1,std::ref(y2_mean),
                       std::ref(y2_sd));
  
  cout << "Test derivatives and integrals:" << endl;
  
  ik.set_covar_di_noise(N2,x2,y2,fp,fpd,fpd2,fpi,1.0e-12);
  t.test_rel(ik.deriv(x2[0]),df(x2[0],y2_sd),1.0e-3,"der 11");
  t.test_rel(ik.deriv(x2[N2-1]),df(x2[N2-1],y2_sd),5.0e-3,"der 12");
  t.test_rel(ik.deriv((x2[0]+x2[1])/2.0),
             df((x2[0]+x2[1])/2.0,y2_sd),1.0e-3,"der 13");
  
  t.test_rel(ik.deriv2(x2[0]),d2f(x2[0],y2_sd),0.05,"der2 11");
  t.test_rel(ik.deriv2(x2[N2-1]),d2f(x2[N2-1],y2_sd),0.05,"der2 12");
  t.test_rel(ik.deriv2((x2[0]+x2[1])/2.0),
             d2f((x2[0]+x2[1])/2.0,y2_sd),0.02,"der2 13");

  t.test_rel(iqg.integ(fptr,0.0,0.5),ik.integ(0.0,0.5),
             5.0e-2,"integ 1");
  t.test_rel(iqg.integ(fptr,0.5,1.0),ik.integ(0.5,1.0),
             5.0e-2,"integ 2");
  t.test_rel(iqg.integ(fptr,1.0,1.5),ik.integ(1.0,1.5),
             5.0e-2,"integ 3");
  cout << endl;

  cout << "Test derivatives and integrals with rescaling:" << endl;
  
  ik.set_covar_di_noise(N2,x2,y2,fp,fpd,fpd2,fpi,1.0e-12,true);
  t.test_rel(ik.deriv(x2[0]),df(x2[0],y2_sd),1.0e-3,"der 11");
  t.test_rel(ik.deriv(x2[N2-1]),df(x2[N2-1],y2_sd),5.0e-3,"der 12");
  t.test_rel(ik.deriv((x2[0]+x2[1])/2.0),
             df((x2[0]+x2[1])/2.0,y2_sd),1.0e-3,"der 13");
  
  t.test_rel(ik.deriv2(x2[0]),d2f(x2[0],y2_sd),0.05,"der2 11");
  t.test_rel(ik.deriv2(x2[N2-1]),d2f(x2[N2-1],y2_sd),0.05,"der2 12");
  t.test_rel(ik.deriv2((x2[0]+x2[1])/2.0),
             d2f((x2[0]+x2[1])/2.0,y2_sd),0.02,"der2 13");
  
  t.test_rel(iqg.integ(fptr,0.0,0.5),ik.integ(0.0,0.5),
             5.0e-2,"integ 1");
  t.test_rel(iqg.integ(fptr,0.5,1.0),ik.integ(0.5,1.0),
             5.0e-2,"integ 2");
  t.test_rel(iqg.integ(fptr,1.0,1.5),ik.integ(1.0,1.5),
             5.0e-2,"integ 3");
  cout << endl;
  
  // ---------------------------------------------------------------
  // Test interp_krige_optim interface

  interp_krige_optim<ubvector> iko;

  cout << "Class interp_krige_optim with simple interface." << endl;
  iko.verbose=1;
  iko.set(N,x,y);
  iko.verbose=0;

  t.test_rel(iko.eval(x[0]),y[0],1.0e-2,"iko 1");
  t.test_rel(iko.eval(x[N-1]),y[N-1],1.0e-4,"iko 2");
  t.test_rel(iko.eval((x[0]+x[1])/2.0),
             (y[0]+y[1])/2.0,1.0e-2,"iko 3");
  cout << endl;

  //cout << iko.deriv(1.5) << " " << cos(1.5) << endl;
  //cout << iko.deriv2(1.5) << " " << -sin(1.5) << endl;

  cout << "Class interp_krige_optim with simple interface, "
       << "rescaled version." << endl;
  
  iko.set(N,x,y,true);
  double len2=iko.get_length();

  t.test_rel(iko.eval(x[0]),y[0],1.0e-2,"ikor 1");
  t.test_rel(iko.eval(x[N-1]),y[N-1],1.0e-4,"ikor 2");
  t.test_rel(iko.eval((x[0]+x[1])/2.0),(y[0]+y[1])/2.0,1.0e-2,"ikor 3");
  cout << endl;

  // ---------------------------------------------------------------
  // Compare minimization functions

  // Start with the optimal length
  iko.verbose=1;
  iko.set(N,x,y,true);
  double len=iko.get_length();

  double mean_abs=0.0;
  for(size_t j=0;j<N;j++) {
    mean_abs+=fabs(y[j]);
  }
  mean_abs/=N;
  
  int success;
  iko.verbose=0;
  for(double ell=len/500.0;ell<len*30.01;ell*=pow(15000.0,0.01)) {
    cout << ell << " ";
    iko.mode=iko.mode_loo_cv;
    double q=iko.qual_fun(ell,mean_abs/1.0e8,success);
    cout.setf(ios::showpos);
    cout << q << " ";
    cout.unsetf(ios::showpos);
    cout << success << " ";
    iko.mode=iko.mode_max_lml;
    q=iko.qual_fun(ell,mean_abs/1.0e8,success);
    cout.setf(ios::showpos);
    cout << q << " ";
    cout.unsetf(ios::showpos);
    cout << success << " ";
    iko.mode=iko.mode_loo_cv_new;
    q=iko.qual_fun(ell,mean_abs/1.0e8,success);
    cout.setf(ios::showpos);
    cout << q << " ";
    cout.unsetf(ios::showpos);
    cout << success << endl;
  }
  cout << endl;

  // ---------------------------------------------------------------
  // Test interp_krige_optim interface with full minimization

  iko.full_min=true;
  iko.mode=iko.mode_loo_cv;
  
  cout << "Class interp_krige_optim with full minimization" << endl;
  iko.verbose=1;
  iko.set(N,x,y);
  iko.verbose=0;

  double exact=f(1.01,y_mean,y_sd);
  double res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-1,"iko 4");
  exact=f(1.0,y_mean,y_sd);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-1,"iko 5");

  iko.set(N,x,y,true);

  cout << "Class interp_krige_optim with full minimization "
       << "and rescaling" << endl;
  exact=f(1.01,y_mean,y_sd);
  res=iko.eval(1.01);
  t.test_rel(exact,res,1.0e-1,"iko 7");
  exact=f(1.0,y_mean,y_sd);
  res=iko.eval(1.0);
  t.test_rel(exact,res,1.0e-1,"iko 8");
  cout << endl;

  // ---------------------------------------------------------------
  // Test interp_krige_optim interface with full minimization
  // and refined data

  for(iko.mode=2;iko.mode<=3;iko.mode++) {

    iko.def_min.verbose=1;
    iko.verbose=1;
    cout << "Class interp_krige_optim with full minimization and "
         << "refined data." << endl;
    iko.set(N2,x2,y2);
    
    exact=f(1.01,y2_mean,y2_sd);
    res=iko.eval(1.01);
    t.test_rel(exact,res,1.0e-6,"iko 9");
    exact=f(1.0,y2_mean,y2_sd);
    res=iko.eval(1.0);
    t.test_rel(exact,res,1.0e-6,"iko 10");
    
    exact=df(1.01,y2_sd);
    res=iko.deriv(1.01);
    t.test_rel(exact,res,1.0e-6,"iko der 1");
    exact=df(1.0,y2_sd);
    res=iko.deriv(1.0);
    t.test_rel(exact,res,1.0e-6,"iko der 2");
    
    exact=d2f(1.01,y2_sd);
    res=iko.deriv2(1.01);
    t.test_rel(exact,res,1.0e-4,"iko der2 1");
    exact=d2f(1.0,y2_sd);
    res=iko.deriv2(1.0);
    t.test_rel(exact,res,1.0e-4,"iko der2 2");
    
    exact=iqg.integ(fptr,0.5,1.0);
    res=iko.integ(0.5,1.0);
    t.test_rel(exact,res,1.0e-5,"iko integ 1");
    exact=iqg.integ(fptr,1.0,1.5);
    res=iko.integ(1.0,1.5);
    t.test_rel(exact,res,1.0e-5,"iko integ 2");
    cout << endl;
    
    iko.set(N2,x2,y2,true);
    
    cout << "Class interp_krige_optim with full minimization "
         << "and rescaling" << endl;
    
    exact=f(1.01,y2_mean,y2_sd);
    res=iko.eval(1.01);
    t.test_rel(exact,res,1.0e-6,"iko 11");
    exact=f(1.0,y2_mean,y2_sd);
    res=iko.eval(1.0);
    t.test_rel(exact,res,1.0e-6,"iko 12");
    
    exact=df(1.01,y2_sd);
    res=iko.deriv(1.01);
    t.test_rel(exact,res,1.0e-6,"iko der 3");
    exact=df(1.0,y2_sd);
    res=iko.deriv(1.0);
    t.test_rel(exact,res,1.0e-6,"iko der 4");
    
    exact=d2f(1.01,y2_sd);
    res=iko.deriv2(1.01);
    t.test_rel(exact,res,1.0e-4,"iko der2 3");
    exact=d2f(1.0,y2_sd);
    res=iko.deriv2(1.0);
    t.test_rel(exact,res,1.0e-4,"iko der2 4");
    
    exact=iqg.integ(fptr,0.5,1.0);
    res=iko.integ(0.5,1.0);
    t.test_rel(exact,res,1.0e-5,"iko integ 3");
    exact=iqg.integ(fptr,1.0,1.5);
    res=iko.integ(1.0,1.5);
    t.test_rel(exact,res,1.0e-5,"iko integ 4");
    cout << endl;

  }
  
  t.report();

  return 0;
}
