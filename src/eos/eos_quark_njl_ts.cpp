/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_quark_njl.h>
#include <o2scl/test_mgr.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/inte_qag_gsl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  eos_quark_njl nj;

  quark &u=nj.def_up;
  quark &d=nj.def_down;
  quark &s=nj.def_strange;
  thermo &th=nj.def_thermo;
  
  mroot_hybrids<mm_funct> mh;
  deriv_gsl<funct> df;

  t.test_gen(nj.set_parameters()==0,"set_parameters().");
  t.test_rel(nj.B0,21.6084,1.0e-4,"bag constant");
  cout << endl;
  
  u.mu=2.5;
  d.mu=2.75;
  s.mu=3.0;
  
  cout << "Feynman-Hellman theorem:" << endl;
  cout << endl;
  cout << "Verify that (partial Omega)/(Partial qq) = 0" << endl;

  // First compute the EOS at the specified values of mu,
  // giving an initial guess for the quark condensates
  nj.from_qq=true;
  u.qq=-1.0;
  d.qq=-1.0;
  s.qq=-1.0;
  int ret=nj.calc_p(u,d,s,th);
  t.test_gen(ret==0,"EOS success.");

  // Check that calc_p() is solving the gap equations
  double gap1, gap2, gap3;
  nj.calc_eq_p(u,d,s,gap1,gap2,gap3,th);
  
  t.test_rel(gap1,0.0,1.0e-9,"gap1a");
  t.test_rel(gap2,0.0,1.0e-9,"gap2a");
  t.test_rel(gap3,0.0,1.0e-9,"gap3a");

  // Compare with the finite-temperature code and make sure
  // the results match
  nj.calc_p(u,d,s,th);
  thermo zt=th;
  
  nj.calc_temp_p(u,d,s,0.1/hc_mev_fm,th);
  thermo ft=th;
  t.test_rel(zt.ed,ft.ed,1.0e-5,"finite and zero T ed");
  t.test_rel(zt.pr,ft.pr,1.0e-5,"finite and zero T pr");
  t.test_abs(zt.en,ft.en,5.0e-2,"finite and zero T en");
  
  // Then construct the functor for the derivative with
  // respect to the up quark condensate
  funct fderiv_u=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,bool)>
     (&eos_quark_njl::f_therm_pot),
     &nj,std::placeholders::_1,d.qq,s.qq,u.ms,d.ms,s.ms,true);
  
  df.h=0.01;
  double der_u=df.deriv(u.qq,fderiv_u);
  t.test_rel(der_u,0.0,1.0e-8,"fh_u");

  // Test the derivative wrt the down quark condensate
  nj.calc_p(u,d,s,th);
  funct fderiv_d=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,bool)>
     (&eos_quark_njl::f_therm_pot),
     &nj,u.qq,std::placeholders::_1,s.qq,u.ms,d.ms,s.ms,true);
  
  double der_d=df.deriv(d.qq,fderiv_d);
  t.test_rel(der_d,0.0,1.0e-9,"fh_d");
  
  // Test the derivative wrt the strange quark condensate
  nj.calc_p(u,d,s,th);
  funct fderiv_s=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,bool)>
     (&eos_quark_njl::f_therm_pot),
     &nj,u.qq,d.qq,std::placeholders::_1,u.ms,d.ms,s.ms,true);
  
  double der_s=df.deriv(s.qq,fderiv_s);
  t.test_rel(der_s,0.0,1.0e-9,"fh_s");
  cout << endl;
    
  cout << "Verify that (partial (Omega-Omega_{vac}))/(Partial m) = qq" 
       << endl;
  nj.from_qq=false;

  u.ms=0.2;
  d.ms=0.3;
  s.ms=2.0;
  ret=nj.calc_p(u,d,s,th);
  t.test_gen(ret==0,"success 2.");
  // Save the up quark condensate for later comparison
  double qqu=u.qq;
  
  funct fderiv2_u=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,bool)>
     (&eos_quark_njl::f_therm_pot),
     &nj,u.qq,d.qq,s.qq,std::placeholders::_1,d.ms,s.ms,false);

  der_u=df.deriv(u.ms,fderiv2_u);
  t.test_rel(der_u,qqu,1.0e-8,"fh2_u");

  u.ms=0.2;
  d.ms=0.3;
  s.ms=2.0;
  ret=nj.calc_p(u,d,s,th);
  // Save the down quark condensate for later comparison
  double qqd=d.qq;
  
  funct fderiv2_d=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,bool)>
     (&eos_quark_njl::f_therm_pot),
     &nj,u.qq,d.qq,s.qq,u.ms,std::placeholders::_1,s.ms,false);

  der_d=df.deriv(d.ms,fderiv2_d);
  t.test_rel(der_d,qqd,1.0e-8,"fh2_d");

  u.ms=0.2;
  d.ms=0.3;
  s.ms=2.0;
  ret=nj.calc_p(u,d,s,th);
  // Save the strange quark condensate for later comparison
  double qqs=s.qq;
  
  funct fderiv2_s=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,bool)>
     (&eos_quark_njl::f_therm_pot),
     &nj,u.qq,d.qq,s.qq,u.ms,d.ms,std::placeholders::_1,false);
  
  der_s=df.deriv(s.ms,fderiv2_s);
  t.test_rel(der_s,qqs,1.0e-8,"fh2_s");
  cout << endl;
  
  // ---------------------------------------------------------
  // Finite temperature portion
  
  cout << "Feynman-Hellman theorem at T>0:" << endl;
  cout << endl;
  cout << "Verify that (partial Omega)/(Partial qq) = 0" << endl;
  
  // First compute the EOS at the specified values of mu,
  // giving an initial guess for the quark condensates
  nj.from_qq=true;
  u.qq=-1.0;
  d.qq=-1.0;
  s.qq=-1.0;
  double T=0.1;
  ret=nj.calc_temp_p(u,d,s,T,th);
  t.test_gen(ret==0,"EOS success 3.");

  // Check that calc_temp_p() is solving the gap equations
  nj.calc_eq_temp_p(u,d,s,gap1,gap2,gap3,th,T);
  
  t.test_rel(gap1,0.0,1.0e-9,"gap1b");
  t.test_rel(gap2,0.0,1.0e-9,"gap2b");
  t.test_rel(gap3,0.0,1.0e-9,"gap3b");
  
  // Then construct the functor for the derivative with
  // respect to the up quark condensate
  funct fderivT_u=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,bool)>
     (&eos_quark_njl::f_therm_pot_T),
     &nj,std::placeholders::_1,d.qq,s.qq,u.ms,d.ms,s.ms,T,true);
  
  der_u=df.deriv(u.qq,fderivT_u);
  t.test_rel(der_u,0.0,1.0e-9,"fhT_u");

  // Test the derivative wrt the down quark condensate
  nj.calc_temp_p(u,d,s,T,th);
  funct fderivT_d=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,bool)>
     (&eos_quark_njl::f_therm_pot_T),
     &nj,u.qq,std::placeholders::_1,s.qq,u.ms,d.ms,s.ms,T,true);
  
  der_d=df.deriv(d.qq,fderivT_d);
  t.test_rel(der_d,0.0,1.0e-9,"fhT_d");
  
  // Test the derivative wrt the strange quark condensate
  nj.calc_temp_p(u,d,s,T,th);
  funct fderivT_s=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,bool)>
     (&eos_quark_njl::f_therm_pot_T),
     &nj,u.qq,d.qq,std::placeholders::_1,u.ms,d.ms,s.ms,T,true);
  
  der_s=df.deriv(s.qq,fderivT_s);
  t.test_rel(der_s,0.0,1.0e-9,"fhT_s");
  cout << endl;
    
  cout << "Verify that (partial (Omega-Omega_{vac}))/(Partial m) = qq" 
       << endl;
  nj.from_qq=false;

  u.ms=0.2;
  d.ms=0.3;
  s.ms=2.0;
  ret=nj.calc_temp_p(u,d,s,T,th);
  t.test_gen(ret==0,"success 4.");
  // Save the up quark condensate for later comparison
  qqu=u.qq;
  
  funct fderivT2_u=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,bool)>
     (&eos_quark_njl::f_therm_pot_T),
     &nj,u.qq,d.qq,s.qq,std::placeholders::_1,d.ms,s.ms,T,false);

  der_u=df.deriv(u.ms,fderivT2_u);
  t.test_rel(der_u,qqu,1.0e-8,"fh2_u");

  u.ms=0.2;
  d.ms=0.3;
  s.ms=2.0;
  ret=nj.calc_temp_p(u,d,s,T,th);
  // Save the down quark condensate for later comparison
  qqd=d.qq;
  
  funct fderivT2_d=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,bool)>
     (&eos_quark_njl::f_therm_pot_T),
     &nj,u.qq,d.qq,s.qq,u.ms,std::placeholders::_1,s.ms,T,false);

  der_d=df.deriv(d.ms,fderivT2_d);
  t.test_rel(der_d,qqd,1.0e-8,"fh2_d");

  u.ms=0.2;
  d.ms=0.3;
  s.ms=2.0;
  ret=nj.calc_temp_p(u,d,s,T,th);
  // Save the strange quark condensate for later comparison
  qqs=s.qq;
  
  funct fderivT2_s=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,bool)>
     (&eos_quark_njl::f_therm_pot_T),
     &nj,u.qq,d.qq,s.qq,u.ms,d.ms,std::placeholders::_1,T,false);
  
  der_s=df.deriv(s.ms,fderivT2_s);
  t.test_rel(der_s,qqs,1.0e-8,"fh2_s");
  cout << endl;
    
  eos_quark_njl_vec njv;

  t.test_gen(njv.set_parameters()==0,"set_parameters().");
  cout << njv.B0 << endl;
  t.test_rel(njv.B0,21.6084,1.0e-4,"bag constant");
  cout << endl;

  nj.from_qq=false;
  u.mu=2.5;
  d.mu=2.75;
  s.mu=3.0;
  u.nu=2.5;
  d.nu=2.75;
  s.nu=3.0;
  u.ms=0.2;
  d.ms=0.3;
  s.ms=2.0;
  ret=nj.calc_p(u,d,s,th);
    
  njv.from_qq=false;
  njv.GV=njv.G;
  u.mu=2.5;
  d.mu=2.75;
  s.mu=3.0;
  u.nu=2.5;
  d.nu=2.75;
  s.nu=3.0;
  u.ms=0.2;
  d.ms=0.3;
  s.ms=2.0;

  ret=njv.calc_p(u,d,s,th);
  zt=th;
  cout << "T=0: " << th.ed << " " << th.pr << " "
       << u.mu << " " << u.n << " " << d.mu << " "
       << u.mu*u.n+d.mu*d.n+s.mu*s.n << endl;
  cout << u.qq << " " << u.ms << endl;
  
  njv.calc_temp_p(u,d,s,0.1/hc_mev_fm,th);
  ft=th;
  cout << "T>0: " << th.ed << " " << th.pr << " " 
       << u.mu << " " << u.n << " " << d.mu << " "
       << u.mu*u.n+d.mu*d.n+s.mu*s.n << endl;
  t.test_rel(zt.ed,ft.ed,1.0e-4,"T=0 vs T>0");
  t.test_rel(zt.pr,ft.pr,1.0e-4,"T=0 vs T>0");
  cout << endl;
  
  ret=njv.calc_p(u,d,s,th);
    
  // Then construct the functor for the derivative with
  // respect to the up quark condensate
  njv.from_qq=true;
    
  funct fderiv_u2=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,double,double,bool)>
     (&eos_quark_njl_vec::f_therm_pot),
     &njv,std::placeholders::_1,d.qq,s.qq,u.ms,d.ms,s.ms,
     u.nu,d.nu,s.nu,true);

  df.h=0.01;
  der_u=df.deriv(u.qq,fderiv_u2);
  cout << "der_u: " << der_u << endl;
  t.test_rel(der_u,0.0,5.0e-9,"fh_u");

  // Then construct the functor for the derivative with
  // respect to the down quark condensate
  funct fderiv_d2=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,double,double,bool)>
     (&eos_quark_njl_vec::f_therm_pot),
     &njv,u.qq,std::placeholders::_1,s.qq,u.ms,d.ms,s.ms,
     u.nu,d.nu,s.nu,true);

  df.h=0.01;
  der_d=df.deriv(d.qq,fderiv_d2);
  cout << "der_d: " << der_d << endl;
  t.test_rel(der_d,0.0,4.0e-6,"fh_d");

  // Then construct the functor for the derivative with
  // respect to the strange quark condensate
  funct fderiv_s2=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,double,double,bool)>
     (&eos_quark_njl_vec::f_therm_pot),
     &njv,u.qq,d.qq,std::placeholders::_1,u.ms,d.ms,s.ms,
     u.nu,d.nu,s.nu,true);

  df.h=0.01;
  der_s=df.deriv(s.qq,fderiv_s2);
  cout << "der_s: " << der_s << endl;
  t.test_rel(der_s,0.0,4.0e-6,"fh_s");

  funct fderiv_u3=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,double,double,bool)>
     (&eos_quark_njl_vec::f_therm_pot),
     &njv,u.qq,d.qq,s.qq,u.ms,d.ms,s.ms,
     std::placeholders::_1,d.nu,s.nu,true);
    
  df.h=0.01;
  double der_u3=df.deriv(u.nu,fderiv_u3);
  cout << "Here: " << der_u3 << endl;
  t.test_rel(der_u3,0.0,5.0e-6,"fh_u");

  funct fderiv_d3=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,double,double,bool)>
     (&eos_quark_njl_vec::f_therm_pot),
     &njv,u.qq,d.qq,s.qq,u.ms,d.ms,s.ms,
     u.nu,std::placeholders::_1,s.nu,true);
    
  df.h=0.01;
  double der_d3=df.deriv(d.nu,fderiv_d3);
  cout << "Here: " << der_d3 << endl;
  t.test_rel(der_d3,0.0,5.0e-6,"fh_d");

  funct fderiv_s3=std::bind
    (std::mem_fn<double(double,double,double,double,double,
                        double,double,double,double,bool)>
     (&eos_quark_njl_vec::f_therm_pot),
     &njv,u.qq,d.qq,s.qq,u.ms,d.ms,s.ms,
     u.nu,d.nu,std::placeholders::_1,true);
    
  df.h=0.01;
  double der_s3=df.deriv(s.nu,fderiv_s3);
  cout << "Here: " << der_s3 << endl;
  t.test_rel(der_s3,0.0,5.0e-5,"fh_s");

  t.report();

  return 0;
}

