/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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

// Prevent warnings about global variables by creating a namespace
namespace eos_quark_njl_ts_ns {

  eos_quark_njl njx;
  eos_quark_njl njt;
  eos_quark_njl nj;

  quark u(njx.up_default_mass,6.0);
  quark d(njx.down_default_mass,6.0);
  quark s(njx.strange_default_mass,6.0);
  thermo th;
  int dtype;
  
  inte_qag_gsl<funct> gl;
  
  double temper=0.01;

  int ftsolve(size_t nv, const ubvector &x, ubvector &y) {
    double gap1, gap2, gap3;
    u.qq=x[0];
    d.qq=x[1];
    s.qq=x[2];
    njt.calc_eq_temp_p(u,d,s,gap1,gap2,gap3,th,temper);
    y[0]=gap1;
    y[1]=gap2;
    y[2]=gap3;
    return 0;
  }

  double omfun(double x);

  double omfun(double x) {
    double L=602.3/hc_mev_fm;
    double g1,g2,g3,G=1.835/L/L,K=12.36/pow(L,5.0);

    if (dtype==1) {
      u.qq=x;
      nj.fromqq=true;
      nj.calc_eq_p(u,d,s,g1,g2,g3,th);
      return -th.pr;
    } else if (dtype==2) {
      u.ms=x;
      nj.fromqq=false;
      nj.calc_eq_p(u,d,s,g1,g2,g3,th);
      th.pr+=2.0*G*(u.qq*u.qq+d.qq*d.qq+s.qq*s.qq)-
	4.0*K*u.qq*d.qq*s.qq;
      return -th.pr;
    }
    return 0.0;
  }

}

using namespace eos_quark_njl_ts_ns;

int main(void) {

  ubvector ax(3);
  test_mgr t;
  t.set_output_level(2);
  int ret=0;
  
  mroot_hybrids<mm_funct> nd;
  deriv_gsl<funct> df;
  int vpx=0;

  cout.setf(ios::scientific);

  nj.set_quarks(u,d,s);
  nj.set_thermo(th);
  t.test_gen(nj.set_parameters()==0,"set_parameters().");
  t.test_rel(nj.B0,21.6084,1.0e-4,"bag constant");
  
  u.mu=2.5; d.mu=2.5; s.mu=2.5;
  
  mm_funct fqq=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_quark_njl::gapfunqq),
     &nj,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  mm_funct fms=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_quark_njl::gapfunms),
     &nj,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  funct fderiv=omfun;
  
  cout << "Feynman-Hellman theorem" << endl;
  cout << "Verify that (partial Omega)/(Partial qqu) = 0" << endl;
  nj.set_quarks(u,d,s);
  nj.fromqq=true;
  ax[0]=-1.0; ax[1]=-1.0; ax[2]=-1.0; 
  nd.msolve(3,ax,fqq);
  dtype=1;
  df.h=0.01;
  ret=nd.msolve(3,ax,fqq);
  t.test_gen(ret==0,"successful solve.");
  t.test_rel(df.deriv(u.qq,fderiv),0.0,1.0e-6,"fh1");
  
  cout << "Verify that (partial (Omega-Omega_{vac}))/(Partial mu) = qqu" 
       << endl;
  nj.fromqq=false;
  ax[0]=0.2; ax[1]=0.2; ax[2]=2.0; 
  nd.msolve(3,ax,fms);
  cout << u.ms << " " << u.qq << " " << th.ed << " " << th.pr << endl;
  dtype=2;
  df.h=0.01;
  nd.msolve(3,ax,fms);
  t.test_rel(df.deriv(u.ms,fderiv),u.qq,1.0e-1,"fh2");

  // ---------------------------------------------------------
  // Finite temperature portion

  double gap1, gap2, gap3;
  ubvector axx(4);
  
  cout.setf(ios::scientific);

  nd.tol_rel/=100.0;
  nd.tol_abs/=100.0;

  nj.set_quarks(u,d,s);
  nj.set_thermo(th);
  nj.set_parameters();
  
  njt.set_quarks(u,d,s);
  njt.set_thermo(th);
  njt.set_parameters();
  
  mm_funct fqq2=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&eos_quark_njl::gapfunqq),
     &nj,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);
  mm_funct fts=ftsolve;

  u.mu=2.5;
  d.mu=2.5;
  s.mu=2.5;
  axx[0]=-1.0; 
  axx[1]=-1.0; 
  axx[2]=-1.0; 
  nj.fromqq=true;
  int r=nd.msolve(3,axx,fqq2);
  double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
  t1=axx[0];
  t2=axx[1];
  t3=axx[2];
  t4=u.n;
  t5=d.n;
  t6=s.n;
  t7=u.ed;
  t8=d.ed;
  t9=s.ed;
  t10=u.pr;
  t11=d.pr;
  t12=s.pr;
  t.test_rel((u.n*u.mu-u.pr-u.ed)/u.ed,0.0,1.0e-10,"therm. ident. u");
  t.test_rel((d.n*d.mu-d.pr-d.ed)/d.ed,0.0,1.0e-10,"therm. ident. d");
  t.test_rel((s.n*s.mu-s.pr-s.ed)/s.ed,0.0,1.0e-10,"therm. ident. s");

  // Check that fqq solves the gap equations for both the zero
  // and finite temperature code

  nj.calc_eq_p(u,d,s,gap1,gap2,gap3,th);
  
  t.test_rel(gap1,0.0,1.0e-6,"gap1");
  t.test_rel(gap2,0.0,1.0e-6,"gap2");
  t.test_rel(gap3,0.0,1.0e-6,"gap3");
  
  njt.calc_eq_temp_p(u,d,s,gap1,gap2,gap3,th,0.01);
  
  t.test_rel(gap1,0.0,2.0e-4,"gap1");
  t.test_rel(gap2,0.0,2.0e-4,"gap2");
  t.test_rel(gap3,0.0,2.0e-2,"gap3");
  
  cout << ": " << gap1 << " " << gap2 << " " << gap3 << endl;

  // Now solve the gap equations at finite temperature, and 
  // compare the quark condensates
  nd.msolve(3,axx,fts);
  t.test_rel(t1,axx[0],5.0e-4,"qqu");
  t.test_rel(t2,axx[1],5.0e-4,"qqd");
  t.test_rel(t3,axx[2],5.0e-3,"qqs");
  t.test_rel(t4,u.n,5.0e-4,"nu");
  t.test_rel(t5,d.n,5.0e-4,"nd");
  t.test_rel(t6,s.n,5.0e-3,"ns");
  t.test_rel(t7,u.ed,5.0e-4,"edu");
  t.test_rel(t8,d.ed,5.0e-4,"edd");
  t.test_rel(t9,s.ed,1.0e-3,"eds");
  t.test_rel(t10,u.pr,5.0e-4,"pru");
  t.test_rel(t11,d.pr,5.0e-4,"prd");
  t.test_rel(t12,s.pr,5.0e-4,"prs");

  t.report();

  return 0;
}

