/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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

#include <o2scl/eos_quark_cfl.h>
#include <o2scl/eos_quark_bag.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

typedef boost::numeric::ublas::vector<double> ubvector;

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  {
    // Various initializations

    double qq1, qq2, qq3, y[3], rr1, rr2, rr3, ss1, ss2, ss3;
    int i;

    eos_quark_njl nj;
    eos_quark_njl njt;
    eos_quark_cfl cfl;

    quark u(nj.up_default_mass,6.0);
    quark d(nj.down_default_mass,6.0);
    quark s(nj.strange_default_mass,6.0);
    quark u2(nj.up_default_mass,6.0);
    quark d2(nj.down_default_mass,6.0);
    quark s2(nj.strange_default_mass,6.0);
    quark u3(nj.up_default_mass,6.0);
    quark d3(nj.down_default_mass,6.0);
    quark s3(nj.strange_default_mass,6.0);
  
    mroot_hybrids<mm_funct> nd;
    inte_qng_gsl<funct> gl, gl2;
    thermo th, th2, th3;
  
    nd.tol_rel/=100.0;
    nd.tol_abs/=100.0;
  
    nj.set_quarks(u,d,s);
    nj.set_thermo(th);
    nj.set_parameters();

    njt.set_quarks(u2,d2,s2);
    njt.set_thermo(th2);
    njt.set_parameters();

    cfl.set_quarks(u3,d3,s3);
    cfl.set_thermo(th3);
    cfl.set_parameters_cfl();
    
    inte_qng_gsl<funct> ngnew;

    cfl.set_inte(ngnew);

    mm_funct fqq=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&eos_quark_njl::eos_quark_njl::gapfunqq),
       &nj,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);
    mm_funct fqq2=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&eos_quark_njl::eos_quark_njl::gapfunqq),
       &njt,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);
    mm_funct fqq3=std::bind
      (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
       (&eos_quark_njl::eos_quark_njl::gapfunqq),
       &cfl,std::placeholders::_1,std::placeholders::_2,
       std::placeholders::_3);

    //mm_funct_mfptr<eos_quark_njl> fqq(&nj,&eos_quark_njl::gapfunqq);
    //mm_funct_mfptr<eos_quark_njl> fqq2(&njt,&eos_quark_njl::gapfunqq);
    //   mm_funct_mfptr<eos_quark_njl> fqq3(&cfl,&eos_quark_njl::gapfunqq);
  
    // Set the quark chemical potentials

    u.mu=2.5;
    d.mu=2.5;
    s.mu=2.5;
    u2.mu=2.5;
    d2.mu=2.5;
    s2.mu=2.5;
  
    // Compute the masses from the quark condensates by 
    // solving the gap equation
  
    ubvector ax(3), ay(3);
    ax[0]=-1.0;
    ax[1]=-1.0;
    ax[2]=-2.0;
  
    nj.fromqq=true;
    nd.msolve(3,ax,fqq);
    cout << ax[0] << " " << ax[1] << " " << ax[2] << endl;

    // Compute the quark condensates from the masses, 
    // and show the solution works
  
    ax[0]=u.ms;
    ax[1]=d.ms;
    ax[2]=s.ms;
    nj.fromqq=false;
    nj.gapfunms(3,ax,ay);
    t.test_rel(ay[0],0.0,5.0e-9,"zero T ungapped gap eqn 1");
    t.test_rel(ay[1],0.0,5.0e-9,"zero T ungapped gap eqn 2");
    t.test_rel(ay[2],0.0,1.0e-6,"zero T ungapped gap eqn 3");

    // Now compute the full EOS at T=0
  
    nj.calc_eq_p(u,d,s,qq1,qq2,qq3,th);

    cout << "Comparing eos_quark_njl, eos_quark_njl, and eos_quark_cfl:" 
	 << endl;
    cout << "qq: " << u.qq << " " << d.qq << " " << s.qq << endl;
    cout << "sol: " << qq1 << " " << qq2 << " " << qq3 << endl;
    cout << "ms: " << u.ms << " " << d.ms << " " << s.ms << endl;
    cout << "n: " << u.n << " " << d.n << " " << s.n << endl;
    cout << "mu: " << u.mu << " " << d.mu << " " << s.mu << endl;
    cout << "ed: " << u.ed << " " << d.ed << " " << s.ed << endl;
    cout << "pr: " << u.pr << " " << d.pr << " " << s.pr << endl;
    cout << "th: " << th.ed << " " << th.pr << " " << th.en << endl;
    cout << endl;

    // Compute the masses from the quark condensates by 
    // solving the gap equation again for the 'njt' object

    ax[0]=-1.0;
    ax[1]=-1.0;
    ax[2]=-2.0;
    int vpy;
    njt.fromqq=true;
    nd.msolve(3,ax,fqq2);
    cout << ax[0] << " " << ax[1] << " " << ax[2] << endl;

    // Now compute the full EOS at T!=0

    u2.qq=ax[0];
    d2.qq=ax[1];
    s2.qq=ax[2];
    njt.fromqq=true;
    njt.calc_eq_temp_p(u2,d2,s2,rr1,rr2,rr3,th2,0.01);
  
    cout << "qq: " << u2.qq << " " << d2.qq << " " << s2.qq << endl;
    cout << "sol: " << rr1 << " " << rr2 << " " << rr3 << endl;
    cout << "ms: " << u2.ms << " " << d2.ms << " " << s2.ms << endl;
    cout << "n: " << u2.n << " " << d2.n << " " << s2.n << endl;
    cout << "mu: " << u2.mu << " " << d2.mu << " " << s2.mu << endl;
    cout << "ed: " << u2.ed << " " << d2.ed << " " << s2.ed << endl;
    cout << "pr: " << u2.pr << " " << d2.pr << " " << s2.pr << endl;
    cout << "th: " << th2.ed << " " << th2.pr << " " << th2.en << endl;
    cout << endl;

    // Compare T=0 and T!=0 results

    t.test_rel(u2.qq,u.qq,1.0e-6,"qqu");
    t.test_rel(d2.qq,d.qq,1.0e-6,"qqd");
    t.test_rel(s2.qq,s.qq,1.0e-6,"qqs");
    t.test_rel(rr1,0.0,1.0e-3,"qq1");
    t.test_rel(rr2,0.0,1.0e-3,"qq2");
    t.test_rel(rr3,0.0,1.0e-1,"qq3");
    t.test_rel(u2.n,u.n,1.0e-2,"nu");
    t.test_rel(d2.n,d.n,1.0e-2,"nd");
    t.test_rel(s2.n,s.n,1.0e-1,"ns");
    t.test_rel(u2.mu,u.mu,1.0e-6,"muu");
    t.test_rel(d2.mu,d.mu,1.0e-6,"mud");
    t.test_rel(s2.mu,s.mu,1.0e-6,"mus");
    t.test_rel(u2.ed,u.ed,1.0e-2,"edu");
    t.test_rel(d2.ed,d.ed,1.0e-2,"edd");
    t.test_rel(s2.ed,s.ed,1.0e-2,"eds");
    t.test_rel(u2.pr,u.pr,1.0e-4,"pru");
    t.test_rel(d2.pr,d.pr,1.0e-4,"prd");
    t.test_rel(s2.pr,s.pr,1.0e-4,"prs");
    t.test_rel(th.ed,th2.ed,1.0e-1,"thed");
    t.test_rel(th.pr,th2.pr,5.0e-3,"thpr");
    t.test_rel(th2.en,th.en,2.0e-1,"then");

    // Now compare the Delta=0 code with the NJL model at T!=0

    double gap1, gap2, gap3, n3, n8;
    quartic_real_coeff_cern<> cq;

    u3.mu=2.5;
    d3.mu=2.5;
    s3.mu=2.5;
    u3.del=0.0;
    d3.del=0.0;
    s3.del=0.0;

    cfl.gap_limit=-0.01;

    // Solve the gap equations
    ax[0]=-1.0;
    ax[1]=-1.0;
    ax[2]=-2.0;
    int vpz;
    nd.msolve(3,ax,fqq3);
    cout << ax[0] << " " << ax[1] << " " << ax[2] << endl;
  
    cfl.calc_eq_temp_p(u3,d3,s3,ss1,ss2,ss3,gap1,gap2,gap3,
		       0.0,0.0,n3,n8,th3,0.01);
  
    cout << "qq: " << u3.qq << " " << d3.qq << " " << s3.qq << endl;
    cout << "sol: " << ss1 << " " << ss2 << " " << ss3 << endl;
    cout << "ms: " << u3.ms << " " << d3.ms << " " << s3.ms << endl;
    cout << "n: " << u3.n << " " << d3.n << " " << s3.n << endl;
    cout << "mu: " << u3.mu << " " << d3.mu << " " << s3.mu << endl;
    cout << "ed: " << u3.ed << " " << d3.ed << " " << s3.ed << endl;
    cout << "pr: " << u3.pr << " " << d3.pr << " " << s3.pr << endl;
    cout << "th: " << th3.ed << " " << th3.pr << " " << th3.en << endl;

    t.test_rel(u2.qq,u3.qq,1.0e-6,"qqu");
    t.test_rel(d2.qq,d3.qq,1.0e-6,"qqd");
    t.test_rel(s2.qq,s3.qq,1.0e-6,"qqs");
    t.test_rel(ss1,0.0,1.0e-3,"qq1");
    t.test_rel(ss2,0.0,1.0e-3,"qq2");
    t.test_rel(ss3,0.0,1.0e-2,"qq3");
    t.test_rel(u2.n,u3.n,1.0e-2,"nu");
    t.test_rel(d2.n,d3.n,1.0e-2,"nd");
    t.test_rel(s2.n,s3.n,1.0e-1,"ns");
    t.test_rel(u2.mu,u3.mu,1.0e-6,"muu");
    t.test_rel(d2.mu,d3.mu,1.0e-6,"mud");
    t.test_rel(s2.mu,s3.mu,1.0e-6,"mus");
    t.test_rel(u2.ed,u3.ed,1.0e-2,"edu");
    t.test_rel(d2.ed,d3.ed,1.0e-2,"edd");
    t.test_rel(s2.ed,s3.ed,1.0e-2,"eds");
    t.test_rel(u2.pr,u3.pr,1.0e-4,"pru");
    t.test_rel(d2.pr,d3.pr,1.0e-4,"prd");
    t.test_rel(s2.pr,s3.pr,1.0e-4,"prs");
    t.test_rel(th3.ed,th2.ed,5.0e-2,"thed");
    t.test_rel(th3.pr,th2.pr,6.0e-4,"thpr");
    t.test_rel(th2.en,th3.en,1.0e-1,"then");

    t.test_rel(th3.pr,1.423674,1.0e-3,"make sure we're not changing.");
    //t.test_rel(th3.en,0.08216656,4.0e-3,"make sure we're not changing2.");

    //cout << "\nTest thd_potential(): " 
    //<< cfl.thd_potential(u3,d3,s3,0.0,0.0,0.01) << endl;
    //t.test_rel(-cfl.thd_potential(u3,d3,s3,0.0,0.0,0.01),
    //th3.pr,1.0e-6,"thd_potential:pr");
  
    double egv[36];

    double dedmuu[36]; 
    double dedmud[36]; 
    double dedmus[36]; 
    double dedmu[36]; 
    double dedmd[36]; 
    double dedms[36]; 
    double dedu[36]; 
    double dedd[36]; 
    double deds[36]; 
    double dedmu3[36];
    double dedmu8[36];

    // Test CFL phase
    u3.mu=2.5;
    d3.mu=2.5;
    s3.mu=2.5;
    u3.ms=0.0;
    d3.ms=0.0;
    s3.ms=0.0;
    u3.del=0.5;
    d3.del=0.5;
    s3.del=0.5;
    cfl.eigenvalues(1.0,0.0,0.0,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,dedms,
		    dedu,dedd,deds,dedmu3,dedmu8);
  
    egv[0]=sqrt(pow(1.0+2.5,2.0)+4.0*0.5*0.5);
    egv[1]=sqrt(pow(1.0-2.5,2.0)+4.0*0.5*0.5);
    egv[2]=-sqrt(pow(1.0-2.5,2.0)+4.0*0.5*0.5);
    egv[3]=-sqrt(pow(1.0+2.5,2.0)+4.0*0.5*0.5);
    for(i=0;i<8;i++) {
      egv[4+4*i]=sqrt(pow(1.0+2.5,2.0)+0.5*0.5);
      egv[5+4*i]=sqrt(pow(1.0-2.5,2.0)+0.5*0.5);
      egv[6+4*i]=-sqrt(pow(1.0-2.5,2.0)+0.5*0.5);
      egv[7+4*i]=-sqrt(pow(1.0+2.5,2.0)+0.5*0.5);
    }

    /*  
	gsl_sort_vector(ev);
	gsl_sort_vector(fv);
	cout << "\nTest eigenvalues of CFL phase: " << endl;
	for(i=0;i<36;i++) {
	t.test_rel(gsl_vector_get(ev,i),gsl_vector_get(fv,i),1.0e-6,"CFL");
	}
    */

    //Test 2SC phase 
    u3.mu=2.5;
    d3.mu=2.5;
    s3.mu=2.5;
    u3.ms=0.0;
    d3.ms=0.0;
    s3.ms=1.5;
    u3.del=0.0;
    d3.del=0.0;
    s3.del=0.5;
    cfl.eigenvalues(1.0,0.0,0.0,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,dedms,
		    dedu,dedd,deds,dedmu3,dedmu8);
  
    for(i=0;i<3;i++) {
      egv[0+4*i]=2.5-sqrt(1.0*1.0+1.5*1.5);
      egv[1+4*i]=2.5+sqrt(1.0*1.0+1.5*1.5);
      egv[2+4*i]=-2.5-sqrt(1.0*1.0+1.5*1.5);
      egv[3+4*i]=-2.5+sqrt(1.0*1.0+1.5*1.5);
    }
    for(i=0;i<2;i++) {
      egv[12+4*i]=-1.0+2.5;
      egv[13+4*i]=-1.0-2.5;
      egv[14+4*i]=1.0+2.5;
      egv[15+4*i]=1.0-2.5;
    }
    for(i=0;i<4;i++) {
      egv[20+4*i]=sqrt(pow(1.0+2.5,2.0)+0.5*0.5);
      egv[21+4*i]=sqrt(pow(1.0-2.5,2.0)+0.5*0.5);
      egv[22+4*i]=-sqrt(pow(1.0-2.5,2.0)+0.5*0.5);
      egv[23+4*i]=-sqrt(pow(1.0+2.5,2.0)+0.5*0.5);
    }

    /*
      cout << "\nTest eigenvalues of 2SC phase: " << endl;
      for(i=0;i<36;i++) {
      t.test_rel(gsl_vector_get(ev,i),gsl_vector_get(fv,i),1.0e-6,"2SC");
      }
    */

    cout << "\nCompare fixed mass and zero gaps with bag model:" << endl;
    eos_quark_bag bag;
    bag.bag_constant=0.0;

    u3.del=0.0;
    d3.del=0.0;
    s3.del=0.0;
    cfl.gap_limit=-1.0;
    cfl.fixed_mass=true;
    cfl.calc_eq_temp_p(u3,d3,s3,ss1,ss2,ss3,gap1,gap2,gap3,
		       0.0,0.0,n3,n8,th3,0.1);

    quark u4(u.m,6.0), d4(d.m,6.0), s4(s.m,6.0);
    u4.mu=u.mu;
    d4.mu=d.mu;
    s4.mu=s.mu;
    thermo th4;

    int ret=bag.calc_temp_p(u4,d4,s4,0.1,th4);
    t.test_gen(ret==0,"bag_temp success");
    t.test_rel(th4.ed-th3.ed,th3.pr-th4.pr,1.0e-2,"bag constant is constant");

    t.test_rel(u3.n,u4.n,1.0e-2,"nu");
    t.test_rel(d3.n,d4.n,1.0e-2,"nd");
    t.test_rel(s3.n,s4.n,1.0e-2,"ns");
    t.test_rel(th3.en,th4.en,4.0e-2,"entropy");

  }

  {
    // Inits
  
    double ss1,ss2,ss3,gap1,gap2,gap3,n3,n8;
    double ss12,ss22,ss32,gap12,gap22,gap32,n32,n82;
    eos_quark_cfl cfl2;
    thermo th, th2;
    quark u2(cfl2.up_default_mass,6.0);
    quark d2(cfl2.down_default_mass,6.0);
    quark s2(cfl2.strange_default_mass,6.0);
    cfl2.set_quarks(u2,d2,s2);
    cfl2.set_thermo(th2);

    // Run test functions
  
    cfl2.test_normal_eigenvalues(t);
    cfl2.test_gapped_eigenvalues(t);

    u2.mu=2.0;
    d2.mu=2.1;
    s2.mu=2.2;
    u2.ms=0.5;
    d2.ms=0.6;
    s2.ms=0.7;
    u2.del=1.0;
    d2.del=1.1;
    s2.del=1.2;
    cfl2.test_derivatives(1.5,0.0,0.0,t);

    cfl2.test_integration(t);

  }

  t.report();

  return 0;
}
  

