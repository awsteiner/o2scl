/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2017, Andrew W. Steiner

  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------
*/
#include <string>
#include <ctime>
#include <o2scl/test_mgr.h>
#include <o2scl/poly.h>
#include <o2scl/misc.h>

/*
  This program evaluates the performance and accuracy of the O2scl
  polynomial solvers.
*/

using namespace std;
using namespace o2scl;

int wid=21;

class stats {
public:
  std::string type;
  std::vector<double> w;
  double s1, s2;
  double m1, m2;
  size_t nok;
  size_t nbad;
  size_t ntests;
  double time;

  int out() {
    cout << type << endl;
    cout << "   Average quality: " 
	 << s1 << " " << s2 << " " << m1 << " " << m2 << endl;
    cout << "   Worst case coefficients: " << endl;
    cout << "   " << w[0] << " " << w[1] << " " << w[2] << " " 
	 << w[3] << " " << w[4] << " " << w[5] << endl;
    cout << "   Nok: " << nok << " Nbad: " << nbad 
	 << " Ntot: " << ntests << endl;
    cout << "   Time: " << time << endl;
    return 0;
  }

};

void test_quadratic_real_coeff(quadratic_real_coeff *po, 
			       string str, double alpha, double e1, 
			       double e2, double e3, double e4,
			       double &time_taken);
void test_quadratic_complex(quadratic_complex *po, string str, 
			    double e1, double e2, double e3, double e4);
void test_cubic_real_coeff(cubic_real_coeff *po, string str, 
			   double alpha, double e1, double e2, double e3, 
			   double e4);
void test_cubic_complex(cubic_complex *po, string str, double e1, 
			double e2, double e3, double e4);
void test_quartic_real(quartic_real *po, string str, double alpha, 
		       double e1, double e2, double e3, double e4);
void test_quartic_real_coeff(quartic_real_coeff *po, string str, 
			     double e1, double e2, double e3, double e4);
void test_quartic_complex(quartic_complex *po, string str, 
			  double e1, double e2, double e3, double e4);

void test_quadratic_real_coeff(quadratic_real_coeff *po, string str, 
			       double alpha, stats &s) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> cr2,cr3,czo2,czo3,cap,cbp,ccp,cdp,cr1,czo1;
  double ca,cb,cc;
  complex<double> i(0.0,1.0);
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  lt1=clock();
  size_t ne=0, nok=0, nbad=0;

  gen_test_number<100> ga, gb, gc;
  for(int j1=0;j1<100;j1++) {
    ca=ga.gen();
    for(int j2=0;j2<100;j2++) {
      cb=gb.gen()*alpha;
      for(int j3=0;j3<100;j3++) {
	cc=gc.gen();

	// Ensure there is a solution
	if (fabs(ca)>0.0) {

	  po->solve_rc(ca,cb,cc,cr1,cr2);
	  
	  cbp=-(cr1+cr2)*ca;
	  ccp=(cr1*cr2)*ca;
	  
	  czo1=(ca*cr1+cb)*cr1+cc;
	  czo2=(ca*cr2+cb)*cr2+cc;
	  q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
		  fabs(cc-ccp.real())*fabs(cc-ccp.real()));
	  q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2));
	  
	  s1+=q1;
	  if (q1>m1) {
	    m1=q1;
	  }
	  s2+=q2;
	  if (q2>m2) {
	    m2=q2;
	    s.w[0]=ca;
	    s.w[1]=cb;
	    s.w[2]=cc;
	  }
	  ne++;
	  
	  if (fabs(q1)<1.0e-8 && fabs(q2)<1.0e-8) nok++;
	  else nbad++;

	}
      }
    }
  }

  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);

  s.time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  s.s1=s1;
  s.s2=s2;
  s.m1=m1;
  s.m2=m2;
  s.type=str;
  s.ntests=ne;
  s.nok=nok;
  s.nbad=nbad;

  return;
}

void test_quadratic_complex(quadratic_complex *po, string str, 
			    stats &s) {

  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> ca,cb,cc,cd,cr1,cr2,cr3,czo1,czo2,czo3,cap,cbp,ccp,cdp;
  complex<double> i(0.0,1.0);
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  size_t ne=0, nok=0, nbad=0;
  gen_test_number<15> ga, gb, gc, gd, ge, gf;
  double rca, rcb, rcc, rcd, rce, rcf;

  lt1=clock();
  
  for(int j1=0;j1<15;j1++) {
    rca=ga.gen();
    for(int j2=0;j2<15;j2++) {
      rcb=gb.gen();
      for(int j3=0;j3<15;j3++) {
	rcc=gc.gen();
	for(int j4=0;j4<15;j4++) {
	  rcd=gd.gen();
	  for(int j5=0;j5<15;j5++) {
	    rce=ge.gen();
	    for(int j6=0;j6<15;j6++) {
	      rcf=gf.gen();
	      
	      if (fabs(rca)>0.0 || fabs(rcb>0.0)) {
		ca=rca+i*rcb;
		cb=rcc+i*rcd;
		cc=rce+i*rcf;
    
		po->solve_c(ca,cb,cc,cr1,cr2);
	      
		cbp=-(cr1+cr2)*ca;
		ccp=(cr1*cr2)*ca;
	      
		czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
		czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
		q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
			abs(cc-ccp)*abs(cc-ccp));
		q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2));
		s1+=q1;
		if (q1>m1) m1=q1;
		s2+=q2;
		if (q2>m2) {
		  m2=q2;
		  s.w[0]=rca;
		  s.w[1]=rcb;
		  s.w[2]=rcc;
		  s.w[3]=rcd;
		  s.w[4]=rce;
		  s.w[5]=rcf;
		}
		ne++;
		if (fabs(q1)<1.0e-8 && fabs(q2)<1.0e-8) nok++;
		else nbad++;
		
	      }

	    }
	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);

  s.time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  s.s1=s1;
  s.s2=s2;
  s.m1=m1;
  s.m2=m2;
  s.type=str;
  s.ntests=ne;
  s.nok=nok;
  s.nbad=nbad;

  return;
}

void test_cubic_real_coeff(cubic_real_coeff *po, string str, 
			   double alpha, stats &s) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> cr2,cr3,czo2,czo3,cap,cbp,ccp,cdp;
  double ca,cb,cc,cd,cr1,czo1;
  complex<double> i(0.0,1.0);
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  size_t ne=0, nok=0, nbad=0;
  gen_test_number<21> ga, gb, gc, gd;

  lt1=clock();

  for(int j1=0;j1<21;j1++) {
    ca=ga.gen()*alpha;
    for(int j2=0;j2<21;j2++) {
      cb=gb.gen();
      for(int j3=0;j3<21;j3++) {
	cc=gc.gen()*alpha;
	for(int j4=0;j4<21;j4++) {
	  cd=gd.gen();
	  
	  if (fabs(ca)>0.0) {
	    po->solve_rc(ca,cb,cc,cd,cr1,cr2,cr3);
	    
	    cbp=-(cr1+cr2+cr3)*ca;
	    ccp=(cr1*cr2+cr1*cr3+cr2*cr3)*ca;
	    cdp=-(cr1*cr2*cr3)*ca;
	    
	    czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
	    czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
	    czo3=((ca*cr3+cb)*cr3+cc)*cr3+cd;
	    q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
		    fabs(cc-ccp.real())*fabs(cc-ccp.real())+
		    fabs(cd-cdp.real())*fabs(cd-cdp.real()));
	    q2=sqrt(fabs(czo1)*fabs(czo1)+abs(czo2)*abs(czo2)+
		    abs(czo3)*abs(czo3));
	    s1+=q1;
	    if (q1>m1) m1=q1;
	    s2+=q2;
	    if (q2>m2) {
	      m2=q2;
	      s.w[0]=ca;
	      s.w[1]=cb;
	      s.w[2]=cc;
	      s.w[3]=cd;
	    }
	    if (q2>m2) m2=q2;
	    ne++;
	    if (fabs(q1)<1.0e-8 && fabs(q2)<1.0e-8) nok++;
	    else nbad++;
	    
	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);

  s.time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  s.s1=s1;
  s.s2=s2;
  s.m1=m1;
  s.m2=m2;
  s.type=str;
  s.ntests=ne;
  s.nok=nok;
  s.nbad=nbad;

  return;
}

void test_cubic_complex(cubic_complex *po, string str, stats &s) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> ca,cb,cc,cd,cr1,cr2,cr3,czo1,czo2,czo3,cap,cbp,ccp,cdp;
  complex<double> i(0.0,1.0);
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  size_t ne=0, nok=0, nbad=0;
  gen_test_number<15> ga, gb, gc, gd;
  double rca, rcb, rcc, rcd;

  lt1=clock();
  for(int it=0;it<2;it++) {
    for(int j1=0;j1<15;j1++) {
      rca=ga.gen();
      for(int j2=0;j2<15;j2++) {
	rcb=gb.gen();
	for(int j3=0;j3<15;j3++) {
	  rcc=gc.gen();
	  for(int j4=0;j4<15;j4++) {
	    rcd=gd.gen();

	    if (it==0) {
	      ca=rca+i;
	      cb=rcb+i;
	      cc=rcc+i;
	      cd=rcd+i;
	    } else {
	      ca=1.0+i*rca;
	      cb=1.0+i*rcb;
	      cc=1.0+i*rcc;
	      cd=1.0+i*rcd;
	    }
    
	    cbp=-(cr1+cr2+cr3)*ca;
	    ccp=(cr1*cr2+cr1*cr3+cr2*cr3)*ca;
	    cdp=-(cr1*cr2*cr3)*ca;
      
	    czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
	    czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
	    czo3=((ca*cr3+cb)*cr3+cc)*cr3+cd;
	    q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
		    abs(cc-ccp)*abs(cc-ccp)+
		    abs(cd-cdp)*abs(cd-cdp));
	    q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
		    abs(czo3)*abs(czo3));
	    s1+=q1;
	    if (q1>m1) m1=q1;
	    s2+=q2;
	    if (q2>m2) {
	      m2=q2;
	      s.w[0]=rca;
	      s.w[1]=rcb;
	      s.w[2]=rcc;
	      s.w[3]=rcd;
	      s.w[4]=it;
	    }
	    ne++;
	    if (fabs(q1)<1.0e-8 && fabs(q2)<1.0e-8) nok++;
	    else nbad++;

	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);

  s.time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  s.s1=s1;
  s.s2=s2;
  s.m1=m1;
  s.m2=m2;
  s.type=str;
  s.ntests=ne;
  s.nok=nok;
  s.nbad=nbad;

  return;
}

void test_quartic_real(quartic_real *po, string str, double alpha, stats &s) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  double cr1,cr2,cr3,cr4;
  double r1,r2,r3,r4;
  double ca,cb,cc,cd,ce,zo1,zo2,zo3,zo4;
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  size_t ne=0, nok=0, nbad=0;
  gen_test_number<21> ga, gb, gc, gd;

  lt1=clock();

  for(int j1=0;j1<21;j1++) {
    r1=ga.gen();
    for(int j2=0;j2<21;j2++) {
      r2=-r1+alpha*gb.gen();
      for(int j3=0;j3<21;j3++) {
	r3=gc.gen();
	for(int j4=0;j4<21;j4++) {
	  r4=-r3+alpha*gd.gen();
	  
	  ca=1.0;
	  cb=-(r1+r2+r3+r4);
	  cc=(r1*r2+r1*r3+r2*r3+r1*r4+r2*r4+r3*r4);
	  cd=-(r1*r2*r3+r1*r2*r4+r1*r3*r4+r2*r3*r4);
	  ce=r1*r2*r3*r4;
	  
	  po->solve_r(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);
	  
	  zo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
	  zo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
	  zo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
	  zo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
	  q1=sqrt(zo1*zo1+zo2*zo2+zo3*zo3+zo4*zo4);
	  q2=0.0;
	  s1+=q1;
	  s2+=q2;
	  if (q1>m1) {
	    m1=q1;
	    s.w[0]=r1;
	    s.w[1]=r2;
	    s.w[2]=r3;
	    s.w[3]=r4;
	  }
	  ne++;
	  if (fabs(q1)<1.0e-8 && fabs(q2)<1.0e-8) nok++;
	  else nbad++;
	}
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);

  s.time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  s.s1=s1;
  s.s2=s2;
  s.m1=m1;
  s.m2=m2;
  s.type=str;
  s.ntests=ne;
  s.nok=nok;
  s.nbad=nbad;

  return;
}

void test_quartic_real_coeff(quartic_real_coeff *po, string str, 
			     stats &s) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> cr1,cr2,cr3,cr4,czo1,czo2,czo3,czo4;  
  complex<double> i(0.0,1.0),cap,cbp,ccp,cdp,cep;
  double ca,cb,cc,cd,ce;
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  size_t ne=0, nok=0, nbad=0;
  gen_test_number<15> ga, gb, gc, gd, ge;

  lt1=clock();

  for(int j1=0;j1<15;j1++) {
    ca=ga.gen();
    for(int j2=0;j2<15;j2++) {
      cb=gb.gen();
      for(int j3=0;j3<15;j3++) {
	cc=gc.gen();
	for(int j4=0;j4<15;j4++) {
	  cd=gd.gen();
	  for(int j5=0;j5<15;j5++) {
	    ce=ge.gen();
	    
	    if (fabs(ca)>0.0) {
	      po->solve_rc(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);
	      
	      cbp=-(cr1+cr2+cr3+cr4)*ca;
	      ccp=(cr1*cr2+cr1*cr3+cr2*cr3+cr1*cr4+cr2*cr4+cr3*cr4)*ca;
	      cdp=-(cr1*cr2*cr3+cr1*cr2*cr4+cr1*cr3*cr4+cr2*cr3*cr4)*ca;
	      cep=cr1*cr2*cr3*cr4*ca;
	      
	      czo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
	      czo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
	      czo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
	      czo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
	      q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
		      abs(cc-ccp)*abs(cc-ccp)+
		      abs(cd-cdp)*abs(cd-cdp)+
		      abs(ce-cep)*abs(ce-cep));
	      q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
		      abs(czo3)*abs(czo3)+abs(czo4)*abs(czo4));
	      s1+=q1;
	      if (q1>m1) m1=q1;
	      s2+=q2;
	      if (q2>m2) {
		s.w[0]=ca;
		s.w[1]=cb;
		s.w[2]=cc;
		s.w[3]=cd;
		s.w[4]=ce;
		m2=q2;
	      }
	      ne++;
	      if (fabs(q1)<1.0e-8 && fabs(q2)<1.0e-8) nok++;
	      else nbad++;
	    }

	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);

  s.time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  s.s1=s1;
  s.s2=s2;
  s.m1=m1;
  s.m2=m2;
  s.type=str;
  s.ntests=ne;
  s.nok=nok;
  s.nbad=nbad;

  return;
}

void test_quartic_complex(quartic_complex *po, string str, 
			  stats &s) {
  double s1,s2,m1,m2;
  clock_t lt1, lt2;
  complex<double> ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4,czo1,czo2,czo3,czo4;  
  complex<double> i(0.0,1.0),cap,cbp,ccp,cdp,cep;
  size_t j;
  double q1,q2;
  s1=0.0;
  s2=0.0;
  m1=0.0;
  m2=0.0;
  size_t ne=0, nok=0, nbad=0;
  gen_test_number<15> ga, gb, gc, gd, ge;
  double rca, rcb, rcc, rcd, rce;
  lt1=clock();

  for(int it=0;it<2;it++) {
    for(int j1=0;j1<15;j1++) {
      rca=ga.gen();
      for(int j2=0;j2<15;j2++) {
	rcb=gb.gen();
	for(int j3=0;j3<15;j3++) {
	  rcc=gc.gen();
	  for(int j4=0;j4<15;j4++) {
	    rcd=gd.gen();
	    for(int j5=0;j5<15;j5++) {
	      rce=ge.gen();

	      if (it==0) {
		ca=rca+i;
		cb=rcb+i;
		cc=rcc+i;
		cd=rcd+i;
		ce=rce+i;
	      } else {
		ca=1.0+i*rca;
		cb=1.0+i*rcb;
		cc=1.0+i*rcc;
		cd=1.0+i*rcd;
		ce=1.0+i*rce;
	      }

	      po->solve_c(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);
	    
	      cbp=-(cr1+cr2+cr3+cr4)*ca;
	      ccp=(cr1*cr2+cr1*cr3+cr2*cr3+cr1*cr4+cr2*cr4+cr3*cr4)*ca;
	      cdp=-(cr1*cr2*cr3+cr1*cr2*cr4+cr1*cr3*cr4+cr2*cr3*cr4)*ca;
	      cep=cr1*cr2*cr3*cr4*ca;
	    
	      czo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
	      czo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
	      czo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
	      czo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
	      q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
		      abs(cc-ccp)*abs(cc-ccp)+
		      abs(cd-cdp)*abs(cd-cdp)+
		      abs(ce-cep)*abs(ce-cep));
	      q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
		      abs(czo3)*abs(czo3)+abs(czo4)*abs(czo4));
	      if (!finite(q1)) {
		cout << cr1 << " " << cr2 << " " << cr3 << " " << cr4 << endl;
		exit(-1);
	      }
	      s1+=q1;
	      if (q1>m1) m1=q1;
	      s2+=q2;
	      if (q2>m2) {
		m2=q2;
		s.w[0]=rca;
		s.w[1]=rcb;
		s.w[2]=rcc;
		s.w[3]=rcd;
		s.w[4]=rce;
		s.w[5]=it;
	      }
	      ne++;
	      if (fabs(q1)<1.0e-8 && fabs(q2)<1.0e-8) nok++;
	      else nbad++;
	    }

	  }
	}
      }
    }
  }
  lt2=clock();
  s1/=((double)ne);
  s2/=((double)ne);

  s.time=((double)(lt2-lt1))/CLOCKS_PER_SEC;
  s.s1=s1;
  s.s2=s2;
  s.m1=m1;
  s.m2=m2;
  s.type=str;
  s.ntests=ne;
  s.nok=nok;
  s.nbad=nbad;

  return;
}

int main(void) {
  
  cout.setf(ios::left | ios::scientific);
  cout.precision(4);

  // Generic polynomial solvers
  poly_real_coeff_gsl p3;
  
  // quadratic solvers
  quadratic_real_coeff_gsl t1;
  quadratic_complex_std t2;
  
  // cubic solvers
  cubic_real_coeff_cern c1;
  cubic_real_coeff_gsl c2;
  cubic_complex_std c3;
  
  // quartic solvers
  quartic_real_coeff_cern q1;
  quartic_real_gsl q2;
  quartic_real_gsl2 q3;
  quartic_real_simple q4;
  quartic_complex_simple q5;

  double tt;

  std::vector<stats> slist(32);
  for(size_t i=0;i<32;i++) slist[i].w.resize(6);

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_quadratic_real_coeff(&t1,"quadratic_real_coeff_gsl quadratic",
			    1.0,slist[0]);
  slist[0].out();
  test_quadratic_real_coeff
    (&t2,"quadratic_complex_std quadratic",1.0,slist[1]);
  slist[1].out();
  test_quadratic_real_coeff(&p3,"poly_real_coeff_gsl quadratic",1.0,slist[2]);
  slist[2].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_quadratic_real_coeff
    (&t1,"quadratic_real_coeff_gsl quadratic (small odd coeffs.)",
     1.0e-5,slist[3]);
  slist[3].out();
  test_quadratic_real_coeff
    (&t2,"quadratic_complex_std quadratic (small odd coeffs.)",
     1.0e-5,slist[4]);
  slist[4].out();
  test_quadratic_real_coeff
    (&p3,"poly_real_coeff_gsl quadratic (small odd coeffs.)",1.0e-5,slist[5]);
  slist[5].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_quadratic_complex
    (&t2,"quadratic_complex_std quadratic (complex coeffs.)",slist[6]);
  slist[6].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_cubic_real_coeff(&c1,"real_coeff_cern cubic",1.0,slist[7]);
  slist[7].out();
  test_cubic_real_coeff(&c2,"cubic_real_coeff_gsl cubic",1.0,slist[8]);
  slist[8].out();
  test_cubic_real_coeff(&c3,"cubic_complex_std cubic",1.0,slist[9]);
  slist[9].out();
  test_cubic_real_coeff(&p3,"poly_real_coeff_gsl cubic",1.0,slist[10]);
  slist[10].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_cubic_real_coeff
    (&c1,"real_coeff_cern cubic (small odd coeffs.)",1.0e-5,slist[11]);
  slist[11].out();
  test_cubic_real_coeff
    (&c2,"cubic_real_coeff_gsl cubic (small odd coeffs.)",1.0e-5,slist[12]);
  slist[12].out();
  test_cubic_real_coeff
    (&c3,"cubic_complex_std cubic (small odd coeffs.)",1.0e-5,slist[13]);
  slist[13].out();
  test_cubic_real_coeff
    (&p3,"poly_real_coeff_gsl cubic (small odd coeffs.)",1.0e-5,slist[14]);
  slist[14].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_cubic_complex
    (&c3,"cubic_complex_std cubic (complex coeffs.)",slist[15]);
  slist[15].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_quartic_real
    (&q1,"cern_real_coeff quartic (real roots)",1.0,slist[16]);
  slist[16].out();
  test_quartic_real
    (&q2,"quartic_real_gsl quartic (real roots)",1.0,slist[17]);
  slist[17].out();
  test_quartic_real
    (&q3,"quartic_real_gsl2 quartic (real roots)",1.0,slist[18]);
  slist[18].out();
  test_quartic_real
    (&q4,"quartic_real_simple quartic (real roots)",1.0,slist[19]);
  slist[19].out();
  test_quartic_real
    (&q5,"quartic_complex_simple quartic (real roots)",1.0,slist[20]);
  slist[20].out();
  test_quartic_real
    (&p3,"poly_real_coeff_gsl quartic (real roots)",1.0,slist[21]);
  slist[21].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_quartic_real
    (&q1,"real_coeff_cern quartic (real roots, small odd coeffs.)",
     1.0e-5,slist[22]);
  slist[22].out();
  test_quartic_real
    (&q2,"quartic_real_gsl quartic (real roots, small odd coeffs.)",
     1.0e-5,slist[23]);
  slist[23].out();
  test_quartic_real
    (&q3,"quartic_real_gsl2 quartic (real roots, small odd coeffs.)",
     1.0e-5,slist[24]);
  slist[24].out();
  test_quartic_real
    (&q4,"quartic_real_simple quartic (real roots, small odd coeffs.)",
     1.0e-5,slist[25]);
  slist[25].out();
  test_quartic_real
    (&q5,"quartic_complex_simple quartic (real roots, small odd coeffs.)",
     1.0e-5,slist[26]);
  slist[26].out();
  test_quartic_real
    (&p3,"poly_real_coeff_gsl quartic (real roots, small odd coeffs.)",
     1.0e-5,slist[27]);
  slist[27].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_quartic_real_coeff
    (&q1,"real_coeff_cern quartic (complex roots)",slist[28]);
  slist[28].out();
  test_quartic_real_coeff
    (&q5,"quartic_complex_simple quartic (complex roots)",slist[29]);
  slist[29].out();
  test_quartic_real_coeff
    (&p3,"poly_real_coeff_gsl quartic (complex roots)",slist[30]);
  slist[30].out();

  cout << "----------------------------------------"
       << "---------------------------------------" << endl;
  test_quartic_complex
    (&q5,"quartic_complex_simple quartic (complex coeffs.)",slist[31]);
  slist[31].out();

  return 0;
}


