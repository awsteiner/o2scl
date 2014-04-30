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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/cfl6_eos.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

cfl6_eos::cfl6_eos() {
  KD=0.0;
  
  iprop6=gsl_matrix_complex_alloc(mat_size,mat_size);
  eivec6=gsl_matrix_complex_alloc(mat_size,mat_size);
  dipdgapu.resize(mat_size,mat_size);
  dipdgapd.resize(mat_size,mat_size);
  dipdgaps.resize(mat_size,mat_size);
  dipdqqu.resize(mat_size,mat_size);
  dipdqqd.resize(mat_size,mat_size);
  dipdqqs.resize(mat_size,mat_size);
  eval6=gsl_vector_alloc(mat_size);
  w6=gsl_eigen_hermv_alloc(mat_size);
  
  kdlimit=1.0e-6;
}

cfl6_eos::~cfl6_eos() {
  gsl_matrix_complex_free(iprop6);
  gsl_matrix_complex_free(eivec6);
  gsl_vector_free(eval6);
  gsl_eigen_hermv_free(w6);

}

int cfl6_eos::set_masses() {
  if (fixed_mass) {
    up->ms=up->m+KD/4.0/GD/GD*up->del*up->del;
    down->ms=down->m+KD/4.0/GD/GD*down->del*down->del;
    strange->ms=strange->m+KD/4.0/GD/GD*strange->del*strange->del;
  } else {
    up->ms=up->m-4.0*G*up->qq+2.0*K*down->qq*strange->qq+
      KD/4.0/GD/GD*up->del*up->del;
    down->ms=down->m-4.0*G*down->qq+2.0*K*up->qq*strange->qq+
      KD/4.0/GD/GD*down->del*down->del;
    strange->ms=strange->m-4.0*G*strange->qq+2.0*K*down->qq*up->qq+
      KD/4.0/GD/GD*strange->del*strange->del;
  }
  if (up->ms<up->m) up->ms=up->m;
  if (down->ms<down->m) down->ms=down->m;
  if (strange->ms<strange->m) strange->ms=strange->m;
  return 0;
}

int cfl6_eos::calc_eq_temp_p(quark &u, quark &d, quark &s,
			      double &qq1, double &qq2, double &qq3, 
			      double &gap1, double &gap2, double &gap3, 
			      double mu3, double mu8, 
			      double &n3, double &n8, 
			      thermo &qb, const double ltemper) {
  
  if (fabs(KD)<kdlimit) {
    return cfl_njl_eos::calc_eq_temp_p(u,d,s,qq1,qq2,qq3,gap1,gap2,gap3,
					mu3,mu8,n3,n8,qb,ltemper);
  }
  
  up=&u;
  down=&d;
  strange=&s;
  
  if (fromqq==false) {
    O2SCL_ERR2_RET("Does not work with fromqq=false ",
		   "in cfl6_eos::calc_eq_temp_p().",exc_efailed);
  }

  set_masses();
  
  // static variables for communication with tpot
  temper=ltemper;
  smu3=mu3;
  smu8=mu8;

  // Compute the integrals

  ubvector res(13);
  double err2;
  integ_err(0.0,L,13,res,err2);
  qb.pr=-res[0];
  qb.en=res[1];
  u.n=-res[2];
  d.n=-res[3];
  s.n=-res[4];
  qq1=res[5];
  qq2=res[6];
  qq3=res[7];
  gap1=res[8];
  gap2=res[9];
  gap3=res[10];
  n3=res[11];
  n8=res[12];

  // Add the remaining terms and their corresponding contributions
  // to the derivatives

  qb.pr-=B0;
  
  if (!fixed_mass) {
    qb.pr-=u.qq*u.qq*2.0*G;
    qb.pr-=d.qq*d.qq*2.0*G;
    qb.pr-=s.qq*s.qq*2.0*G;
    qb.pr+=u.qq*d.qq*s.qq*4.0*K;
    qb.pr-=KD/GD/GD/2.0*(u.del*u.del*u.qq+d.del*d.del*d.qq+
			 s.del*s.del*s.qq);
  }
  qb.pr-=u.del*u.del/4.0/GD;
  qb.pr-=d.del*d.del/4.0/GD;
  qb.pr-=s.del*s.del/4.0/GD;

  qq1+=4.0*G*u.qq-4.0*K*d.qq*s.qq+KD/GD/GD/2.0*u.del*u.del;
  qq2+=4.0*G*d.qq-4.0*K*u.qq*s.qq+KD/GD/GD/2.0*d.del*d.del;
  qq3+=4.0*G*s.qq-4.0*K*u.qq*d.qq+KD/GD/GD/2.0*s.del*s.del;
  
  gap1+=u.del/2.0/GD+KD/GD/GD*u.del*u.qq;
  gap2+=d.del/2.0/GD+KD/GD/GD*d.del*d.qq;
  gap3+=s.del/2.0/GD+KD/GD/GD*s.del*s.qq;

  // Compute the energy densities 
  
  u.ed=-u.pr+u.n*u.mu;
  d.ed=-d.pr+d.n*d.mu;
  s.ed=-s.pr+s.n*s.mu;
  qb.ed=qb.en*temper-qb.pr+u.mu*u.n+d.mu*d.n+s.mu*s.n;

  return 0;
}

int cfl6_eos::integrands(double p, double res[]) {
  
  if (fabs(KD)<kdlimit) {
    return cfl_njl_eos::integrands(p,res);
  }

  if (integ_test) {
    res[0]=sin(p);
    res[1]=sin(p);
    res[2]=cos(p);
    return 0;
  }

  int k;
  double egv[36];
  double dedmuu[36], dedmud[36], dedmus[36];
  double dedqqu[36], dedqqd[36], dedqqs[36], deds[36];
  double dedd[36], dedu[36], dedmu3[36], dedmu8[36];
  
  eigenvalues6(p,smu3,smu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,dedqqs,
	       dedu,dedd,deds,dedmu3,dedmu8);
  
  res[0]=0.0;
  res[1]=0.0;
  res[2]=0.0;
  res[3]=0.0;
  res[4]=0.0;
  res[5]=0.0;
  res[6]=0.0;
  res[7]=0.0;
  res[8]=0.0;
  res[9]=0.0;
  res[10]=0.0;
  res[11]=0.0;
  res[12]=0.0;
  for (k = 0; k < 36; k++) {
    if (zerot) {
      res[0]+=fabs(egv[k]/2.0);
      if (egv[k]>0.0) {
	res[2]+=dedmuu[k]/2.0;
	res[3]+=dedmud[k]/2.0;
	res[4]+=dedmus[k]/2.0;
	res[5]+=dedqqu[k]/2.0;
	res[6]+=dedqqd[k]/2.0;
	res[7]+=dedqqs[k]/2.0;
	res[8]+=dedu[k]/2.0;
	res[9]+=dedd[k]/2.0;
	res[10]+=deds[k]/2.0;
	res[11]+=dedmu3[k]/2.0;
	res[12]+=dedmu8[k]/2.0;
      } else {
	res[2]-=dedmuu[k]/2.0;
	res[3]-=dedmud[k]/2.0;
	res[4]-=dedmus[k]/2.0;
	res[5]-=dedqqu[k]/2.0;
	res[6]-=dedqqd[k]/2.0;
	res[7]-=dedqqs[k]/2.0;
	res[8]-=dedu[k]/2.0;
	res[9]-=dedd[k]/2.0;
	res[10]-=deds[k]/2.0;
	res[11]-=dedmu3[k]/2.0;
	res[12]-=dedmu8[k]/2.0;
      }
    } else {
      res[0]+=egv[k]/2.0;
      res[1]+=egv[k];
      res[2]+=dedmuu[k]/2.0;
      res[3]+=dedmud[k]/2.0;
      res[4]+=dedmus[k]/2.0;
      res[5]+=dedqqu[k]/2.0;
      res[6]+=dedqqd[k]/2.0;
      res[7]+=dedqqs[k]/2.0;
      res[8]+=dedu[k]/2.0;
      res[9]+=dedd[k]/2.0;
      res[10]+=deds[k]/2.0;
      res[11]+=dedmu3[k]/2.0;
      res[12]+=dedmu8[k]/2.0;
      if ((temper==0.0 && egv[k]<0.0) || egv[k]/temper<-30.0) {
	res[0]-=egv[k];
	res[2]-=dedmuu[k];
	res[3]-=dedmud[k];
	res[4]-=dedmus[k];
	res[5]-=dedqqu[k];
	res[6]-=dedqqd[k];
	res[7]-=dedqqs[k];
	res[8]-=dedu[k];
	res[9]-=dedd[k];
	res[10]-=deds[k];
	res[11]-=dedmu3[k];
	res[12]-=dedmu8[k];
      } else if (temper!=0.0 && egv[k]/temper<30.0) {
	res[0]+=temper*log(1.0+exp(-egv[k]/temper));
	res[2]-=1.0/(1.0+exp(egv[k]/temper))*dedmuu[k];
	res[3]-=1.0/(1.0+exp(egv[k]/temper))*dedmud[k];
	res[4]-=1.0/(1.0+exp(egv[k]/temper))*dedmus[k];
	res[5]-=1.0/(1.0+exp(egv[k]/temper))*dedqqu[k];
	res[6]-=1.0/(1.0+exp(egv[k]/temper))*dedqqd[k];
	res[7]-=1.0/(1.0+exp(egv[k]/temper))*dedqqs[k];
	res[8]-=1.0/(1.0+exp(egv[k]/temper))*dedu[k];
	res[9]-=1.0/(1.0+exp(egv[k]/temper))*dedd[k];
	res[10]-=1.0/(1.0+exp(egv[k]/temper))*deds[k];
	res[11]-=1.0/(1.0+exp(egv[k]/temper))*dedmu3[k];
	res[12]-=1.0/(1.0+exp(egv[k]/temper))*dedmu8[k];
      }
      if (temper!=0.0 && egv[k]/temper<30.0) {
	res[1]+=log(1.0+exp(-egv[k]/temper))+
	  egv[k]/temper/(1.0+exp(egv[k]/temper));
      }
    }
  }
  
  // Multiply by the integration weights and the 
  // Jacobian (momentum^2)
  res[0]*=-p*p*0.5/pi2;
  res[1]*=p*p*0.5/pi2;
  res[2]*=-p*p*0.5/pi2;
  res[3]*=-p*p*0.5/pi2;
  res[4]*=-p*p*0.5/pi2;
  res[5]*=-p*p*0.5/pi2;
  res[6]*=-p*p*0.5/pi2;
  res[7]*=-p*p*0.5/pi2;
  res[8]*=-p*p*0.5/pi2;
  res[9]*=-p*p*0.5/pi2;
  res[10]*=-p*p*0.5/pi2;
  res[11]*=-p*p*0.5/pi2;
  res[12]*=-p*p*0.5/pi2;

  return 0;
}

int cfl6_eos::test_derivatives(double lmom, double mu3, double mu8,
				test_mgr &t) {
  double egv[36];
  double dedmuu[36], dedmud[36], dedmus[36];
  double dedqqu[36], dedqqd[36], dedqqs[36], deds[36];
  double dedd[36], dedu[36], dedmu3[36], dedmu8[36];
  double h=1.0e-7;
  double d1[36], d2[36], tmpt[36];
  
  set_masses();

  // We have to do a lot of sorting here since the ordering of the
  // eigenvalues is not necessarily the same for subsequent calls to
  // the function eigenvalues6(). The same sorting should probably
  // be added to cfl_njl_eos::test_derivatives().

  // Check muu
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmuu,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    up->mu+=h;
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],3.0e-4,"muu");
      }
    }
    up->mu-=h;
  }
  // Check mud
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmud,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    down->mu+=h;
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],3.0e-4,"mud");
      }
    }
    down->mu-=h;
  }
  // Check mus
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmus,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    strange->mu+=h;
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],1.0e-3,"mus");
      }
    }
    strange->mu-=h;
  }
  // Check qqu
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedqqu,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    up->qq+=h;
    set_masses();
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],1.0e-4,"qqu");
      }
    }
    up->qq-=h;
    set_masses();
  }
  // Check qqd
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedqqd,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    down->qq+=h;
    set_masses();
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],1.0e-4,"qqd");
      }
    }
    down->qq-=h;
    set_masses();
  }
  // Check qqs
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedqqs,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    strange->qq+=h;
    set_masses();
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],1.5e-4,"qqs");
      }
    }
    strange->qq-=h;
    set_masses();
  }

  // Check del matrices
  ubmatrix_complex gu1(36,36), gu2(36,36), dgg(36,36);

  make_matrices(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      gsl_complex tmp=gsl_matrix_complex_get(iprop6,i,j);
      std::complex<double> tmp2(GSL_REAL(tmp),GSL_IMAG(tmp));
      gu1(i,j)=tmp2;
      dgg(i,j)=dipdgapu(i,j);
    }
  }
  up->del+=h;
  set_masses();
  make_matrices(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      gsl_complex tmp=gsl_matrix_complex_get(iprop6,i,j);
      std::complex<double> tmp2(GSL_REAL(tmp),GSL_IMAG(tmp));
      gu2(i,j)=tmp2;
    }
  }
  up->del-=h;
  set_masses();
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      if (i>j) {
	t.test_rel(dgg(i,j).real(),(gu2(i,j).real()-gu1(i,j).real())/h,
		   1.0e-5,"gapud");
	t.test_rel(dgg(i,j).imag(),(gu2(i,j).imag()-gu1(i,j).imag())/h,
		   1.0e-5,"gapud");
      }
    }
  }

  make_matrices(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      gsl_complex tmp=gsl_matrix_complex_get(iprop6,i,j);
      std::complex<double> tmp2(GSL_REAL(tmp),GSL_IMAG(tmp));
      gu1(i,j)=tmp2;
      dgg(i,j)=dipdgapd(i,j);
    }
  }
  down->del+=h;
  set_masses();
  make_matrices(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      gsl_complex tmp=gsl_matrix_complex_get(iprop6,i,j);
      std::complex<double> tmp2(GSL_REAL(tmp),GSL_IMAG(tmp));
      gu2(i,j)=tmp2;
    }
  }
  down->del-=h;
  set_masses();
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      if (i>j) {
	t.test_rel(dgg(i,j).real(),(gu2(i,j).real()-gu1(i,j).real())/h,
		   1.0e-5,"gapud");
	t.test_rel(dgg(i,j).imag(),(gu2(i,j).imag()-gu1(i,j).imag())/h,
		   1.0e-5,"gapud");
      }
    }
  }

  make_matrices(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      gsl_complex tmp=gsl_matrix_complex_get(iprop6,i,j);
      std::complex<double> tmp2(GSL_REAL(tmp),GSL_IMAG(tmp));
      gu1(i,j)=tmp2;
      dgg(i,j)=dipdgaps(i,j);
    }
  }
  strange->del+=h;
  set_masses();
  make_matrices(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      gsl_complex tmp=gsl_matrix_complex_get(iprop6,i,j);
      std::complex<double> tmp2(GSL_REAL(tmp),GSL_IMAG(tmp));
      gu2(i,j)=tmp2;
    }
  }
  strange->del-=h;
  set_masses();
  for(size_t i=0;i<36;i++) {
    for(size_t j=0;j<36;j++) {
      if (i>j) {
	t.test_rel(dgg(i,j).real(),(gu2(i,j).real()-gu1(i,j).real())/h,
		   1.0e-5,"gapud");
	t.test_rel(dgg(i,j).imag(),(gu2(i,j).imag()-gu1(i,j).imag())/h,
		   1.0e-5,"gapud");
      }
    }
  }

  // Check delu
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedu,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    up->del+=h;
    set_masses();
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],1.0e-4,"delu");
      }
    }
    up->del-=h;
    set_masses();
  }
  // Check deld
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedd,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    down->del+=h;
    set_masses();
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],1.0e-4,"deld");
      }
    }
    down->del-=h;
    set_masses();
  }
  // Check dels
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,deds,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    strange->del+=h;
    set_masses();
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],1.0e-4,"dels");
      }
    }
    strange->del-=h;
    set_masses();
  }
  // Check mu3
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmu3,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    mu3+=h;
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],5.0e-3,"mu3");
      }
    }
    mu3-=h;
  }
  // Check mu8
  {
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmu8,d2);
    vector_sort<double[36],double>(36,d1);
    vector_sort<double[36],double>(36,d2);
    mu8+=h;
    eigenvalues6(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedqqu,dedqqd,
		 dedqqs,dedu,dedd,deds,dedmu3,dedmu8);
    vector_sort<double[36],double>(36,egv);
    for(size_t i=0;i<36;i++) {
      tmpt[i]=(egv[i]-d1[i])/h;
    }
    vector_sort<double[36],double>(36,tmpt);
    for(size_t i=0;i<36;i++) {
      if (fabs(d2[i])>1.0e-4) {
	t.test_rel(d2[i],tmpt[i],4.0e-4,"mu8");
      }
    }
    mu8-=h;
  }

  double res[13];
  double x1,x2,x3,x4,x5,x6,n3,n8;
  // Check derivative wrt muu in integrands() and calc_eq_temp_p()
  {
    integrands(lmom,res);
    double om1=res[0];
    double dmuu=res[2];
    up->mu+=h;
    integrands(lmom,res);
    double om2=res[0];
    up->mu-=h;
    t.test_rel((om2-om1)/h,dmuu,1.0e-6,"dmuu");
    
    calc_eq_temp_p(*up,*down,*strange,x1,x2,x3,x4,x5,x6,mu3,mu8,
		   n3,n8,*eos_thermo,4.0/hc_mev_fm);
    double pr1=eos_thermo->pr;
    double nu=up->n;
    up->mu+=h;
    calc_eq_temp_p(*up,*down,*strange,x1,x2,x3,x4,x5,x6,mu3,mu8,
		   n3,n8,*eos_thermo,4.0/hc_mev_fm);
    double pr2=eos_thermo->pr;
    up->mu-=h;
    t.test_rel((pr2-pr1)/h,nu,1.0e-5,"nup");
  }
    
  return 0;
}

int cfl6_eos::eigenvalues6(double lmom, double mu3, 
			    double mu8, double egv[36],
			    double dedmuu[36], double dedmud[36],
			    double dedmus[36], double dedqqu[36], 
			    double dedqqd[36], double dedqqs[36],
			    double dedu[36], double dedd[36],
			    double deds[36], double dedmu3[36],
			    double dedmu8[36]) {
  
  int k;
  const double mu=up->ms, md=down->ms, ms=strange->ms;
  const double muu=up->mu, mud=down->mu, mus=strange->mu;
  const double du=up->del, dd=down->del, ds=strange->del;
  const double qqu=up->qq, qqd=down->qq, qqs=strange->qq;

  double mur, mug, mub;
  mur=mu3+mu8/sqrt(3.0);
  mug=-mu3+mu8/sqrt(3.0);
  mub=-2.0*mu8/sqrt(3.0);
  
  for(k=0;k<36;k++) {
    egv[k]=0.0;
    dedmuu[k]=0.0;
    dedmud[k]=0.0;
    dedmus[k]=0.0;
    dedqqu[k]=0.0;
    dedqqd[k]=0.0;
    dedqqs[k]=0.0;
    dedu[k]=0.0;
    dedd[k]=0.0;
    deds[k]=0.0;
    dedmu3[k]=0.0;
    dedmu8[k]=0.0;
  }
  
  ubvector mmuu(mat_size), mmud(mat_size), mmus(mat_size);
  ubvector mmu(mat_size), mmd(mat_size), mms(mat_size);
  ubvector mmur(mat_size), mmug(mat_size), mmub(mat_size);
  for(k=0;k<36;k++) {
    mmuu[k]=0.0;
    mmud[k]=0.0;
    mmus[k]=0.0;
    mmu[k]=0.0;
    mmd[k]=0.0;
    mms[k]=0.0;
    mmur[k]=0.0;
    mmug[k]=0.0;
    mmub[k]=0.0;
  }

  gsl_complex zero={{0.0,0.0}};
  std::complex<double> zero2(0.0,0.0);
  for(int i=0;i<mat_size;i++) {
    for(int j=0;j<mat_size;j++) {
      gsl_matrix_complex_set(iprop6,i,j,zero);
      dipdgapu(i,j)=zero2;
      dipdgapd(i,j)=zero2;
      dipdgaps(i,j)=zero2;
      dipdqqu(i,j)=zero2;
      dipdqqd(i,j)=zero2;
      dipdqqs(i,j)=zero2;
    }
  }
  
  /// Set diagonal elements

  gsl_complex tmp={{0.0,0.0}};

  tmp=gsl_matrix_complex_get(iprop6,0,0);
  GSL_REAL(tmp)=-mu-muu-mur;
  gsl_matrix_complex_set(iprop6,0,0,tmp);

  tmp=gsl_matrix_complex_get(iprop6,1,1);
  GSL_REAL(tmp)=mu-muu-mur;
  gsl_matrix_complex_set(iprop6,1,1,tmp);

  tmp=gsl_matrix_complex_get(iprop6,2,2);
  GSL_REAL(tmp)=-mu+muu+mur;
  gsl_matrix_complex_set(iprop6,2,2,tmp);

  tmp=gsl_matrix_complex_get(iprop6,3,3);
  GSL_REAL(tmp)=mu+muu+mur;
  gsl_matrix_complex_set(iprop6,3,3,tmp);

  mmuu[0]=-1.0;
  mmuu[1]=-1.0;
  mmuu[2]=1.0;
  mmuu[3]=1.0;
  mmur[0]=-1.0;
  mmur[1]=-1.0;
  mmur[2]=1.0;
  mmur[3]=1.0;
  mmu[0]=-1.0;
  mmu[1]=1.0;
  mmu[2]=-1.0;
  mmu[3]=1.0;

  tmp=gsl_matrix_complex_get(iprop6,4,4);
  GSL_REAL(tmp)=-mu-muu-mug;
  gsl_matrix_complex_set(iprop6,4,4,tmp);

  tmp=gsl_matrix_complex_get(iprop6,5,5);
  GSL_REAL(tmp)=mu-muu-mug;
  gsl_matrix_complex_set(iprop6,5,5,tmp);

  tmp=gsl_matrix_complex_get(iprop6,6,6);
  GSL_REAL(tmp)=-mu+muu+mug;
  gsl_matrix_complex_set(iprop6,6,6,tmp);

  tmp=gsl_matrix_complex_get(iprop6,7,7);
  GSL_REAL(tmp)=mu+muu+mug;
  gsl_matrix_complex_set(iprop6,7,7,tmp);

  mmuu[4]=-1.0;
  mmuu[5]=-1.0;
  mmuu[6]=1.0;
  mmuu[7]=1.0;
  mmug[4]=-1.0;
  mmug[5]=-1.0;
  mmug[6]=1.0;
  mmug[7]=1.0;
  mmu[4]=-1.0;
  mmu[5]=1.0;
  mmu[6]=-1.0;
  mmu[7]=1.0;

  
  tmp=gsl_matrix_complex_get(iprop6,8,8);
  GSL_REAL(tmp)=-mu-muu-mub;
  gsl_matrix_complex_set(iprop6,8,8,tmp);

  tmp=gsl_matrix_complex_get(iprop6,9,9);
  GSL_REAL(tmp)=mu-muu-mub;
  gsl_matrix_complex_set(iprop6,9,9,tmp);

  tmp=gsl_matrix_complex_get(iprop6,10,10);
  GSL_REAL(tmp)=-mu+muu+mub;
  gsl_matrix_complex_set(iprop6,10,10,tmp);

  tmp=gsl_matrix_complex_get(iprop6,11,11);
  GSL_REAL(tmp)=mu+muu+mub;
  gsl_matrix_complex_set(iprop6,11,11,tmp);
  
  mmuu[8]=-1.0;
  mmuu[9]=-1.0;
  mmuu[10]=1.0;
  mmuu[11]=1.0;
  mmub[8]=-1.0;
  mmub[9]=-1.0;
  mmub[10]=1.0;
  mmub[11]=1.0;
  mmu[8]=-1.0;
  mmu[9]=1.0;
  mmu[10]=-1.0;
  mmu[11]=1.0;

  tmp=gsl_matrix_complex_get(iprop6,12,12);
  GSL_REAL(tmp)=-md-mud-mur;
  gsl_matrix_complex_set(iprop6,12,12,tmp);

  tmp=gsl_matrix_complex_get(iprop6,13,13);
  GSL_REAL(tmp)=md-mud-mur;
  gsl_matrix_complex_set(iprop6,13,13,tmp);

  tmp=gsl_matrix_complex_get(iprop6,14,14);
  GSL_REAL(tmp)=-md+mud+mur;
  gsl_matrix_complex_set(iprop6,14,14,tmp);

  tmp=gsl_matrix_complex_get(iprop6,15,15);
  GSL_REAL(tmp)=md+mud+mur;
  gsl_matrix_complex_set(iprop6,15,15,tmp);

  mmud[12]=-1.0;
  mmud[13]=-1.0;
  mmud[14]=1.0;
  mmud[15]=1.0;
  mmur[12]=-1.0;
  mmur[13]=-1.0;
  mmur[14]=1.0;
  mmur[15]=1.0;
  mmd[12]=-1.0;
  mmd[13]=1.0;
  mmd[14]=-1.0;
  mmd[15]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,16,16);
  GSL_REAL(tmp)=-md-mud-mug;
  gsl_matrix_complex_set(iprop6,16,16,tmp);

  tmp=gsl_matrix_complex_get(iprop6,17,17);
  GSL_REAL(tmp)=md-mud-mug;
  gsl_matrix_complex_set(iprop6,17,17,tmp);

  tmp=gsl_matrix_complex_get(iprop6,18,18);
  GSL_REAL(tmp)=-md+mud+mug;
  gsl_matrix_complex_set(iprop6,18,18,tmp);

  tmp=gsl_matrix_complex_get(iprop6,19,19);
  GSL_REAL(tmp)=md+mud+mug;
  gsl_matrix_complex_set(iprop6,19,19,tmp);
  
  mmud[16]=-1.0;
  mmud[17]=-1.0;
  mmud[18]=1.0;
  mmud[19]=1.0;
  mmug[16]=-1.0;
  mmug[17]=-1.0;
  mmug[18]=1.0;
  mmug[19]=1.0;
  mmd[16]=-1.0;
  mmd[17]=1.0;
  mmd[18]=-1.0;
  mmd[19]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,20,20);
  GSL_REAL(tmp)=-md-mud-mub;
  gsl_matrix_complex_set(iprop6,20,20,tmp);

  tmp=gsl_matrix_complex_get(iprop6,21,21);
  GSL_REAL(tmp)=md-mud-mub;
  gsl_matrix_complex_set(iprop6,21,21,tmp);

  tmp=gsl_matrix_complex_get(iprop6,22,22);
  GSL_REAL(tmp)=-md+mud+mub;
  gsl_matrix_complex_set(iprop6,22,22,tmp);

  tmp=gsl_matrix_complex_get(iprop6,23,23);
  GSL_REAL(tmp)=md+mud+mub;
  gsl_matrix_complex_set(iprop6,23,23,tmp);

  mmud[20]=-1.0;
  mmud[21]=-1.0;
  mmud[22]=1.0;
  mmud[23]=1.0;
  mmub[20]=-1.0;
  mmub[21]=-1.0;
  mmub[22]=1.0;
  mmub[23]=1.0;
  mmd[20]=-1.0;
  mmd[21]=1.0;
  mmd[22]=-1.0;
  mmd[23]=1.0;

  
  tmp=gsl_matrix_complex_get(iprop6,24,24);
  GSL_REAL(tmp)=-ms-mus-mur;
  gsl_matrix_complex_set(iprop6,24,24,tmp);

  tmp=gsl_matrix_complex_get(iprop6,25,25);
  GSL_REAL(tmp)=ms-mus-mur;
  gsl_matrix_complex_set(iprop6,25,25,tmp);

  tmp=gsl_matrix_complex_get(iprop6,26,26);
  GSL_REAL(tmp)=-ms+mus+mur;
  gsl_matrix_complex_set(iprop6,26,26,tmp);

  tmp=gsl_matrix_complex_get(iprop6,27,27);
  GSL_REAL(tmp)=ms+mus+mur;
  gsl_matrix_complex_set(iprop6,27,27,tmp);

  mmus[24]=-1.0;
  mmus[25]=-1.0;
  mmus[26]=1.0;
  mmus[27]=1.0;
  mmur[24]=-1.0;
  mmur[25]=-1.0;
  mmur[26]=1.0;
  mmur[27]=1.0;
  mms[24]=-1.0;
  mms[25]=1.0;
  mms[26]=-1.0;
  mms[27]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,28,28);
  GSL_REAL(tmp)=-ms-mus-mug;
  gsl_matrix_complex_set(iprop6,28,28,tmp);

  tmp=gsl_matrix_complex_get(iprop6,29,29);
  GSL_REAL(tmp)=ms-mus-mug;
  gsl_matrix_complex_set(iprop6,29,29,tmp);

  tmp=gsl_matrix_complex_get(iprop6,30,30);
  GSL_REAL(tmp)=-ms+mus+mug;
  gsl_matrix_complex_set(iprop6,30,30,tmp);

  tmp=gsl_matrix_complex_get(iprop6,31,31);
  GSL_REAL(tmp)=ms+mus+mug;
  gsl_matrix_complex_set(iprop6,31,31,tmp);

  mmus[28]=-1.0;
  mmus[29]=-1.0;
  mmus[30]=1.0;
  mmus[31]=1.0;
  mmug[28]=-1.0;
  mmug[29]=-1.0;
  mmug[30]=1.0;
  mmug[31]=1.0;
  mms[28]=-1.0;
  mms[29]=1.0;
  mms[30]=-1.0;
  mms[31]=1.0;

  tmp=gsl_matrix_complex_get(iprop6,32,32);
  GSL_REAL(tmp)=-ms-mus-mub;
  gsl_matrix_complex_set(iprop6,32,32,tmp);

  tmp=gsl_matrix_complex_get(iprop6,33,33);
  GSL_REAL(tmp)=ms-mus-mub;
  gsl_matrix_complex_set(iprop6,33,33,tmp);

  tmp=gsl_matrix_complex_get(iprop6,34,34);
  GSL_REAL(tmp)=-ms+mus+mub;
  gsl_matrix_complex_set(iprop6,34,34,tmp);

  tmp=gsl_matrix_complex_get(iprop6,35,35);
  GSL_REAL(tmp)=ms+mus+mub;
  gsl_matrix_complex_set(iprop6,35,35,tmp);

  mmus[32]=-1.0;
  mmus[33]=-1.0;
  mmus[34]=1.0;
  mmus[35]=1.0;
  mmub[32]=-1.0;
  mmub[33]=-1.0;
  mmub[34]=1.0;
  mmub[35]=1.0;
  mms[32]=-1.0;
  mms[33]=1.0;
  mms[34]=-1.0;
  mms[35]=1.0;
  
  // Set momentum
  
  for(int kk=0;kk<35;kk+=2) {
    tmp=gsl_matrix_complex_get(iprop6,kk+1,kk);
    GSL_REAL(tmp)=-lmom;
    gsl_matrix_complex_set(iprop6,kk+1,kk,tmp);
  }
  
  // Set gaps

  if (fixed_mass) {

    // Set gap terms

    tmp=gsl_matrix_complex_get(iprop6,28,23);
    GSL_IMAG(tmp)=du;
    gsl_matrix_complex_set(iprop6,28,23,tmp);

    tmp=gsl_matrix_complex_get(iprop6,30,21);
    GSL_IMAG(tmp)=du;
    gsl_matrix_complex_set(iprop6,30,21,tmp);

    tmp=gsl_matrix_complex_get(iprop6,33,18);
    GSL_IMAG(tmp)=du;
    gsl_matrix_complex_set(iprop6,33,18,tmp);

    tmp=gsl_matrix_complex_get(iprop6,35,16);
    GSL_IMAG(tmp)=du;
    gsl_matrix_complex_set(iprop6,35,16,tmp);

    // --

    tmp=gsl_matrix_complex_get(iprop6,29,22);
    GSL_IMAG(tmp)=-du;
    gsl_matrix_complex_set(iprop6,29,22,tmp);

    tmp=gsl_matrix_complex_get(iprop6,31,20);
    GSL_IMAG(tmp)=-du;
    gsl_matrix_complex_set(iprop6,31,20,tmp);

    tmp=gsl_matrix_complex_get(iprop6,32,19);
    GSL_IMAG(tmp)=-du;
    gsl_matrix_complex_set(iprop6,32,19,tmp);

    tmp=gsl_matrix_complex_get(iprop6,34,17);
    GSL_IMAG(tmp)=-du;
    gsl_matrix_complex_set(iprop6,34,17,tmp);

    // --

    tmp=gsl_matrix_complex_get(iprop6,24,11);
    GSL_IMAG(tmp)=dd;
    gsl_matrix_complex_set(iprop6,24,11,tmp);

    tmp=gsl_matrix_complex_get(iprop6,26,9);
    GSL_IMAG(tmp)=dd;
    gsl_matrix_complex_set(iprop6,26,9,tmp);

    tmp=gsl_matrix_complex_get(iprop6,33,2);
    GSL_IMAG(tmp)=dd;
    gsl_matrix_complex_set(iprop6,33,2,tmp);

    tmp=gsl_matrix_complex_get(iprop6,35,0);
    GSL_IMAG(tmp)=dd;
    gsl_matrix_complex_set(iprop6,35,0,tmp);

    // --			   

    tmp=gsl_matrix_complex_get(iprop6,25,10);
    GSL_IMAG(tmp)=-dd;
    gsl_matrix_complex_set(iprop6,25,10,tmp);

    tmp=gsl_matrix_complex_get(iprop6,27,8);
    GSL_IMAG(tmp)=-dd;
    gsl_matrix_complex_set(iprop6,27,8,tmp);

    tmp=gsl_matrix_complex_get(iprop6,32,3);
    GSL_IMAG(tmp)=-dd;
    gsl_matrix_complex_set(iprop6,32,3,tmp);

    tmp=gsl_matrix_complex_get(iprop6,34,1);
    GSL_IMAG(tmp)=-dd;
    gsl_matrix_complex_set(iprop6,34,1,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,12,7);
    GSL_IMAG(tmp)=ds;
    gsl_matrix_complex_set(iprop6,12,7,tmp);

    tmp=gsl_matrix_complex_get(iprop6,14,5);
    GSL_IMAG(tmp)=ds;
    gsl_matrix_complex_set(iprop6,14,5,tmp);

    tmp=gsl_matrix_complex_get(iprop6,17,2);
    GSL_IMAG(tmp)=ds;
    gsl_matrix_complex_set(iprop6,17,2,tmp);

    tmp=gsl_matrix_complex_get(iprop6,19,0);
    GSL_IMAG(tmp)=ds;
    gsl_matrix_complex_set(iprop6,19,0,tmp);

    // --

    tmp=gsl_matrix_complex_get(iprop6,13,6);
    GSL_IMAG(tmp)=-ds;
    gsl_matrix_complex_set(iprop6,13,6,tmp);

    tmp=gsl_matrix_complex_get(iprop6,15,4);
    GSL_IMAG(tmp)=-ds;
    gsl_matrix_complex_set(iprop6,15,4,tmp);

    tmp=gsl_matrix_complex_get(iprop6,16,3);
    GSL_IMAG(tmp)=-ds;
    gsl_matrix_complex_set(iprop6,16,3,tmp);

    tmp=gsl_matrix_complex_get(iprop6,18,1);
    GSL_IMAG(tmp)=-ds;
    gsl_matrix_complex_set(iprop6,18,1,tmp);

    // Set lower-left derivatives of gap terms

    std::complex<double> littlei(0.0,1.0);
    std::complex<double> mlittlei(0.0,-1.0);

    dipdgapu(28,23)=littlei;
    dipdgapu(30,21)=littlei;
    dipdgapu(33,18)=littlei;
    dipdgapu(35,16)=littlei;

    dipdgapu(29,22)=mlittlei;
    dipdgapu(31,20)=mlittlei;
    dipdgapu(32,19)=mlittlei;
    dipdgapu(34,17)=mlittlei;
    
    dipdgapd(24,11)=littlei;
    dipdgapd(26,9)=littlei;
    dipdgapd(33,2)=littlei;
    dipdgapd(35,0)=littlei;

    dipdgapd(25,10)=mlittlei;
    dipdgapd(27,8)=mlittlei;
    dipdgapd(32,3)=mlittlei;
    dipdgapd(34,1)=mlittlei;
    
    dipdgaps(12,7)=littlei;
    dipdgaps(14,5)=littlei;
    dipdgaps(17,2)=littlei;
    dipdgaps(19,0)=littlei;

    dipdgaps(13,6)=mlittlei;
    dipdgaps(15,4)=mlittlei;
    dipdgaps(16,3)=mlittlei;
    dipdgaps(18,1)=mlittlei;

    // Set upper-right derivatives of gap terms

    dipdgapu(23,28)=mlittlei;
    dipdgapu(21,30)=mlittlei;
    dipdgapu(18,33)=mlittlei;
    dipdgapu(16,35)=mlittlei;
    
    dipdgapu(22,29)=littlei;
    dipdgapu(20,31)=littlei;
    dipdgapu(19,32)=littlei;
    dipdgapu(17,34)=littlei;
    
    dipdgapd(11,24)=mlittlei;
    dipdgapd(9,26)=mlittlei;
    dipdgapd(2,33)=mlittlei;
    dipdgapd(0,35)=mlittlei;
    
    dipdgapd(10,25)=littlei;
    dipdgapd(8,27)=littlei;
    dipdgapd(3,32)=littlei;
    dipdgapd(1,34)=littlei;
    
    dipdgaps(7,12)=mlittlei;
    dipdgaps(5,14)=mlittlei;
    dipdgaps(2,17)=mlittlei;
    dipdgaps(0,19)=mlittlei;
    
    dipdgaps(6,13)=littlei;
    dipdgaps(4,15)=littlei;
    dipdgaps(3,16)=littlei;
    dipdgaps(1,18)=littlei;
    
  } else {

    // set gap terms
    double betau=1.0+KD/GD*qqu/3.0;
    double betad=1.0+KD/GD*qqd/3.0;
    double betas=1.0+KD/GD*qqs/3.0;
    double gammau=du*KD/GD/3.0;
    double gammad=dd*KD/GD/3.0;
    double gammas=ds*KD/GD/3.0;

    tmp=gsl_matrix_complex_get(iprop6,28,23);
    GSL_IMAG(tmp)=du*betau;
    gsl_matrix_complex_set(iprop6,28,23,tmp);

    tmp=gsl_matrix_complex_get(iprop6,30,21);
    GSL_IMAG(tmp)=du*betau;
    gsl_matrix_complex_set(iprop6,30,21,tmp);

    tmp=gsl_matrix_complex_get(iprop6,33,18);
    GSL_IMAG(tmp)=du*betau;
    gsl_matrix_complex_set(iprop6,33,18,tmp);

    tmp=gsl_matrix_complex_get(iprop6,35,16);
    GSL_IMAG(tmp)=du*betau;
    gsl_matrix_complex_set(iprop6,35,16,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,29,22);
    GSL_IMAG(tmp)=-du*betau;
    gsl_matrix_complex_set(iprop6,29,22,tmp);

    tmp=gsl_matrix_complex_get(iprop6,31,20);
    GSL_IMAG(tmp)=-du*betau;
    gsl_matrix_complex_set(iprop6,31,20,tmp);

    tmp=gsl_matrix_complex_get(iprop6,32,19);
    GSL_IMAG(tmp)=-du*betau;
    gsl_matrix_complex_set(iprop6,32,19,tmp);

    tmp=gsl_matrix_complex_get(iprop6,34,17);
    GSL_IMAG(tmp)=-du*betau;
    gsl_matrix_complex_set(iprop6,34,17,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,24,11);
    GSL_IMAG(tmp)=dd*betad;
    gsl_matrix_complex_set(iprop6,24,11,tmp);

    tmp=gsl_matrix_complex_get(iprop6,26,9);
    GSL_IMAG(tmp)=dd*betad;
    gsl_matrix_complex_set(iprop6,26,9,tmp);

    tmp=gsl_matrix_complex_get(iprop6,33,2);
    GSL_IMAG(tmp)=dd*betad;
    gsl_matrix_complex_set(iprop6,33,2,tmp);

    tmp=gsl_matrix_complex_get(iprop6,35,0);
    GSL_IMAG(tmp)=dd*betad;
    gsl_matrix_complex_set(iprop6,35,0,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,25,10);
    GSL_IMAG(tmp)=-dd*betad;
    gsl_matrix_complex_set(iprop6,25,10,tmp);

    tmp=gsl_matrix_complex_get(iprop6,27,8);
    GSL_IMAG(tmp)=-dd*betad;
    gsl_matrix_complex_set(iprop6,27,8,tmp);

    tmp=gsl_matrix_complex_get(iprop6,32,3);
    GSL_IMAG(tmp)=-dd*betad;
    gsl_matrix_complex_set(iprop6,32,3,tmp);

    tmp=gsl_matrix_complex_get(iprop6,34,1);
    GSL_IMAG(tmp)=-dd*betad;
    gsl_matrix_complex_set(iprop6,34,1,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,12,7);
    GSL_IMAG(tmp)=ds*betas;
    gsl_matrix_complex_set(iprop6,12,7,tmp);

    tmp=gsl_matrix_complex_get(iprop6,14,5);
    GSL_IMAG(tmp)=ds*betas;
    gsl_matrix_complex_set(iprop6,14,5,tmp);

    tmp=gsl_matrix_complex_get(iprop6,17,2);
    GSL_IMAG(tmp)=ds*betas;
    gsl_matrix_complex_set(iprop6,17,2,tmp);

    tmp=gsl_matrix_complex_get(iprop6,19,0);
    GSL_IMAG(tmp)=ds*betas;
    gsl_matrix_complex_set(iprop6,19,0,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,13,6);
    GSL_IMAG(tmp)=-ds*betas;
    gsl_matrix_complex_set(iprop6,13,6,tmp);
    
    tmp=gsl_matrix_complex_get(iprop6,15,4);
    GSL_IMAG(tmp)=-ds*betas;
    gsl_matrix_complex_set(iprop6,15,4,tmp);

    tmp=gsl_matrix_complex_get(iprop6,16,3);
    GSL_IMAG(tmp)=-ds*betas;
    gsl_matrix_complex_set(iprop6,16,3,tmp);

    tmp=gsl_matrix_complex_get(iprop6,18,1);
    GSL_IMAG(tmp)=-ds*betas;
    gsl_matrix_complex_set(iprop6,18,1,tmp);
    
    // Set lower-left derivatives of gap terms

    std::complex<double> littlei(0.0,1.0);

    dipdgapu(28,23)=littlei*betau;
    dipdgapu(30,21)=littlei*betau;
    dipdgapu(33,18)=littlei*betau;
    dipdgapu(35,16)=littlei*betau;

    dipdgapu(29,22)=littlei*(-betau);
    dipdgapu(31,20)=littlei*(-betau);
    dipdgapu(32,19)=littlei*(-betau);
    dipdgapu(34,17)=littlei*(-betau);
    
    dipdgapd(24,11)=littlei*betad;
    dipdgapd(26,9)=littlei*betad;
    dipdgapd(33,2)=littlei*betad;
    dipdgapd(35,0)=littlei*betad;

    dipdgapd(25,10)=littlei*(-betad);
    dipdgapd(27,8)=littlei*(-betad);
    dipdgapd(32,3)=littlei*(-betad);
    dipdgapd(34,1)=littlei*(-betad);
    
    dipdgaps(12,7)=littlei*betas;
    dipdgaps(14,5)=littlei*betas;
    dipdgaps(17,2)=littlei*betas;
    dipdgaps(19,0)=littlei*betas;
    
    dipdgaps(13,6)=littlei*(-betas);
    dipdgaps(15,4)=littlei*(-betas);
    dipdgaps(16,3)=littlei*(-betas);
    dipdgaps(18,1)=littlei*(-betas);

    dipdqqu(28,23)=littlei*gammau;
    dipdqqu(30,21)=littlei*gammau;
    dipdqqu(33,18)=littlei*gammau;
    dipdqqu(35,16)=littlei*gammau;

    dipdqqu(29,22)=littlei*(-gammau);
    dipdqqu(31,20)=littlei*(-gammau);
    dipdqqu(32,19)=littlei*(-gammau);
    dipdqqu(34,17)=littlei*(-gammau);
    
    dipdqqd(24,11)=littlei*gammad;
    dipdqqd(26,9)=littlei*gammad;
    dipdqqd(33,2)=littlei*gammad;
    dipdqqd(35,0)=littlei*gammad;

    dipdqqd(25,10)=littlei*(-gammad);
    dipdqqd(27,8)=littlei*(-gammad);
    dipdqqd(32,3)=littlei*(-gammad);
    dipdqqd(34,1)=littlei*(-gammad);
    
    dipdqqs(12,7)=littlei*gammas;
    dipdqqs(14,5)=littlei*gammas;
    dipdqqs(17,2)=littlei*gammas;
    dipdqqs(19,0)=littlei*gammas;
    
    dipdqqs(13,6)=littlei*(-gammas);
    dipdqqs(15,4)=littlei*(-gammas);
    dipdqqs(16,3)=littlei*(-gammas);
    dipdqqs(18,1)=littlei*(-gammas);

    // Set upper-right derivatives of gap terms

    dipdgapu(23,28)=littlei*(-betau);
    dipdgapu(21,30)=littlei*(-betau);
    dipdgapu(18,33)=littlei*(-betau);
    dipdgapu(16,35)=littlei*(-betau);
    
    dipdgapu(22,29)=littlei*betau;
    dipdgapu(20,31)=littlei*betau;
    dipdgapu(19,32)=littlei*betau;
    dipdgapu(17,34)=littlei*betau;
    
    dipdgapd(11,24)=littlei*(-betad);
    dipdgapd(9,26)=littlei*(-betad);
    dipdgapd(2,33)=littlei*(-betad);
    dipdgapd(0,35)=littlei*(-betad);
    
    dipdgapd(10,25)=littlei*betad;
    dipdgapd(8,27)=littlei*betad;
    dipdgapd(3,32)=littlei*betad;
    dipdgapd(1,34)=littlei*betad;
    
    dipdgaps(7,12)=littlei*(-betas);
    dipdgaps(5,14)=littlei*(-betas);
    dipdgaps(2,17)=littlei*(-betas);
    dipdgaps(0,19)=littlei*(-betas);
    
    dipdgaps(6,13)=littlei*betas;
    dipdgaps(4,15)=littlei*betas;
    dipdgaps(3,16)=littlei*betas;
    dipdgaps(1,18)=littlei*betas;

    dipdqqu(23,28)=littlei*(-gammau);
    dipdqqu(21,30)=littlei*(-gammau);
    dipdqqu(18,33)=littlei*(-gammau);
    dipdqqu(16,35)=littlei*(-gammau);
    
    dipdqqu(22,29)=littlei*gammau;
    dipdqqu(20,31)=littlei*gammau;
    dipdqqu(19,32)=littlei*gammau;
    dipdqqu(17,34)=littlei*gammau;
    
    dipdqqd(11,24)=littlei*(-gammad);
    dipdqqd(9,26)=littlei*(-gammad);
    dipdqqd(2,33)=littlei*(-gammad);
    dipdqqd(0,35)=littlei*(-gammad);
    
    dipdqqd(10,25)=littlei*gammad;
    dipdqqd(8,27)=littlei*gammad;
    dipdqqd(3,32)=littlei*gammad;
    dipdqqd(1,34)=littlei*gammad;
    
    dipdqqs(7,12)=littlei*(-gammas);
    dipdqqs(5,14)=littlei*(-gammas);
    dipdqqs(2,17)=littlei*(-gammas);
    dipdqqs(0,19)=littlei*(-gammas);
    
    dipdqqs(6,13)=littlei*gammas;
    dipdqqs(4,15)=littlei*gammas;
    dipdqqs(3,16)=littlei*gammas;
    dipdqqs(1,18)=littlei*gammas;
  }
  
  double acons=KD/4.0/GD/GD;
  double ua=du*acons;
  double da=dd*acons;
  double sa=ds*acons;
  double uda=du*dd*acons;
  double usa=du*ds*acons;
  double dsa=dd*ds*acons;

  // Set off-diagonal mass-like terms
  
  GSL_REAL(tmp)=0.0;
  GSL_IMAG(tmp)=0.0;

  tmp=gsl_matrix_complex_get(iprop6,12,0);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,12,0,tmp);

  tmp=gsl_matrix_complex_get(iprop6,15,3);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,15,3,tmp);

  tmp=gsl_matrix_complex_get(iprop6,16,4);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,16,4,tmp);

  tmp=gsl_matrix_complex_get(iprop6,19,7);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,19,7,tmp);

  tmp=gsl_matrix_complex_get(iprop6,20,8);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,20,8,tmp);

  tmp=gsl_matrix_complex_get(iprop6,23,11);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,23,11,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,13,1);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,13,1,tmp);

  tmp=gsl_matrix_complex_get(iprop6,14,2);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,14,2,tmp);

  tmp=gsl_matrix_complex_get(iprop6,17,5);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,17,5,tmp);

  tmp=gsl_matrix_complex_get(iprop6,18,6);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,18,6,tmp);

  tmp=gsl_matrix_complex_get(iprop6,21,9);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,21,9,tmp);

  tmp=gsl_matrix_complex_get(iprop6,22,10);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,22,10,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,24,0);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,24,0,tmp);

  tmp=gsl_matrix_complex_get(iprop6,27,3);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,27,3,tmp);

  tmp=gsl_matrix_complex_get(iprop6,28,4);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,28,4,tmp);

  tmp=gsl_matrix_complex_get(iprop6,31,7);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,31,7,tmp);

  tmp=gsl_matrix_complex_get(iprop6,32,8);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,32,8,tmp);

  tmp=gsl_matrix_complex_get(iprop6,35,11);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,35,11,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,25,1);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,25,1,tmp);

  tmp=gsl_matrix_complex_get(iprop6,26,2);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,26,2,tmp);

  tmp=gsl_matrix_complex_get(iprop6,29,5);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,29,5,tmp);

  tmp=gsl_matrix_complex_get(iprop6,30,6);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,30,6,tmp);

  tmp=gsl_matrix_complex_get(iprop6,33,9);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,33,9,tmp);

  tmp=gsl_matrix_complex_get(iprop6,34,10);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,34,10,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,24,12);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,24,12,tmp);

  tmp=gsl_matrix_complex_get(iprop6,27,15);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,27,15,tmp);

  tmp=gsl_matrix_complex_get(iprop6,28,16);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,28,16,tmp);

  tmp=gsl_matrix_complex_get(iprop6,31,19);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,31,19,tmp);

  tmp=gsl_matrix_complex_get(iprop6,32,20);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,32,20,tmp);

  tmp=gsl_matrix_complex_get(iprop6,35,23);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,35,23,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,25,13);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,25,13,tmp);

  tmp=gsl_matrix_complex_get(iprop6,26,14);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,26,14,tmp);

  tmp=gsl_matrix_complex_get(iprop6,29,17);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,29,17,tmp);

  tmp=gsl_matrix_complex_get(iprop6,30,18);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,30,18,tmp);

  tmp=gsl_matrix_complex_get(iprop6,33,21);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,33,21,tmp);

  tmp=gsl_matrix_complex_get(iprop6,34,22);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,34,22,tmp);
  
  // Set lower-left derivatives of off-diagonal mass-like terms

  dipdgapu(12,0)=da;
  dipdgapu(15,3)=da;
  dipdgapu(16,4)=da;
  dipdgapu(19,7)=da;
  dipdgapu(20,8)=da;
  dipdgapu(23,11)=da;
  
  dipdgapd(12,0)=ua;
  dipdgapd(15,3)=ua;
  dipdgapd(16,4)=ua;
  dipdgapd(19,7)=ua;
  dipdgapd(20,8)=ua;
  dipdgapd(23,11)=ua;
  
  dipdgapu(13,1)=-da;
  dipdgapu(14,2)=-da;
  dipdgapu(17,5)=-da;
  dipdgapu(18,6)=-da;
  dipdgapu(21,9)=-da;
  dipdgapu(22,10)=-da;
  
  dipdgapd(13,1)=-ua;
  dipdgapd(14,2)=-ua;
  dipdgapd(17,5)=-ua;
  dipdgapd(18,6)=-ua;
  dipdgapd(21,9)=-ua;
  dipdgapd(22,10)=-ua;
  
  dipdgapu(24,0)=sa;
  dipdgapu(27,3)=sa;
  dipdgapu(28,4)=sa;
  dipdgapu(31,7)=sa;
  dipdgapu(32,8)=sa;
  dipdgapu(35,11)=sa;
  
  dipdgaps(24,0)=ua;
  dipdgaps(27,3)=ua;
  dipdgaps(28,4)=ua;
  dipdgaps(31,7)=ua;
  dipdgaps(32,8)=ua;
  dipdgaps(35,11)=ua;
  
  dipdgapu(25,1)=-sa;
  dipdgapu(26,2)=-sa;
  dipdgapu(29,5)=-sa;
  dipdgapu(30,6)=-sa;
  dipdgapu(33,9)=-sa;
  dipdgapu(34,10)=-sa;
  
  dipdgaps(25,1)=-ua;
  dipdgaps(26,2)=-ua;
  dipdgaps(29,5)=-ua;
  dipdgaps(30,6)=-ua;
  dipdgaps(33,9)=-ua;
  dipdgaps(34,10)=-ua;
  
  dipdgapd(24,12)=sa;
  dipdgapd(27,15)=sa;
  dipdgapd(28,16)=sa;
  dipdgapd(31,19)=sa;
  dipdgapd(32,20)=sa;
  dipdgapd(35,23)=sa;
  
  dipdgaps(24,12)=da;
  dipdgaps(27,15)=da;
  dipdgaps(28,16)=da;
  dipdgaps(31,19)=da;
  dipdgaps(32,20)=da;
  dipdgaps(35,23)=da;
  
  dipdgapd(25,13)=-sa;
  dipdgapd(26,14)=-sa;
  dipdgapd(29,17)=-sa;
  dipdgapd(30,18)=-sa;
  dipdgapd(33,21)=-sa;
  dipdgapd(34,22)=-sa;
  
  dipdgaps(25,13)=-da;
  dipdgaps(26,14)=-da;
  dipdgaps(29,17)=-da;
  dipdgaps(30,18)=-da;
  dipdgaps(33,21)=-da;
  dipdgaps(34,22)=-da;
  
  // Set upper-right derivatives of off-diagonal mass-like terms

  dipdgapu(0,12)=da;
  dipdgapu(3,15)=da;
  dipdgapu(4,16)=da;
  dipdgapu(7,19)=da;
  dipdgapu(8,20)=da;
  dipdgapu(11,23)=da;
  
  dipdgapd(0,12)=ua;
  dipdgapd(3,15)=ua;
  dipdgapd(4,16)=ua;
  dipdgapd(7,19)=ua;
  dipdgapd(8,20)=ua;
  dipdgapd(11,23)=ua;
  
  dipdgapu(1,13)=-da;
  dipdgapu(2,14)=-da;
  dipdgapu(5,17)=-da;
  dipdgapu(6,18)=-da;
  dipdgapu(9,21)=-da;
  dipdgapu(10,22)=-da;
  
  dipdgapd(1,13)=-ua;
  dipdgapd(2,14)=-ua;
  dipdgapd(5,17)=-ua;
  dipdgapd(6,18)=-ua;
  dipdgapd(9,21)=-ua;
  dipdgapd(10,22)=-ua;
  
  dipdgapu(0,24)=sa;
  dipdgapu(3,27)=sa;
  dipdgapu(4,28)=sa;
  dipdgapu(7,31)=sa;
  dipdgapu(8,32)=sa;
  dipdgapu(11,35)=sa;
  
  dipdgaps(0,24)=ua;
  dipdgaps(3,27)=ua;
  dipdgaps(4,28)=ua;
  dipdgaps(7,31)=ua;
  dipdgaps(8,32)=ua;
  dipdgaps(11,35)=ua;
  
  dipdgapu(1,25)=-sa;
  dipdgapu(2,26)=-sa;
  dipdgapu(5,29)=-sa;
  dipdgapu(6,30)=-sa;
  dipdgapu(9,33)=-sa;
  dipdgapu(10,34)=-sa;
  
  dipdgaps(1,25)=-ua;
  dipdgaps(2,26)=-ua;
  dipdgaps(5,29)=-ua;
  dipdgaps(6,30)=-ua;
  dipdgaps(9,33)=-ua;
  dipdgaps(10,34)=-ua;
  
  dipdgapd(12,24)=sa;
  dipdgapd(15,27)=sa;
  dipdgapd(16,28)=sa;
  dipdgapd(19,31)=sa;
  dipdgapd(20,32)=sa;
  dipdgapd(23,35)=sa;
  
  dipdgaps(12,24)=da;
  dipdgaps(15,27)=da;
  dipdgaps(16,28)=da;
  dipdgaps(19,31)=da;
  dipdgaps(20,32)=da;
  dipdgaps(23,35)=da;
  
  dipdgapd(13,25)=-sa;
  dipdgapd(14,26)=-sa;
  dipdgapd(17,29)=-sa;
  dipdgapd(18,30)=-sa;
  dipdgapd(21,33)=-sa;
  dipdgapd(22,34)=-sa;
  
  dipdgaps(13,25)=-da;
  dipdgaps(14,26)=-da;
  dipdgaps(17,29)=-da;
  dipdgaps(18,30)=-da;
  dipdgaps(21,33)=-da;
  dipdgaps(22,34)=-da;
  
  gsl_eigen_hermv(iprop6,eval6,eivec6,w6);
  
  double dmudqqu=-4.0*G;
  double dmudqqd=2.0*K*strange->qq;
  double dmudqqs=2.0*K*down->qq;
  double dmddqqu=2.0*K*strange->qq;
  double dmddqqd=-4.0*G;
  double dmddqqs=2.0*K*up->qq;
  double dmsdqqu=2.0*K*down->qq;
  double dmsdqqd=2.0*K*up->qq;
  double dmsdqqs=-4.0*G;
  double dmudgu=KD/2.0/GD/GD*up->del;
  double dmddgd=KD/2.0/GD/GD*down->del;
  double dmsdgs=KD/2.0/GD/GD*strange->del;
  
  for(k=0;k<36;k++) {
    egv[k]=gsl_vector_get(eval6,k);

    // Select the eigenvector and compute its conjugate
    ubvector_complex eicol(36), eiconj(36);
    for(size_t ij=0;ij<36;ij++) {
      gsl_complex tmp2;
      tmp2=gsl_matrix_complex_get(eivec6,ij,k);
      eicol[ij]=std::complex<double>(GSL_REAL(tmp2),GSL_IMAG(tmp2));
      eiconj[ij]=conj(eicol[ij]);
    }
    
    // Compute v+ M v
    std::complex<double> dp;
    std::complex<double> dedmu, dedmd, dedms;
    ubvector_complex tmpe;
    std::complex<double> zero3(0.0,0.0);
    
    // Contribution from explicit gap terms
    dp=zero3;
    tmpe=prod(dipdgapu,eicol);
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*tmpe[i];
    }
    dedu[k]=dp.real();
    dp=zero3;
    tmpe=prod(dipdgapd,eicol);
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*tmpe[i];
    }
    dedd[k]=dp.real();
    dp=zero3;
    tmpe=prod(dipdgaps,eicol);
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*tmpe[i];
    }
    deds[k]=dp.real();
    
    // Contribution from explicit <qq> entries
    dp=zero3;
    tmpe=prod(dipdqqu,eicol);
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*tmpe[i];
    }
    dedqqu[k]=dp.real();
    dp=zero3;
    tmpe=prod(dipdqqd,eicol);
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*tmpe[i];
    }
    dedqqd[k]=dp.real();
    dp=zero3;
    tmpe=prod(dipdqqs,eicol);
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*tmpe[i];
    }
    dedqqs[k]=dp.real();

    // mass derivatives
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mmu[i]*eicol[i];
    }
    dedmu=dp.real();
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mmd[i]*eicol[i];
    }
    dedmd=dp.real();
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mms[i]*eicol[i];
    }
    dedms=dp.real();

    // Contribution from mass terms to quark condensate derivatives
    dedqqu[k]+=dedmu.real()*dmudqqu+
      dedmd.real()*dmddqqu+dedms.real()*dmsdqqu;
    dedqqd[k]+=dedmu.real()*dmudqqd+
      dedmd.real()*dmddqqd+dedms.real()*dmsdqqd;
    dedqqs[k]+=dedmu.real()*dmudqqs+
      dedmd.real()*dmddqqs+dedms.real()*dmsdqqs;

    // Contribution from mass terms to gap derivatives
    dedu[k]+=dedmu.real()*dmudgu;
    dedd[k]+=dedmd.real()*dmddgd;
    deds[k]+=dedms.real()*dmsdgs;

    // Flavor chemical potential derivatives
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mmuu[i]*eicol[i];
    }
    dedmuu[k]=dp.real();
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mmud[i]*eicol[i];
    }
    dedmud[k]=dp.real();
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mmus[i]*eicol[i];
    }
    dedmus[k]=dp.real();
    
    // Color chemical potential derivatives
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mmur[i]*eicol[i];
    }
    double dedmur=dp.real();
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mmug[i]*eicol[i];
    }
    double dedmug=dp.real();
    dp=zero3;
    for(size_t i=0;i<36;i++) {
      dp+=eiconj[i]*mmub[i]*eicol[i];
    }
    double dedmub=dp.real();

    dedmu3[k]=dedmur-dedmug;
    dedmu8[k]=(dedmur+dedmug-2.0*dedmub)/sqrt(3.0);

  } 
  
  return 0;
}

int cfl6_eos::make_matrices(double lmom, double mu3, 
			     double mu8, double egv[36],
			     double dedmuu[36], double dedmud[36],
			     double dedmus[36], double dedqqu[36], 
			     double dedqqd[36], double dedqqs[36],
			     double dedu[36], double dedd[36],
			     double deds[36], double dedmu3[36],
			     double dedmu8[36]) {
  
  int k;
  const double mu=up->ms, md=down->ms, ms=strange->ms;
  const double muu=up->mu, mud=down->mu, mus=strange->mu;
  const double du=up->del, dd=down->del, ds=strange->del;
  const double qqu=up->qq, qqd=down->qq, qqs=strange->qq;

  double mur, mug, mub;
  mur=mu3+mu8/sqrt(3.0);
  mug=-mu3+mu8/sqrt(3.0);
  mub=-2.0*mu8/sqrt(3.0);
  
  for(k=0;k<36;k++) {
    egv[k]=0.0;
    dedmuu[k]=0.0;
    dedmud[k]=0.0;
    dedmus[k]=0.0;
    dedqqu[k]=0.0;
    dedqqd[k]=0.0;
    dedqqs[k]=0.0;
    dedu[k]=0.0;
    dedd[k]=0.0;
    deds[k]=0.0;
    dedmu3[k]=0.0;
    dedmu8[k]=0.0;
  }
  
  ubvector mmuu(mat_size), mmud(mat_size), mmus(mat_size);
  ubvector mmu(mat_size), mmd(mat_size), mms(mat_size);
  ubvector mmur(mat_size), mmug(mat_size), mmub(mat_size);
  for(k=0;k<36;k++) {
    mmuu[k]=0.0;
    mmud[k]=0.0;
    mmus[k]=0.0;
    mmu[k]=0.0;
    mmd[k]=0.0;
    mms[k]=0.0;
    mmur[k]=0.0;
    mmug[k]=0.0;
    mmub[k]=0.0;
  }

  gsl_complex zero={{0.0,0.0}};
  std::complex<double> zero2(0.0,0.0);
  for(int i=0;i<mat_size;i++) {
    for(int j=0;j<mat_size;j++) {
      gsl_matrix_complex_set(iprop6,i,j,zero);
      dipdgapu(i,j)=zero2;
      dipdgapd(i,j)=zero2;
      dipdgaps(i,j)=zero2;
      dipdqqu(i,j)=zero2;
      dipdqqd(i,j)=zero2;
      dipdqqs(i,j)=zero2;
    }
  }
  
  /// Set diagonal elements

  gsl_complex tmp={{0.0,0.0}};
  
  tmp=gsl_matrix_complex_get(iprop6,0,0);
  GSL_REAL(tmp)=-mu-muu-mur;
  gsl_matrix_complex_set(iprop6,0,0,tmp);

  tmp=gsl_matrix_complex_get(iprop6,1,1);
  GSL_REAL(tmp)=mu-muu-mur;
  gsl_matrix_complex_set(iprop6,1,1,tmp);

  tmp=gsl_matrix_complex_get(iprop6,2,2);
  GSL_REAL(tmp)=-mu+muu+mur;
  gsl_matrix_complex_set(iprop6,2,2,tmp);

  tmp=gsl_matrix_complex_get(iprop6,3,3);
  GSL_REAL(tmp)=mu+muu+mur;
  gsl_matrix_complex_set(iprop6,3,3,tmp);

  mmuu[0]=-1.0;
  mmuu[1]=-1.0;
  mmuu[2]=1.0;
  mmuu[3]=1.0;
  mmur[0]=-1.0;
  mmur[1]=-1.0;
  mmur[2]=1.0;
  mmur[3]=1.0;
  mmu[0]=-1.0;
  mmu[1]=1.0;
  mmu[2]=-1.0;
  mmu[3]=1.0;

  tmp=gsl_matrix_complex_get(iprop6,4,4);
  GSL_REAL(tmp)=-mu-muu-mug;
  gsl_matrix_complex_set(iprop6,4,4,tmp);

  tmp=gsl_matrix_complex_get(iprop6,5,5);
  GSL_REAL(tmp)=mu-muu-mug;
  gsl_matrix_complex_set(iprop6,5,5,tmp);

  tmp=gsl_matrix_complex_get(iprop6,6,6);
  GSL_REAL(tmp)=-mu+muu+mug;
  gsl_matrix_complex_set(iprop6,6,6,tmp);

  tmp=gsl_matrix_complex_get(iprop6,7,7);
  GSL_REAL(tmp)=mu+muu+mug;
  gsl_matrix_complex_set(iprop6,7,7,tmp);

  mmuu[4]=-1.0;
  mmuu[5]=-1.0;
  mmuu[6]=1.0;
  mmuu[7]=1.0;
  mmug[4]=-1.0;
  mmug[5]=-1.0;
  mmug[6]=1.0;
  mmug[7]=1.0;
  mmu[4]=-1.0;
  mmu[5]=1.0;
  mmu[6]=-1.0;
  mmu[7]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,8,8);
  GSL_REAL(tmp)=-mu-muu-mub;
  gsl_matrix_complex_set(iprop6,8,8,tmp);

  tmp=gsl_matrix_complex_get(iprop6,9,9);
  GSL_REAL(tmp)=mu-muu-mub;
  gsl_matrix_complex_set(iprop6,9,9,tmp);

  tmp=gsl_matrix_complex_get(iprop6,10,10);
  GSL_REAL(tmp)=-mu+muu+mub;
  gsl_matrix_complex_set(iprop6,10,10,tmp);

  tmp=gsl_matrix_complex_get(iprop6,11,11);
  GSL_REAL(tmp)=mu+muu+mub;
  gsl_matrix_complex_set(iprop6,11,11,tmp);

  mmuu[8]=-1.0;
  mmuu[9]=-1.0;
  mmuu[10]=1.0;
  mmuu[11]=1.0;
  mmub[8]=-1.0;
  mmub[9]=-1.0;
  mmub[10]=1.0;
  mmub[11]=1.0;
  mmu[8]=-1.0;
  mmu[9]=1.0;
  mmu[10]=-1.0;
  mmu[11]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,12,12);
  GSL_REAL(tmp)=-md-mud-mur;
  gsl_matrix_complex_set(iprop6,12,12,tmp);

  tmp=gsl_matrix_complex_get(iprop6,13,13);
  GSL_REAL(tmp)=md-mud-mur;
  gsl_matrix_complex_set(iprop6,13,13,tmp);

  tmp=gsl_matrix_complex_get(iprop6,14,14);
  GSL_REAL(tmp)=-md+mud+mur;
  gsl_matrix_complex_set(iprop6,14,14,tmp);

  tmp=gsl_matrix_complex_get(iprop6,15,15);
  GSL_REAL(tmp)=md+mud+mur;
  gsl_matrix_complex_set(iprop6,15,15,tmp);

  mmud[12]=-1.0;
  mmud[13]=-1.0;
  mmud[14]=1.0;
  mmud[15]=1.0;
  mmur[12]=-1.0;
  mmur[13]=-1.0;
  mmur[14]=1.0;
  mmur[15]=1.0;
  mmd[12]=-1.0;
  mmd[13]=1.0;
  mmd[14]=-1.0;
  mmd[15]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,16,16);
  GSL_REAL(tmp)=-md-mud-mug;
  gsl_matrix_complex_set(iprop6,16,16,tmp);

  tmp=gsl_matrix_complex_get(iprop6,17,17);
  GSL_REAL(tmp)=md-mud-mug;
  gsl_matrix_complex_set(iprop6,17,17,tmp);

  tmp=gsl_matrix_complex_get(iprop6,18,18);
  GSL_REAL(tmp)=-md+mud+mug;
  gsl_matrix_complex_set(iprop6,18,18,tmp);

  tmp=gsl_matrix_complex_get(iprop6,19,19);
  GSL_REAL(tmp)=md+mud+mug;
  gsl_matrix_complex_set(iprop6,19,19,tmp);
  
  mmud[16]=-1.0;
  mmud[17]=-1.0;
  mmud[18]=1.0;
  mmud[19]=1.0;
  mmug[16]=-1.0;
  mmug[17]=-1.0;
  mmug[18]=1.0;
  mmug[19]=1.0;
  mmd[16]=-1.0;
  mmd[17]=1.0;
  mmd[18]=-1.0;
  mmd[19]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,20,20);
  GSL_REAL(tmp)=-md-mud-mub;
  gsl_matrix_complex_set(iprop6,20,20,tmp);

  tmp=gsl_matrix_complex_get(iprop6,21,21);
  GSL_REAL(tmp)=md-mud-mub;
  gsl_matrix_complex_set(iprop6,21,21,tmp);

  tmp=gsl_matrix_complex_get(iprop6,22,22);
  GSL_REAL(tmp)=-md+mud+mub;
  gsl_matrix_complex_set(iprop6,22,22,tmp);

  tmp=gsl_matrix_complex_get(iprop6,23,23);
  GSL_REAL(tmp)=md+mud+mub;
  gsl_matrix_complex_set(iprop6,23,23,tmp);

  mmud[20]=-1.0;
  mmud[21]=-1.0;
  mmud[22]=1.0;
  mmud[23]=1.0;
  mmub[20]=-1.0;
  mmub[21]=-1.0;
  mmub[22]=1.0;
  mmub[23]=1.0;
  mmd[20]=-1.0;
  mmd[21]=1.0;
  mmd[22]=-1.0;
  mmd[23]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,24,24);
  GSL_REAL(tmp)=-ms-mus-mur;
  gsl_matrix_complex_set(iprop6,24,24,tmp);

  tmp=gsl_matrix_complex_get(iprop6,25,25);
  GSL_REAL(tmp)=ms-mus-mur;
  gsl_matrix_complex_set(iprop6,25,25,tmp);

  tmp=gsl_matrix_complex_get(iprop6,26,26);
  GSL_REAL(tmp)=-ms+mus+mur;
  gsl_matrix_complex_set(iprop6,26,26,tmp);

  tmp=gsl_matrix_complex_get(iprop6,27,27);
  GSL_REAL(tmp)=ms+mus+mur;
  gsl_matrix_complex_set(iprop6,27,27,tmp);

  mmus[24]=-1.0;
  mmus[25]=-1.0;
  mmus[26]=1.0;
  mmus[27]=1.0;
  mmur[24]=-1.0;
  mmur[25]=-1.0;
  mmur[26]=1.0;
  mmur[27]=1.0;
  mms[24]=-1.0;
  mms[25]=1.0;
  mms[26]=-1.0;
  mms[27]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,28,28);
  GSL_REAL(tmp)=-ms-mus-mug;
  gsl_matrix_complex_set(iprop6,28,28,tmp);

  tmp=gsl_matrix_complex_get(iprop6,29,29);
  GSL_REAL(tmp)=ms-mus-mug;
  gsl_matrix_complex_set(iprop6,29,29,tmp);

  tmp=gsl_matrix_complex_get(iprop6,30,30);
  GSL_REAL(tmp)=-ms+mus+mug;
  gsl_matrix_complex_set(iprop6,30,30,tmp);

  tmp=gsl_matrix_complex_get(iprop6,31,31);
  GSL_REAL(tmp)=ms+mus+mug;
  gsl_matrix_complex_set(iprop6,31,31,tmp);

  mmus[28]=-1.0;
  mmus[29]=-1.0;
  mmus[30]=1.0;
  mmus[31]=1.0;
  mmug[28]=-1.0;
  mmug[29]=-1.0;
  mmug[30]=1.0;
  mmug[31]=1.0;
  mms[28]=-1.0;
  mms[29]=1.0;
  mms[30]=-1.0;
  mms[31]=1.0;
  
  tmp=gsl_matrix_complex_get(iprop6,32,32);
  GSL_REAL(tmp)=-ms-mus-mub;
  gsl_matrix_complex_set(iprop6,32,32,tmp);

  tmp=gsl_matrix_complex_get(iprop6,33,33);
  GSL_REAL(tmp)=ms-mus-mub;
  gsl_matrix_complex_set(iprop6,33,33,tmp);

  tmp=gsl_matrix_complex_get(iprop6,34,34);
  GSL_REAL(tmp)=-ms+mus+mub;
  gsl_matrix_complex_set(iprop6,34,34,tmp);

  tmp=gsl_matrix_complex_get(iprop6,35,35);
  GSL_REAL(tmp)=ms+mus+mub;
  gsl_matrix_complex_set(iprop6,35,35,tmp);

  mmus[32]=-1.0;
  mmus[33]=-1.0;
  mmus[34]=1.0;
  mmus[35]=1.0;
  mmub[32]=-1.0;
  mmub[33]=-1.0;
  mmub[34]=1.0;
  mmub[35]=1.0;
  mms[32]=-1.0;
  mms[33]=1.0;
  mms[34]=-1.0;
  mms[35]=1.0;
  
  // Set momentum
  
  for(int kk=0;kk<35;kk+=2) {
    tmp=gsl_matrix_complex_get(iprop6,kk+1,kk);
    GSL_REAL(tmp)=-lmom;
    gsl_matrix_complex_set(iprop6,kk+1,kk,tmp);
  }

  GSL_REAL(tmp)=0.0;
  
  // Set gaps

  if (fixed_mass) {

    // Set gap terms

    tmp=gsl_matrix_complex_get(iprop6,28,23);
    GSL_IMAG(tmp)=du;
    gsl_matrix_complex_set(iprop6,28,23,tmp);

    tmp=gsl_matrix_complex_get(iprop6,30,21);
    GSL_IMAG(tmp)=du;
    gsl_matrix_complex_set(iprop6,30,21,tmp);

    tmp=gsl_matrix_complex_get(iprop6,33,18);
    GSL_IMAG(tmp)=du;
    gsl_matrix_complex_set(iprop6,33,18,tmp);

    tmp=gsl_matrix_complex_get(iprop6,35,16);
    GSL_IMAG(tmp)=du;
    gsl_matrix_complex_set(iprop6,35,16,tmp);

    // --

    tmp=gsl_matrix_complex_get(iprop6,29,22);
    GSL_IMAG(tmp)=-du;
    gsl_matrix_complex_set(iprop6,29,22,tmp);

    tmp=gsl_matrix_complex_get(iprop6,31,20);
    GSL_IMAG(tmp)=-du;
    gsl_matrix_complex_set(iprop6,31,20,tmp);

    tmp=gsl_matrix_complex_get(iprop6,32,19);
    GSL_IMAG(tmp)=-du;
    gsl_matrix_complex_set(iprop6,32,19,tmp);

    tmp=gsl_matrix_complex_get(iprop6,34,17);
    GSL_IMAG(tmp)=-du;
    gsl_matrix_complex_set(iprop6,34,17,tmp);

    // --

    tmp=gsl_matrix_complex_get(iprop6,24,11);
    GSL_IMAG(tmp)=dd;
    gsl_matrix_complex_set(iprop6,24,11,tmp);

    tmp=gsl_matrix_complex_get(iprop6,26,9);
    GSL_IMAG(tmp)=dd;
    gsl_matrix_complex_set(iprop6,26,9,tmp);

    tmp=gsl_matrix_complex_get(iprop6,33,2);
    GSL_IMAG(tmp)=dd;
    gsl_matrix_complex_set(iprop6,33,2,tmp);

    tmp=gsl_matrix_complex_get(iprop6,35,0);
    GSL_IMAG(tmp)=dd;
    gsl_matrix_complex_set(iprop6,35,0,tmp);

    // --			   

    tmp=gsl_matrix_complex_get(iprop6,25,10);
    GSL_IMAG(tmp)=-dd;
    gsl_matrix_complex_set(iprop6,25,10,tmp);

    tmp=gsl_matrix_complex_get(iprop6,27,8);
    GSL_IMAG(tmp)=-dd;
    gsl_matrix_complex_set(iprop6,27,8,tmp);

    tmp=gsl_matrix_complex_get(iprop6,32,3);
    GSL_IMAG(tmp)=-dd;
    gsl_matrix_complex_set(iprop6,32,3,tmp);

    tmp=gsl_matrix_complex_get(iprop6,34,1);
    GSL_IMAG(tmp)=-dd;
    gsl_matrix_complex_set(iprop6,34,1,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,12,7);
    GSL_IMAG(tmp)=ds;
    gsl_matrix_complex_set(iprop6,12,7,tmp);

    tmp=gsl_matrix_complex_get(iprop6,14,5);
    GSL_IMAG(tmp)=ds;
    gsl_matrix_complex_set(iprop6,14,5,tmp);

    tmp=gsl_matrix_complex_get(iprop6,17,2);
    GSL_IMAG(tmp)=ds;
    gsl_matrix_complex_set(iprop6,17,2,tmp);

    tmp=gsl_matrix_complex_get(iprop6,19,0);
    GSL_IMAG(tmp)=ds;
    gsl_matrix_complex_set(iprop6,19,0,tmp);

    // --

    tmp=gsl_matrix_complex_get(iprop6,13,6);
    GSL_IMAG(tmp)=-ds;
    gsl_matrix_complex_set(iprop6,13,6,tmp);

    tmp=gsl_matrix_complex_get(iprop6,15,4);
    GSL_IMAG(tmp)=-ds;
    gsl_matrix_complex_set(iprop6,15,4,tmp);

    tmp=gsl_matrix_complex_get(iprop6,16,3);
    GSL_IMAG(tmp)=-ds;
    gsl_matrix_complex_set(iprop6,16,3,tmp);

    tmp=gsl_matrix_complex_get(iprop6,18,1);
    GSL_IMAG(tmp)=-ds;
    gsl_matrix_complex_set(iprop6,18,1,tmp);
    
    // Set lower-left derivatives of gap terms

    std::complex<double> littlei(0.0,1.0);
    std::complex<double> mlittlei(0.0,-1.0);

    dipdgapu(28,23)=littlei;
    dipdgapu(30,21)=littlei;
    dipdgapu(33,18)=littlei;
    dipdgapu(35,16)=littlei;

    dipdgapu(29,22)=mlittlei;
    dipdgapu(31,20)=mlittlei;
    dipdgapu(32,19)=mlittlei;
    dipdgapu(34,17)=mlittlei;
    
    dipdgapd(24,11)=littlei;
    dipdgapd(26,9)=littlei;
    dipdgapd(33,2)=littlei;
    dipdgapd(35,0)=littlei;

    dipdgapd(25,10)=mlittlei;
    dipdgapd(27,8)=mlittlei;
    dipdgapd(32,3)=mlittlei;
    dipdgapd(34,1)=mlittlei;
    
    dipdgaps(12,7)=littlei;
    dipdgaps(14,5)=littlei;
    dipdgaps(17,2)=littlei;
    dipdgaps(19,0)=littlei;

    dipdgaps(13,6)=mlittlei;
    dipdgaps(15,4)=mlittlei;
    dipdgaps(16,3)=mlittlei;
    dipdgaps(18,1)=mlittlei;

    // Set upper-right derivatives of gap terms

    dipdgapu(23,28)=mlittlei;
    dipdgapu(21,30)=mlittlei;
    dipdgapu(18,33)=mlittlei;
    dipdgapu(16,35)=mlittlei;
    
    dipdgapu(22,29)=littlei;
    dipdgapu(20,31)=littlei;
    dipdgapu(19,32)=littlei;
    dipdgapu(17,34)=littlei;
    
    dipdgapd(11,24)=mlittlei;
    dipdgapd(9,26)=mlittlei;
    dipdgapd(2,33)=mlittlei;
    dipdgapd(0,35)=mlittlei;
    
    dipdgapd(10,25)=littlei;
    dipdgapd(8,27)=littlei;
    dipdgapd(3,32)=littlei;
    dipdgapd(1,34)=littlei;
    
    dipdgaps(7,12)=mlittlei;
    dipdgaps(5,14)=mlittlei;
    dipdgaps(2,17)=mlittlei;
    dipdgaps(0,19)=mlittlei;
    
    dipdgaps(6,13)=littlei;
    dipdgaps(4,15)=littlei;
    dipdgaps(3,16)=littlei;
    dipdgaps(1,18)=littlei;
    
  } else {

    // set gap terms
    double betau=1.0+KD/GD*qqu/3.0;
    double betad=1.0+KD/GD*qqd/3.0;
    double betas=1.0+KD/GD*qqs/3.0;
    double gammau=du*KD/GD/3.0;
    double gammad=dd*KD/GD/3.0;
    double gammas=ds*KD/GD/3.0;
    
    tmp=gsl_matrix_complex_get(iprop6,28,23);
    GSL_IMAG(tmp)=du*betau;
    gsl_matrix_complex_set(iprop6,28,23,tmp);

    tmp=gsl_matrix_complex_get(iprop6,30,21);
    GSL_IMAG(tmp)=du*betau;
    gsl_matrix_complex_set(iprop6,30,21,tmp);

    tmp=gsl_matrix_complex_get(iprop6,33,18);
    GSL_IMAG(tmp)=du*betau;
    gsl_matrix_complex_set(iprop6,33,18,tmp);

    tmp=gsl_matrix_complex_get(iprop6,35,16);
    GSL_IMAG(tmp)=du*betau;
    gsl_matrix_complex_set(iprop6,35,16,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,29,22);
    GSL_IMAG(tmp)=-du*betau;
    gsl_matrix_complex_set(iprop6,29,22,tmp);

    tmp=gsl_matrix_complex_get(iprop6,31,20);
    GSL_IMAG(tmp)=-du*betau;
    gsl_matrix_complex_set(iprop6,31,20,tmp);

    tmp=gsl_matrix_complex_get(iprop6,32,19);
    GSL_IMAG(tmp)=-du*betau;
    gsl_matrix_complex_set(iprop6,32,19,tmp);

    tmp=gsl_matrix_complex_get(iprop6,34,17);
    GSL_IMAG(tmp)=-du*betau;
    gsl_matrix_complex_set(iprop6,34,17,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,24,11);
    GSL_IMAG(tmp)=dd*betad;
    gsl_matrix_complex_set(iprop6,24,11,tmp);

    tmp=gsl_matrix_complex_get(iprop6,26,9);
    GSL_IMAG(tmp)=dd*betad;
    gsl_matrix_complex_set(iprop6,26,9,tmp);

    tmp=gsl_matrix_complex_get(iprop6,33,2);
    GSL_IMAG(tmp)=dd*betad;
    gsl_matrix_complex_set(iprop6,33,2,tmp);

    tmp=gsl_matrix_complex_get(iprop6,35,0);
    GSL_IMAG(tmp)=dd*betad;
    gsl_matrix_complex_set(iprop6,35,0,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,25,10);
    GSL_IMAG(tmp)=-dd*betad;
    gsl_matrix_complex_set(iprop6,25,10,tmp);

    tmp=gsl_matrix_complex_get(iprop6,27,8);
    GSL_IMAG(tmp)=-dd*betad;
    gsl_matrix_complex_set(iprop6,27,8,tmp);

    tmp=gsl_matrix_complex_get(iprop6,32,3);
    GSL_IMAG(tmp)=-dd*betad;
    gsl_matrix_complex_set(iprop6,32,3,tmp);
    
    tmp=gsl_matrix_complex_get(iprop6,34,1);
    GSL_IMAG(tmp)=-dd*betad;
    gsl_matrix_complex_set(iprop6,34,1,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,12,7);
    GSL_IMAG(tmp)=ds*betas;
    gsl_matrix_complex_set(iprop6,12,7,tmp);

    tmp=gsl_matrix_complex_get(iprop6,14,5);
    GSL_IMAG(tmp)=ds*betas;
    gsl_matrix_complex_set(iprop6,14,5,tmp);

    tmp=gsl_matrix_complex_get(iprop6,17,2);
    GSL_IMAG(tmp)=ds*betas;
    gsl_matrix_complex_set(iprop6,17,2,tmp);

    tmp=gsl_matrix_complex_get(iprop6,19,0);
    GSL_IMAG(tmp)=ds*betas;
    gsl_matrix_complex_set(iprop6,19,0,tmp);

    // --
    
    tmp=gsl_matrix_complex_get(iprop6,13,6);
    GSL_IMAG(tmp)=-ds*betas;
    gsl_matrix_complex_set(iprop6,13,6,tmp);
    
    tmp=gsl_matrix_complex_get(iprop6,15,4);
    GSL_IMAG(tmp)=-ds*betas;
    gsl_matrix_complex_set(iprop6,15,4,tmp);

    tmp=gsl_matrix_complex_get(iprop6,16,3);
    GSL_IMAG(tmp)=-ds*betas;
    gsl_matrix_complex_set(iprop6,16,3,tmp);

    tmp=gsl_matrix_complex_get(iprop6,18,1);
    GSL_IMAG(tmp)=-ds*betas;
    gsl_matrix_complex_set(iprop6,18,1,tmp);
    
    // Set lower-left derivatives of gap terms

    std::complex<double> littlei(0.0,1.0);

    dipdgapu(28,23)=littlei*betau;
    dipdgapu(30,21)=littlei*betau;
    dipdgapu(33,18)=littlei*betau;
    dipdgapu(35,16)=littlei*betau;

    dipdgapu(29,22)=littlei*-betau;
    dipdgapu(31,20)=littlei*-betau;
    dipdgapu(32,19)=littlei*-betau;
    dipdgapu(34,17)=littlei*-betau;
    
    dipdgapd(24,11)=littlei*betad;
    dipdgapd(26,9)=littlei*betad;
    dipdgapd(33,2)=littlei*betad;
    dipdgapd(35,0)=littlei*betad;

    dipdgapd(25,10)=littlei*-betad;
    dipdgapd(27,8)=littlei*-betad;
    dipdgapd(32,3)=littlei*-betad;
    dipdgapd(34,1)=littlei*-betad;
    
    dipdgaps(12,7)=littlei*betas;
    dipdgaps(14,5)=littlei*betas;
    dipdgaps(17,2)=littlei*betas;
    dipdgaps(19,0)=littlei*betas;
    
    dipdgaps(13,6)=littlei*-betas;
    dipdgaps(15,4)=littlei*-betas;
    dipdgaps(16,3)=littlei*-betas;
    dipdgaps(18,1)=littlei*-betas;

    dipdqqu(28,23)=littlei*gammau;
    dipdqqu(30,21)=littlei*gammau;
    dipdqqu(33,18)=littlei*gammau;
    dipdqqu(35,16)=littlei*gammau;

    dipdqqu(29,22)=littlei*-gammau;
    dipdqqu(31,20)=littlei*-gammau;
    dipdqqu(32,19)=littlei*-gammau;
    dipdqqu(34,17)=littlei*-gammau;
    
    dipdqqd(24,11)=littlei*gammad;
    dipdqqd(26,9)=littlei*gammad;
    dipdqqd(33,2)=littlei*gammad;
    dipdqqd(35,0)=littlei*gammad;

    dipdqqd(25,10)=littlei*-gammad;
    dipdqqd(27,8)=littlei*-gammad;
    dipdqqd(32,3)=littlei*-gammad;
    dipdqqd(34,1)=littlei*-gammad;
    
    dipdqqs(12,7)=littlei*gammas;
    dipdqqs(14,5)=littlei*gammas;
    dipdqqs(17,2)=littlei*gammas;
    dipdqqs(19,0)=littlei*gammas;
    
    dipdqqs(13,6)=littlei*-gammas;
    dipdqqs(15,4)=littlei*-gammas;
    dipdqqs(16,3)=littlei*-gammas;
    dipdqqs(18,1)=littlei*-gammas;

    // Set upper-right derivatives of gap terms

    dipdgapu(23,28)=littlei*-betau;
    dipdgapu(21,30)=littlei*-betau;
    dipdgapu(18,33)=littlei*-betau;
    dipdgapu(16,35)=littlei*-betau;
    
    dipdgapu(22,29)=littlei*betau;
    dipdgapu(20,31)=littlei*betau;
    dipdgapu(19,32)=littlei*betau;
    dipdgapu(17,34)=littlei*betau;
    
    dipdgapd(11,24)=littlei*-betad;
    dipdgapd(9,26)=littlei*-betad;
    dipdgapd(2,33)=littlei*-betad;
    dipdgapd(0,35)=littlei*-betad;
    
    dipdgapd(10,25)=littlei*betad;
    dipdgapd(8,27)=littlei*betad;
    dipdgapd(3,32)=littlei*betad;
    dipdgapd(1,34)=littlei*betad;
    
    dipdgaps(7,12)=littlei*-betas;
    dipdgaps(5,14)=littlei*-betas;
    dipdgaps(2,17)=littlei*-betas;
    dipdgaps(0,19)=littlei*-betas;
    
    dipdgaps(6,13)=littlei*betas;
    dipdgaps(4,15)=littlei*betas;
    dipdgaps(3,16)=littlei*betas;
    dipdgaps(1,18)=littlei*betas;

    dipdqqu(23,28)=littlei*-gammau;
    dipdqqu(21,30)=littlei*-gammau;
    dipdqqu(18,33)=littlei*-gammau;
    dipdqqu(16,35)=littlei*-gammau;
    
    dipdqqu(22,29)=littlei*gammau;
    dipdqqu(20,31)=littlei*gammau;
    dipdqqu(19,32)=littlei*gammau;
    dipdqqu(17,34)=littlei*gammau;
    
    dipdqqd(11,24)=littlei*-gammad;
    dipdqqd(9,26)=littlei*-gammad;
    dipdqqd(2,33)=littlei*-gammad;
    dipdqqd(0,35)=littlei*-gammad;
    
    dipdqqd(10,25)=littlei*gammad;
    dipdqqd(8,27)=littlei*gammad;
    dipdqqd(3,32)=littlei*gammad;
    dipdqqd(1,34)=littlei*gammad;
    
    dipdqqs(7,12)=littlei*-gammas;
    dipdqqs(5,14)=littlei*-gammas;
    dipdqqs(2,17)=littlei*-gammas;
    dipdqqs(0,19)=littlei*-gammas;
    
    dipdqqs(6,13)=littlei*gammas;
    dipdqqs(4,15)=littlei*gammas;
    dipdqqs(3,16)=littlei*gammas;
    dipdqqs(1,18)=littlei*gammas;

  }
  
  double acons=KD/4.0/GD/GD;
  double ua=du*acons;
  double da=dd*acons;
  double sa=ds*acons;
  double uda=du*dd*acons;
  double usa=du*ds*acons;
  double dsa=dd*ds*acons;

  // Set off-diagonal mass-like terms

  tmp=gsl_matrix_complex_get(iprop6,12,0);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,12,0,tmp);

  tmp=gsl_matrix_complex_get(iprop6,15,3);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,15,3,tmp);

  tmp=gsl_matrix_complex_get(iprop6,16,4);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,16,4,tmp);

  tmp=gsl_matrix_complex_get(iprop6,19,7);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,19,7,tmp);

  tmp=gsl_matrix_complex_get(iprop6,20,8);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,20,8,tmp);

  tmp=gsl_matrix_complex_get(iprop6,23,11);
  GSL_REAL(tmp)=uda;
  gsl_matrix_complex_set(iprop6,23,11,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,13,1);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,13,1,tmp);

  tmp=gsl_matrix_complex_get(iprop6,14,2);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,14,2,tmp);

  tmp=gsl_matrix_complex_get(iprop6,17,5);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,17,5,tmp);

  tmp=gsl_matrix_complex_get(iprop6,18,6);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,18,6,tmp);

  tmp=gsl_matrix_complex_get(iprop6,21,9);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,21,9,tmp);

  tmp=gsl_matrix_complex_get(iprop6,22,10);
  GSL_REAL(tmp)=-uda;
  gsl_matrix_complex_set(iprop6,22,10,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,24,0);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,24,0,tmp);

  tmp=gsl_matrix_complex_get(iprop6,27,3);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,27,3,tmp);

  tmp=gsl_matrix_complex_get(iprop6,28,4);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,28,4,tmp);

  tmp=gsl_matrix_complex_get(iprop6,31,7);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,31,7,tmp);

  tmp=gsl_matrix_complex_get(iprop6,32,8);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,32,8,tmp);

  tmp=gsl_matrix_complex_get(iprop6,35,11);
  GSL_REAL(tmp)=usa;
  gsl_matrix_complex_set(iprop6,35,11,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,25,1);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,25,1,tmp);

  tmp=gsl_matrix_complex_get(iprop6,26,2);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,26,2,tmp);

  tmp=gsl_matrix_complex_get(iprop6,29,5);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,29,5,tmp);

  tmp=gsl_matrix_complex_get(iprop6,30,6);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,30,6,tmp);

  tmp=gsl_matrix_complex_get(iprop6,33,9);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,33,9,tmp);

  tmp=gsl_matrix_complex_get(iprop6,34,10);
  GSL_REAL(tmp)=-usa;
  gsl_matrix_complex_set(iprop6,34,10,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,24,12);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,24,12,tmp);

  tmp=gsl_matrix_complex_get(iprop6,27,15);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,27,15,tmp);

  tmp=gsl_matrix_complex_get(iprop6,28,16);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,28,16,tmp);

  tmp=gsl_matrix_complex_get(iprop6,31,19);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,31,19,tmp);

  tmp=gsl_matrix_complex_get(iprop6,32,20);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,32,20,tmp);

  tmp=gsl_matrix_complex_get(iprop6,35,23);
  GSL_REAL(tmp)=dsa;
  gsl_matrix_complex_set(iprop6,35,23,tmp);

  
  tmp=gsl_matrix_complex_get(iprop6,25,13);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,25,13,tmp);

  tmp=gsl_matrix_complex_get(iprop6,26,14);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,26,14,tmp);

  tmp=gsl_matrix_complex_get(iprop6,29,17);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,29,17,tmp);

  tmp=gsl_matrix_complex_get(iprop6,30,18);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,30,18,tmp);

  tmp=gsl_matrix_complex_get(iprop6,33,21);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,33,21,tmp);

  tmp=gsl_matrix_complex_get(iprop6,34,22);
  GSL_REAL(tmp)=-dsa;
  gsl_matrix_complex_set(iprop6,34,22,tmp);

  // Set lower-left derivatives of off-diagonal mass-like terms

  dipdgapu(12,0)=da;
  dipdgapu(15,3)=da;
  dipdgapu(16,4)=da;
  dipdgapu(19,7)=da;
  dipdgapu(20,8)=da;
  dipdgapu(23,11)=da;
  
  dipdgapd(12,0)=ua;
  dipdgapd(15,3)=ua;
  dipdgapd(16,4)=ua;
  dipdgapd(19,7)=ua;
  dipdgapd(20,8)=ua;
  dipdgapd(23,11)=ua;
  
  dipdgapu(13,1)=-da;
  dipdgapu(14,2)=-da;
  dipdgapu(17,5)=-da;
  dipdgapu(18,6)=-da;
  dipdgapu(21,9)=-da;
  dipdgapu(22,10)=-da;
  
  dipdgapd(13,1)=-ua;
  dipdgapd(14,2)=-ua;
  dipdgapd(17,5)=-ua;
  dipdgapd(18,6)=-ua;
  dipdgapd(21,9)=-ua;
  dipdgapd(22,10)=-ua;
  
  dipdgapu(24,0)=sa;
  dipdgapu(27,3)=sa;
  dipdgapu(28,4)=sa;
  dipdgapu(31,7)=sa;
  dipdgapu(32,8)=sa;
  dipdgapu(35,11)=sa;
  
  dipdgaps(24,0)=ua;
  dipdgaps(27,3)=ua;
  dipdgaps(28,4)=ua;
  dipdgaps(31,7)=ua;
  dipdgaps(32,8)=ua;
  dipdgaps(35,11)=ua;
  
  dipdgapu(25,1)=-sa;
  dipdgapu(26,2)=-sa;
  dipdgapu(29,5)=-sa;
  dipdgapu(30,6)=-sa;
  dipdgapu(33,9)=-sa;
  dipdgapu(34,10)=-sa;
  
  dipdgaps(25,1)=-ua;
  dipdgaps(26,2)=-ua;
  dipdgaps(29,5)=-ua;
  dipdgaps(30,6)=-ua;
  dipdgaps(33,9)=-ua;
  dipdgaps(34,10)=-ua;
  
  dipdgapd(24,12)=sa;
  dipdgapd(27,15)=sa;
  dipdgapd(28,16)=sa;
  dipdgapd(31,19)=sa;
  dipdgapd(32,20)=sa;
  dipdgapd(35,23)=sa;
  
  dipdgaps(24,12)=da;
  dipdgaps(27,15)=da;
  dipdgaps(28,16)=da;
  dipdgaps(31,19)=da;
  dipdgaps(32,20)=da;
  dipdgaps(35,23)=da;
  
  dipdgapd(25,13)=-sa;
  dipdgapd(26,14)=-sa;
  dipdgapd(29,17)=-sa;
  dipdgapd(30,18)=-sa;
  dipdgapd(33,21)=-sa;
  dipdgapd(34,22)=-sa;
  
  dipdgaps(25,13)=-da;
  dipdgaps(26,14)=-da;
  dipdgaps(29,17)=-da;
  dipdgaps(30,18)=-da;
  dipdgaps(33,21)=-da;
  dipdgaps(34,22)=-da;
  
  // Set upper-right derivatives of off-diagonal mass-like terms

  dipdgapu(0,12)=da;
  dipdgapu(3,15)=da;
  dipdgapu(4,16)=da;
  dipdgapu(7,19)=da;
  dipdgapu(8,20)=da;
  dipdgapu(11,23)=da;
  
  dipdgapd(0,12)=ua;
  dipdgapd(3,15)=ua;
  dipdgapd(4,16)=ua;
  dipdgapd(7,19)=ua;
  dipdgapd(8,20)=ua;
  dipdgapd(11,23)=ua;
  
  dipdgapu(1,13)=-da;
  dipdgapu(2,14)=-da;
  dipdgapu(5,17)=-da;
  dipdgapu(6,18)=-da;
  dipdgapu(9,21)=-da;
  dipdgapu(10,22)=-da;
  
  dipdgapd(1,13)=-ua;
  dipdgapd(2,14)=-ua;
  dipdgapd(5,17)=-ua;
  dipdgapd(6,18)=-ua;
  dipdgapd(9,21)=-ua;
  dipdgapd(10,22)=-ua;
  
  dipdgapu(0,24)=sa;
  dipdgapu(3,27)=sa;
  dipdgapu(4,28)=sa;
  dipdgapu(7,31)=sa;
  dipdgapu(8,32)=sa;
  dipdgapu(11,35)=sa;
  
  dipdgaps(0,24)=ua;
  dipdgaps(3,27)=ua;
  dipdgaps(4,28)=ua;
  dipdgaps(7,31)=ua;
  dipdgaps(8,32)=ua;
  dipdgaps(11,35)=ua;
  
  dipdgapu(1,25)=-sa;
  dipdgapu(2,26)=-sa;
  dipdgapu(5,29)=-sa;
  dipdgapu(6,30)=-sa;
  dipdgapu(9,33)=-sa;
  dipdgapu(10,34)=-sa;
  
  dipdgaps(1,25)=-ua;
  dipdgaps(2,26)=-ua;
  dipdgaps(5,29)=-ua;
  dipdgaps(6,30)=-ua;
  dipdgaps(9,33)=-ua;
  dipdgaps(10,34)=-ua;
  
  dipdgapd(12,24)=sa;
  dipdgapd(15,27)=sa;
  dipdgapd(16,28)=sa;
  dipdgapd(19,31)=sa;
  dipdgapd(20,32)=sa;
  dipdgapd(23,35)=sa;
  
  dipdgaps(12,24)=da;
  dipdgaps(15,27)=da;
  dipdgaps(16,28)=da;
  dipdgaps(19,31)=da;
  dipdgaps(20,32)=da;
  dipdgaps(23,35)=da;
  
  dipdgapd(13,25)=-sa;
  dipdgapd(14,26)=-sa;
  dipdgapd(17,29)=-sa;
  dipdgapd(18,30)=-sa;
  dipdgapd(21,33)=-sa;
  dipdgapd(22,34)=-sa;
  
  dipdgaps(13,25)=-da;
  dipdgaps(14,26)=-da;
  dipdgaps(17,29)=-da;
  dipdgaps(18,30)=-da;
  dipdgaps(21,33)=-da;
  dipdgaps(22,34)=-da;
  
  return 0;
}
