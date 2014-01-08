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

#include <o2scl/cfl_njl_eos.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

cfl_njl_eos::cfl_njl_eos() {

  // Allocate storage for the eigenproblem
  
  iprop=gsl_matrix_complex_alloc(12,12);
  eivec=gsl_matrix_complex_alloc(12,12);
  eval=gsl_vector_alloc(12);

  dipdgapu.resize(12,12);
  dipdgapd.resize(12,12);
  dipdgaps.resize(12,12);
  w=gsl_eigen_hermv_alloc(12);
  
  // Other initializations
  quartic=&def_quartic;
  
  gap_limit=1.0e-6;
  zerot=false;
  fixed_mass=false;
  fromqq=true;
  eq_limit=1.0e-6;

  smu3=0.0;
  smu8=0.0;

  integ_test=false;
  inte_epsabs=1.0e-4;
  inte_epsrel=1.0e-4;
  inte_npoints=0;
}

cfl_njl_eos::~cfl_njl_eos() {
  gsl_matrix_complex_free(iprop);
  gsl_matrix_complex_free(eivec);
  gsl_vector_free(eval);
  gsl_eigen_hermv_free(w);
}

double cfl_njl_eos::rescale_error(double err, double result_abs, 
				  double result_asc) {
  err = fabs(err);
  
  if (result_asc != 0 && err != 0) {
    
    double scale = pow((200 * err / result_asc), 1.5);
    
    if (scale < 1) {
      err = result_asc * scale;
    } else {
      err = result_asc;
    }
  }
  
#ifndef O2SCL_NO_CPP11
      double dbl_eps=std::numeric_limits<double>::epsilon();
#else 
      double dbl_eps=GSL_DBL_EPSILON;
#endif

  if (result_abs > GSL_DBL_MIN / (50 * dbl_eps)) {
    
    double min_err = 50 * dbl_eps * result_abs;
    
    if (min_err > err) {
      err = min_err;
    }
  }
  
  return err;
}

int cfl_njl_eos::test_integration(test_mgr &t) {
  integ_test=true;
  ubvector res(3);
  double err;
  integ_err(0.0,pi,3,res,err);
  t.test_rel(res[0],2.0,1.0e-12,"ti1");
  t.test_rel(res[1],2.0,1.0e-12,"ti2");
  t.test_rel(res[2]+1.0,1.0,1.0e-12,"ti3");
  integ_test=false;
  return 0;
}

int cfl_njl_eos::integ_err(double a, double b, const size_t nr,
			    ubvector &res, double &err2) {

  //double fval1[nr], fval2[nr], fval[nr];
  double fv1[5], fv2[5], fv3[5], fv4[5];
  //double savfun[21][nr];  
  //double res10[nr], res21[nr], res43[nr], res87[nr];    
  //double result_kronrod[nr];
  double err;
  double resabs; 
  double resasc; 

  double *fval1=new double[nr];
  double *fval2=new double[nr];
  double *fval=new double[nr];
  double **savfun=new double *[21];
  for(size_t i=0;i<21;i++) savfun[i]=new double[nr];
  double *res10=new double[nr];
  double *res21=new double[nr];
  double *res43=new double[nr];
  double *res87=new double[nr];
  double *result_kronrod=new double[nr];
  double *f_center=new double[nr];

  const double half_length= 0.5*(b-a);
  const double abs_half_length=fabs (half_length);
  const double center=0.5*(b+a);

  //double f_center[nr];
  integrands(center,f_center);
    
  int k;
  
#ifndef O2SCL_NO_CPP11
      double dbl_eps=std::numeric_limits<double>::epsilon();
#else 
      double dbl_eps=GSL_DBL_EPSILON;
#endif

  if (inte_epsabs<=0 && (inte_epsrel<50*dbl_eps || 
			 inte_epsrel<0.5e-28)) {
    for(size_t j=0;j<res.size();j++) res[j]=0.0;
    err2=0;
    O2SCL_ERR2_RET("Tolerance cannot be acheived with given epsabs and epsrel",
		   " in cfl_njl_eos::integ_err().",exc_ebadtol);
  };
  
  for (size_t j=0;j<nr;j++) {
    res10[j]=0;
    res21[j]=o2scl_inte_qng_coeffs::w21b[5]*f_center[j];
  }
  resabs=o2scl_inte_qng_coeffs::w21b[5]*fabs(f_center[0]);
    
  for (k=0;k<5;k++) {
    const double abscissa=half_length*o2scl_inte_qng_coeffs::x1[k];
    integrands(center+abscissa,fval1);
    integrands(center-abscissa,fval2);
    for(size_t j=0;j<nr;j++) {
      fval[j]=fval1[j]+fval2[j];
      res10[j]+=o2scl_inte_qng_coeffs::w10[k]*fval[j];
      res21[j]+=o2scl_inte_qng_coeffs::w21a[k]*fval[j];
      savfun[k][j]=fval[j];
    }
    resabs+=o2scl_inte_qng_coeffs::w21a[k]*(fabs(fval1[0])+fabs(fval2[0]));
    fv1[k]=fval1[0];
    fv2[k]=fval2[0];
  }
    
  for (k=0;k<5;k++) {
    const double abscissa=half_length*o2scl_inte_qng_coeffs::x2[k];
    integrands(center+abscissa,fval1);
    integrands(center-abscissa,fval2);
    for(size_t j=0;j<nr;j++) {
      fval[j]=fval1[j]+fval2[j];
      res21[j]+=o2scl_inte_qng_coeffs::w21b[k]*fval[j];
      savfun[k+5][j]=fval[j];
    }
    resabs+=o2scl_inte_qng_coeffs::w21b[k]*(fabs(fval1[0])+fabs(fval2[0]));
    fv3[k]=fval1[0];
    fv4[k]=fval2[0];
  }
  
  resabs*=abs_half_length;
  const double mean=0.5*res21[0];

  resasc=o2scl_inte_qng_coeffs::w21b[5]*fabs(f_center[0]-mean);
  for (k=0;k<5;k++) {
    resasc+=(o2scl_inte_qng_coeffs::w21a[k]*(fabs(fv1[k]-mean)+
					   fabs(fv2[k]-mean))
	     +o2scl_inte_qng_coeffs::w21b[k]*(fabs(fv3[k]-mean)+
					    fabs(fv4[k]-mean)));
  }
  resasc*=abs_half_length;
  for(size_t j=0;j<nr;j++) result_kronrod[j]=res21[j]*half_length;
  err=rescale_error((res21[0]-res10[0])*half_length,resabs,resasc);
    
  if (err < inte_epsabs || err < inte_epsrel*fabs (result_kronrod[0])) {
    for(size_t j=0;j<nr;j++) res[j]=result_kronrod[j];
    err2=err;
    inte_npoints=21;
  { 
    for(size_t i=0;i<21;i++) delete[] savfun[i];
    delete[] fval1;
    delete[] fval2;
    delete[] fval;
    delete[] savfun;
    delete[] res10;
    delete[] res21;
    delete[] res43;
    delete[] res87;
    delete[] result_kronrod;
    delete[] f_center;
  }
    
    return success;
  }
      
  for(size_t j=0;j<nr;j++) {
    res43[j]=o2scl_inte_qng_coeffs::w43b[11]*f_center[j];
    for (k=0;k<10;k++) {
      res43[j]+=savfun[k][j]*o2scl_inte_qng_coeffs::w43a[k];
    }
  }
      
  for (k=0; k < 11; k++) {
    const double abscissa=half_length*o2scl_inte_qng_coeffs::x3[k];
    integrands(center+abscissa,fval1);
    integrands(center-abscissa,fval2);
    for(size_t j=0;j<nr;j++) {
      fval[j]=fval1[j]+fval2[j];
      res43[j]+=fval[j]*o2scl_inte_qng_coeffs::w43b[k];
      savfun[k+10][j]=fval[j];
    }
  }
     
  for(size_t j=0;j<nr;j++) {
    result_kronrod[j]=res43[j]*half_length;
  }
  err=rescale_error((res43[0]-res21[0])*half_length,resabs,resasc);
  
  if (err < inte_epsabs || err < inte_epsrel*fabs (result_kronrod[0])) {
    for(size_t j=0;j<nr;j++) res[j]=result_kronrod[j];
    err2=err;
    inte_npoints=43;
  { 
    for(size_t i=0;i<21;i++) delete[] savfun[i];
    delete[] fval1;
    delete[] fval2;
    delete[] fval;
    delete[] savfun;
    delete[] res10;
    delete[] res21;
    delete[] res43;
    delete[] res87;
    delete[] result_kronrod;
    delete[] f_center;
  }
    
    return success;
  }
      
  for(size_t j=0;j<nr;j++) {
    res87[j]=o2scl_inte_qng_coeffs::w87b[22]*f_center[j];
    for (k=0;k<21;k++) {
      res87[j]+=savfun[k][j]*o2scl_inte_qng_coeffs::w87a[k];
    }
  }
  
  for (k=0;k<22;k++) {
    const double abscissa=half_length*o2scl_inte_qng_coeffs::x4[k];
    integrands(center+abscissa,fval1);
    integrands(center-abscissa,fval2);
    for(size_t j=0;j<nr;j++) {
      res87[j]+=o2scl_inte_qng_coeffs::w87b[k]*(fval1[j]+fval2[j]);
    }
  }
  
  for(size_t j=0;j<nr;j++) result_kronrod[j]=res87[j]*half_length;
  err=rescale_error((res87[0]-res43[0])*half_length,resabs,resasc);
  
  for(size_t j=0;j<nr;j++) res[j]=result_kronrod[j];
  err2=err;

  if (err < inte_epsabs || err < inte_epsrel*fabs (result_kronrod[0])) {
    inte_npoints=87;
  { 
    for(size_t i=0;i<21;i++) delete[] savfun[i];
    delete[] fval1;
    delete[] fval2;
    delete[] fval;
    delete[] savfun;
    delete[] res10;
    delete[] res21;
    delete[] res43;
    delete[] res87;
    delete[] result_kronrod;
    delete[] f_center;
  }
    
    return success;
  }
      
  { 
    for(size_t i=0;i<21;i++) delete[] savfun[i];
    delete[] fval1;
    delete[] fval2;
    delete[] fval;
    delete[] savfun;
    delete[] res10;
    delete[] res21;
    delete[] res43;
    delete[] res87;
    delete[] result_kronrod;
    delete[] f_center;
  }
  inte_npoints=88;
  O2SCL_ERR_RET("failed to reach tolerance with highest-order rule",
		exc_etol);
}
  
int cfl_njl_eos::set_parameters(double lambda, double fourferm, 
				 double sixferm, double fourgap) {
  int ret;

  if (fixed_mass) {

    L=lambda;
    G=0.0;
    K=0.0;
    GD=fourgap;
    
  } else {

    // Run the parent method first
    ret=nambujl_eos::set_parameters(lambda,fourferm,sixferm);

    if (fourgap==0.0) {
      GD=3.0*G/4.0;
    } else {
      GD=fourgap;
    }
  }

  return 0;
}

int cfl_njl_eos::calc_eq_temp_p(quark &u, quark &d, quark &s,
				double &qq1, double &qq2, double &qq3, 
				double &gap1, double &gap2, double &gap3, 
				double mu3, double mu8, 
				double &n3, double &n8, 
				thermo &qb, const double ltemper) {

  if (fromqq==false) {
    O2SCL_ERR_RET("fromqq=false in cfl_njl_eos::calc_eq_temp_p()",exc_efailed);
  }
  
  up=&u;
  down=&d;
  strange=&s;

  // Store the temperature and color chemical potentials
  temper=ltemper;
  smu3=mu3;
  smu8=mu8;

  // In the case that the gaps are small enough, run the ungapped version
  
  if ((u.del)<gap_limit && (d.del)<gap_limit && (s.del)<gap_limit) {
    int ret;
    ret=nambujl_eos::calc_eq_temp_p(u,d,s,qq1,qq2,qq3,qb,temper);
    gap1=u.del;
    gap2=d.del;
    gap3=s.del;
    return ret;
  }

  // Compute the masses

  if (fixed_mass) {
    u.ms=u.m;
    d.ms=d.m;
    s.ms=s.m;
  } else {
    u.ms=u.m-4.0*G*u.qq+2.0*K*d.qq*s.qq;
    d.ms=d.m-4.0*G*d.qq+2.0*K*u.qq*s.qq;
    s.ms=s.m-4.0*G*s.qq+2.0*K*d.qq*u.qq;
  }
  double dmudqqu=-4.0*G;
  double dmudqqd=2.0*K*s.qq;
  double dmudqqs=2.0*K*d.qq;
  double dmddqqu=2.0*K*s.qq;
  double dmddqqd=-4.0*G;
  double dmddqqs=2.0*K*u.qq;
  double dmsdqqu=2.0*K*d.qq;
  double dmsdqqd=2.0*K*u.qq;
  double dmsdqqs=-4.0*G;
  
  // Compute the integrals

  ubvector res(13);
  double err2;

  integ_err(0.0,L,13,res,err2);
  qb.pr=-res[0];
  qb.en=res[1];
  u.n=-res[2];
  d.n=-res[3];
  s.n=-res[4];
  qq1=res[5]*dmudqqu+res[6]*dmddqqu+res[7]*dmsdqqu;
  qq2=res[5]*dmudqqd+res[6]*dmddqqd+res[7]*dmsdqqd;
  qq3=res[5]*dmudqqs+res[6]*dmddqqs+res[7]*dmsdqqs;
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
  }

  qq1+=4.0*G*u.qq-4.0*K*d.qq*s.qq;
  qq2+=4.0*G*d.qq-4.0*K*u.qq*s.qq;
  qq3+=4.0*G*s.qq-4.0*K*u.qq*d.qq;

  qb.pr-=u.del*u.del/4.0/GD;
  qb.pr-=d.del*d.del/4.0/GD;
  qb.pr-=s.del*s.del/4.0/GD;

  gap1+=u.del/2.0/GD;
  gap2+=d.del/2.0/GD;
  gap3+=s.del/2.0/GD;

  // Compute the energy densities 
  
  u.ed=-u.pr+u.n*u.mu;
  d.ed=-d.pr+d.n*d.mu;
  s.ed=-s.pr+s.n*s.mu;
  qb.ed=qb.en*temper-qb.pr+u.mu*u.n+d.mu*d.n+s.mu*s.n;

  return 0;

}

int cfl_njl_eos::integrands(double p, double res[]) {

  if (integ_test) {
    res[0]=sin(p);
    res[1]=sin(p);
    res[2]=cos(p);
    return 0;
  }

  int k;
  double egv[36];
  double dedmuu[36], dedmud[36], dedmus[36];
  double dedmu[36], dedmd[36], dedms[36], deds[36];
  double dedd[36], dedu[36], dedmu3[36], dedmu8[36];
  
  eigenvalues(p,smu3,smu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,dedms,
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
	res[5]+=dedmu[k]/2.0;
	res[6]+=dedmd[k]/2.0;
	res[7]+=dedms[k]/2.0;
	res[8]+=dedu[k]/2.0;
	res[9]+=dedd[k]/2.0;
	res[10]+=deds[k]/2.0;
	res[11]+=dedmu3[k]/2.0;
	res[12]+=dedmu8[k]/2.0;
      } else {
	res[2]-=dedmuu[k]/2.0;
	res[3]-=dedmud[k]/2.0;
	res[4]-=dedmus[k]/2.0;
	res[5]-=dedmu[k]/2.0;
	res[6]-=dedmd[k]/2.0;
	res[7]-=dedms[k]/2.0;
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
      res[5]+=dedmu[k]/2.0;
      res[6]+=dedmd[k]/2.0;
      res[7]+=dedms[k]/2.0;
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
	res[5]-=dedmu[k];
	res[6]-=dedmd[k];
	res[7]-=dedms[k];
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
	res[5]-=1.0/(1.0+exp(egv[k]/temper))*dedmu[k];
	res[6]-=1.0/(1.0+exp(egv[k]/temper))*dedmd[k];
	res[7]-=1.0/(1.0+exp(egv[k]/temper))*dedms[k];
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

int cfl_njl_eos::test_normal_eigenvalues(test_mgr &t) {

  double h=1.0e-7;
  double lam[2], dldmu[2], dldm[2];
  double mu1[2], m1[2], mu2[2], m2[2];
  normal_eigenvalues(0.5,1.0,2.0,lam,dldmu,dldm);
  m1[0]=lam[0];
  m1[1]=lam[1];
  mu1[0]=lam[0];
  mu1[1]=lam[1];
  normal_eigenvalues(0.5+h,1.0,2.0,lam,dldmu,dldm);
  m2[0]=lam[0];
  m2[1]=lam[1];
  normal_eigenvalues(0.5,1.0,2.0+h,lam,dldmu,dldm);
  mu2[0]=lam[0];
  mu2[1]=lam[1];
  normal_eigenvalues(0.5,1.0,2.0,lam,dldmu,dldm);

  t.test_rel((m2[0]-m1[0])/h,dldm[0],1.0e-5,"");
  t.test_rel((m2[1]-m1[1])/h,dldm[1],1.0e-5,"");
  t.test_rel((mu2[0]-mu1[0])/h,dldmu[0],1.0e-5,"");
  t.test_rel((mu2[1]-mu1[1])/h,dldmu[1],1.0e-5,"");

  return 0;
}

int cfl_njl_eos::normal_eigenvalues(double ms, double lmom, double mu, 
				     double lam[2], double dldmu[2], 
				     double dldm[2]) {
  double E;
  
  // ungapped limit

  if (ms==0.0) {
    E=lmom;
  } else {
    E=sqrt(ms*ms+lmom*lmom);
  }
  
  lam[0]=-mu-E;
  lam[1]=-mu+E;
  dldmu[0]=-1.0;
  dldmu[1]=-1.0;
  dldm[0]=-ms/E;
  dldm[1]=ms/E;

  return 0;
}

int cfl_njl_eos::test_gapped_eigenvalues(test_mgr &t) {

  double h=1.0e-7;
  double lam[4], dldmu1[4], dldm1[4], dldmu2[4], dldm2[4], dldg[4];
  double start[4], mu1b[4], mu2b[4], m1b[4], m2b[4], gb[4], lamb[4];
  
  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0,0.0,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0,0.01,
		     lamb,mu1b,mu2b,m1b,m2b,gb);
  for(size_t i=0;i<4;i++) {
    t.test_rel(lam[i],lamb[i],1.0e-4,"lam");
  }
  
  gapped_eigenvalues(0.5,0.5001,1.0,2.0,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  gapped_eigenvalues(0.5,0.5,1.0,2.0,3.0,0.6,
		     lamb,mu1b,mu2b,m1b,m2b,gb);

  for(size_t i=0;i<4;i++) {
    t.test_rel(lam[i],lamb[i],1.0e-3,"lam");
    t.test_rel(dldmu1[i],mu1b[i],1.0e-3,"mu1");
    t.test_rel(dldmu2[i],mu2b[i],1.0e-3,"mu2");
    t.test_rel(dldg[i],gb[i],1.0e-3,"g");
  }

  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,start);

  gapped_eigenvalues(0.5+h,0.7,1.0,2.0,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,m1b);

  gapped_eigenvalues(0.5,0.7+h,1.0,2.0,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,m2b);
  
  gapped_eigenvalues(0.5,0.7,1.0,2.0+h,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,mu1b);

  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0+h,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,mu2b);

  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0,0.6+h,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,gb);

  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  
  for(size_t i=0;i<4;i++) {
    t.test_rel((m1b[i]-start[i])/h,dldm1[i],1.0e-5,"");
    t.test_rel((m2b[i]-start[i])/h,dldm2[i],1.0e-5,"");
    t.test_rel((mu1b[i]-start[i])/h,dldmu1[i],1.0e-5,"");
    t.test_rel((mu2b[i]-start[i])/h,dldmu2[i],1.0e-5,"");
    t.test_rel((gb[i]-start[i])/h,dldg[i],1.0e-5,"");
  }

  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0,0.0,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,start);

  gapped_eigenvalues(0.5+h,0.7,1.0,2.0,3.0,0.0,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,m1b);

  gapped_eigenvalues(0.5,0.7+h,1.0,2.0,3.0,0.0,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,m2b);
  
  gapped_eigenvalues(0.5,0.7,1.0,2.0+h,3.0,0.0,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,mu1b);

  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0+h,0.0,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,mu2b);

  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0,0.0+h,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,gb);

  gapped_eigenvalues(0.5,0.7,1.0,2.0,3.0,0.0,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  
  for(size_t i=0;i<4;i++) {
    t.test_rel((m1b[i]-start[i])/h,dldm1[i],1.0e-5,"");
    t.test_rel((m2b[i]-start[i])/h,dldm2[i],1.0e-5,"");
    t.test_rel((mu1b[i]-start[i])/h,dldmu1[i],1.0e-5,"");
    t.test_rel((mu2b[i]-start[i])/h,dldmu2[i],1.0e-5,"");
    t.test_rel((gb[i]-start[i])/h,dldg[i],1.0e-5,"");
  }

  double tmp=eq_limit;
  eq_limit=1.0;
  gapped_eigenvalues(0.5,0.5,1.0,2.0,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,start);

  gapped_eigenvalues(0.5+h,0.5,1.0,2.0,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,m1b);

  gapped_eigenvalues(0.5,0.5,1.0,2.0+h,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,mu1b);

  gapped_eigenvalues(0.5,0.5,1.0,2.0,3.0+h,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,mu2b);

  gapped_eigenvalues(0.5,0.5,1.0,2.0,3.0,0.6+h,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  vector_copy(4,lam,gb);

  gapped_eigenvalues(0.5,0.5,1.0,2.0,3.0,0.6,
		     lam,dldmu1,dldmu2,dldm1,dldm2,dldg);
  eq_limit=tmp;
  
  for(size_t i=0;i<4;i++) {
    t.test_rel((mu1b[i]-start[i])/h,dldmu1[i],1.0e-5,"mu1");
    t.test_rel((mu2b[i]-start[i])/h,dldmu2[i],1.0e-5,"mu2");
    t.test_rel((gb[i]-start[i])/h,dldg[i],1.0e-3,"g");
  }

  return 0;
}


/// Treat the simply gapped quarks in all cases gracefully
int cfl_njl_eos::gapped_eigenvalues(double ms1, double ms2, double lmom,
				     double mu1, double mu2, double tdelta,
				     double lam[4], double dldmu1[4], 
				     double dldmu2[4], double dldm1[4],
				     double dldm2[4], double dldg[4]) {
  double E1, E2;

  if (tdelta<gap_limit) {
    
    // ungapped limit

    if (ms1==0.0 && ms2==0.0) {
      E1=lmom;
      E2=lmom;
    } else {
      E1=sqrt(ms1*ms1+lmom*lmom);
      E2=sqrt(ms2*ms2+lmom*lmom);
    }
    
    lam[0]=mu2+E2;
    lam[1]=-mu1-E1;
    lam[2]=-mu1+E1;
    lam[3]=mu2-E2;

    dldmu1[0]=0.0;
    dldmu1[1]=-1.0;
    dldmu1[2]=-1.0;
    dldmu1[3]=0.0;

    dldmu2[0]=1.0;
    dldmu2[1]=0.0;
    dldmu2[2]=0.0;
    dldmu2[3]=1.0;

    dldm1[0]=0.0;
    dldm1[1]=-ms1/E1;
    dldm1[2]=ms1/E1;
    dldm1[3]=0.0;

    dldm2[0]=ms2/E2;
    dldm2[1]=0.0;
    dldm2[2]=0.0;
    dldm2[3]=-ms2/E2;

    dldg[2]=0.0;
    dldg[0]=0.0;
    dldg[1]=0.0;
    dldg[3]=0.0;
    
  } else if (false && fabs(ms1-ms2)<eq_limit) {
    
    // equal mass limit

    if (ms1==0.0 && ms2==0.0) {
      E1=lmom;
    } else {
      E1=sqrt(ms1*ms1+lmom*lmom);
    }
    
    double smu2=(mu1+mu2)/2.0;
    double dmu2=(mu2-mu1)/2.0;
    double dEdm=ms1/E1;
    double dm=sqrt(pow(smu2-E1,2.0)+tdelta*tdelta);
    double dp=sqrt(pow(smu2+E1,2.0)+tdelta*tdelta);
    lam[0]=dmu2+dp;
    lam[1]=dmu2-dp;
    lam[2]=dmu2-dm;
    lam[3]=dmu2+dm;
    dldmu1[0]=-0.5+(smu2+E1)/dp/2.0;
    dldmu1[1]=-0.5-(smu2+E1)/dp/2.0;
    dldmu1[2]=-0.5-(smu2-E1)/dm/2.0;
    dldmu1[3]=-0.5+(smu2-E1)/dm/2.0;
    dldmu2[0]=0.5+(smu2+E1)/dp/2.0;
    dldmu2[1]=0.5-(smu2+E1)/dp/2.0;
    dldmu2[2]=0.5-(smu2-E1)/dm/2.0;
    dldmu2[3]=0.5+(smu2-E1)/dm/2.0;
    dldm1[0]=(smu2+E1)/dp*dEdm;
    dldm1[1]=-(smu2+E1)/dp*dEdm;
    dldm1[2]=(smu2-E1)/dm*dEdm;
    dldm1[3]=-(smu2-E1)/dm*dEdm;
    dldm2[0]=(smu2+E1)/dp*dEdm;
    dldm2[1]=-(smu2+E1)/dp*dEdm;
    dldm2[2]=(smu2-E1)/dm*dEdm;
    dldm2[3]=-(smu2-E1)/dm*dEdm;
    dldg[0]=tdelta/dp;
    dldg[1]=-tdelta/dp;
    dldg[2]=-tdelta/dm;
    dldg[3]=tdelta/dm;
    
  } else {

    // generic case
    
    double coef1, coef2, coef3, coef4, coef5;

    coef1=1.0;
    coef2=2.0*(mu1-mu2);
    coef3=mu1*mu1-4.0*mu1*mu2+mu2*mu2-2.0*lmom*lmom-2.0*tdelta*tdelta-
      ms1*ms1-ms2*ms2;
    coef4=(mu1-mu2)*(-2.0*tdelta*tdelta-2.0*mu1*mu2-
		     2.0*lmom*lmom)+2.0*ms1*ms1*mu2-2.0*ms2*ms2*mu1;
    coef5=mu1*mu1*mu2*mu2+lmom*lmom*(ms1*ms1+ms2*ms2+2.0*tdelta*tdelta-
				     mu1*mu1-mu2*mu2)-mu1*mu1*ms2*ms2-
      mu2*mu2*ms1*ms1+pow(tdelta,4.0)+pow(lmom,4.0)+
      2.0*ms1*ms2*tdelta*tdelta+ms1*ms1*ms2*ms2+2.0*mu1*mu2*tdelta*tdelta;

    double d2dmu1=2.0;
    double d2dmu2=-2.0;
    double d3dmu1=2.0*mu1-4.0*mu2;
    double d3dmu2=2.0*mu2-4.0*mu1;
    double d3dm1=-2.0*ms1;
    double d3dm2=-2.0*ms2;
    double d3dg=-4.0*tdelta;
    double d4dmu1=-2.0*ms2*ms2-2.0*mu2*(mu1-mu2)+
      (-2.0*lmom*lmom-2.0*tdelta*tdelta-2.0*mu1*mu2);
    double d4dmu2=2.0*ms1*ms1-2.0*mu1*(mu1-mu2)-
      (-2.0*lmom*lmom-2.0*tdelta*tdelta-2.0*mu1*mu2);
    double d4dm1=4.0*ms1*mu2;
    double d4dm2=-4.0*ms2*mu1;
    double d4dg=-4.0*(mu1-mu2)*tdelta;
    double d5dmu1=-2.0*ms2*ms2*mu1+2.0*tdelta*tdelta*mu2+2.0*mu1*mu2*mu2-
      2.0*mu1*lmom*lmom;
    double d5dmu2=-2.0*ms1*ms1*mu2+2.0*tdelta*tdelta*mu1+2.0*mu1*mu1*mu2-
      2.0*mu2*lmom*lmom;
    double d5dm1=2.0*ms1*ms2*ms2+2.0*ms2*tdelta*tdelta-2.0*ms1*mu2*mu2+
      2.0*lmom*lmom*ms1;
    double d5dm2=2.0*ms2*ms1*ms1+2.0*ms1*tdelta*tdelta-2.0*ms2*mu1*mu1+
      2.0*lmom*lmom*ms2;
    double d5dg=4.0*tdelta*tdelta*tdelta+4.0*tdelta*lmom*lmom+
      4.0*ms1*ms2*tdelta+4.0*mu1*mu2*tdelta;

    std::complex<double> a1, a2, a3, a4;
    quartic->solve_rc(coef1,coef2,coef3,coef4,coef5,a1,a2,a3,a4);
    lam[0]=a1.real();
    lam[1]=a2.real();
    lam[2]=a3.real();
    lam[3]=a4.real();

    double lam1, lam2, lam3;
    for(size_t i=0;i<4;i++) {
      lam1=lam[i];
      lam2=lam1*lam1;
      lam3=lam1*lam2;
      double den=4.0*coef1*lam3+3.0*coef2*lam2+2.0*coef3*lam1+coef4;
      dldmu1[i]=(-lam3*d2dmu1-lam2*d3dmu1-lam1*d4dmu1-d5dmu1)/den;
      dldmu2[i]=(-lam3*d2dmu2-lam2*d3dmu2-lam1*d4dmu2-d5dmu2)/den;
      dldm1[i]=(-lam2*d3dm1-lam1*d4dm1-d5dm1)/den;
      dldm2[i]=(-lam2*d3dm2-lam1*d4dm2-d5dm2)/den;
      dldg[i]=(-lam2*d3dg-lam1*d4dg-d5dg)/den;
    }
  }
  
  return 0;
}

int cfl_njl_eos::test_derivatives(double lmom, double mu3, double mu8,
				    test_mgr &t) {
  double egv[36];
  double dedmuu[36], dedmud[36], dedmus[36];
  double dedmu[36], dedmd[36], dedms[36], deds[36];
  double dedd[36], dedu[36], dedmu3[36], dedmu8[36];
  double h=1.0e-6;
  double d1[36], d2[36];
  // Check muu
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmuu,d2);
    up->mu+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"muu");
    }
    up->mu-=h;
  }
  // Check mud
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmud,d2);
    down->mu+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"mud");
    }
    down->mu-=h;
  }
  // Check mus
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmus,d2);
    strange->mu+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"mus");
    }
    strange->mu-=h;
  }
  // Check mu
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmu,d2);
    up->ms+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"mu");
    }
    up->ms-=h;
  }
  // Check md
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmd,d2);
    down->ms+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"md");
    }
    down->ms-=h;
  }
  // Check ms
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedms,d2);
    strange->ms+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"ms");
    }
    strange->ms-=h;
  }
  // Check delu
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedu,d2);
    up->del+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"delu");
    }
    up->del-=h;
  }
  // Check deld
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedd,d2);
    down->del+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"deld");
    }
    down->del-=h;
  }
  // Check dels
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,deds,d2);
    strange->del+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"dels");
    }
    strange->del-=h;
  }
  // Check mu3
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmu3,d2);
    mu3+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"mu3");
    }
    mu3-=h;
  }
  // Check mu8
  {
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    vector_copy(36,egv,d1);
    vector_copy(36,dedmu8,d2);
    mu8+=h;
    eigenvalues(lmom,mu3,mu8,egv,dedmuu,dedmud,dedmus,dedmu,dedmd,
		dedms,dedu,dedd,deds,dedmu3,dedmu8);
    for(size_t i=0;i<36;i++) {
      t.test_rel(d2[i],(egv[i]-d1[i])/h,1.0e-4,"mu8");
    }
    mu8-=h;
  }
  return 0;
}


int cfl_njl_eos::eigenvalues(double lmom, double mu3, 
			      double mu8, double egv[36],
			      double dedmuu[36], double dedmud[36],
			      double dedmus[36], double dedmu[36], 
			      double dedmd[36], double dedms[36],
			      double dedu[36], double dedd[36],
			      double deds[36], double dedmu3[36],
			      double dedmu8[36]) {

  gsl_vector *gvp;
  int k;
  double coef1, coef2, coef3, coef4, coef5, mur, mug, mub;
  quark *u=up, *d=down, *s=strange;
  double s3=sqrt(3.0);
  
  mur=mu3+mu8/s3;
  mug=-mu3+mu8/s3;
  mub=-2.0*mu8/s3;
  double drd3=1.0, drd8=1.0/s3, dgd3=-1.0, dgd8=1.0/s3;
  double dbd3=0.0, dbd8=-2.0/s3;
  
  for(k=0;k<36;k++) {
    egv[k]=0.0;
    dedmuu[k]=0.0;
    dedmud[k]=0.0;
    dedmus[k]=0.0;
    dedmu[k]=0.0;
    dedmd[k]=0.0;
    dedms[k]=0.0;
    dedu[k]=0.0;
    dedd[k]=0.0;
    deds[k]=0.0;
    dedmu3[k]=0.0;
    dedmu8[k]=0.0;
  }

  bool add_neg=true;
  if (u->del<gap_limit && d->del<gap_limit) {
    // The standard 2sc phase:
    
    gapped_eigenvalues(u->ms,d->ms,lmom,u->mu+mur,d->mu+mug,s->del,
		       &egv[0],&dedmuu[0],&dedmud[0],&dedmu[0],&dedmd[0],
		       &deds[0]);
    normal_eigenvalues(s->ms,lmom,s->mu+mub,&egv[4],&dedmus[4],&dedms[4]);
    
  } else if (s->del<gap_limit && d->del<gap_limit) {
    // The 2SCds phase:
    
    gapped_eigenvalues(d->ms,s->ms,lmom,d->mu+mug,s->mu+mub,u->del,
		       &egv[0],&dedmud[0],&dedmus[0],&dedmd[0],&dedms[0],
		       &dedu[0]);
    normal_eigenvalues(u->ms,lmom,u->mu+mur,&egv[4],&dedmuu[4],&dedmu[4]);
    
  } else if (u->del<gap_limit && s->del<gap_limit) {
    // The 2SCus phase:

    gapped_eigenvalues(u->ms,s->ms,lmom,u->mu+mur,s->mu+mub,u->del,
		       &egv[0],&dedmuu[0],&dedmus[0],&dedmu[0],&dedms[0],
		       &dedd[0]);
    normal_eigenvalues(d->ms,lmom,d->mu+mug,&egv[4],&dedmud[4],&dedmd[4]);

  } else {
    
    gsl_complex zero={{0.0,0.0}};
    std::complex<double> zero2(0.0,0.0);
    std::complex<double> littlei2(0.0,1.0);

    for(size_t i=0;i<12;i++) {
      for(size_t j=0;j<12;j++) {
	gsl_matrix_complex_set(iprop,i,j,zero);
	dipdgapu(i,j)=zero2;
	dipdgapd(i,j)=zero2;
	dipdgaps(i,j)=zero2;
      }
    }
    
    ubvector mmuu(12), mmud(12), mmus(12), mmu(12), mmd(12), mms(12);
    for(k=0;k<12;k++) {
      mmuu[k]=0.0;
      mmud[k]=0.0;
      mmus[k]=0.0;
      mmu[k]=0.0;
      mmd[k]=0.0;
      mms[k]=0.0;
    }
    
    // Since the matrix is destroyed, we have to fill the entries 
    // every time. Some of the entries are commented out, since 
    // they represent the part of the matrix that isn't needed.

    gsl_complex tmp={{0.0,0.0}};

    // Set diagonal entries
    tmp=gsl_matrix_complex_get(iprop,0,0);
    GSL_REAL(tmp)=(u->mu+mur)-u->ms;
    gsl_matrix_complex_set(iprop,0,0,tmp);

    tmp=gsl_matrix_complex_get(iprop,10,10);
    GSL_REAL(tmp)=-(u->mu+mur)-u->ms;
    gsl_matrix_complex_set(iprop,10,10,tmp);

    tmp=gsl_matrix_complex_get(iprop,1,1);
    GSL_REAL(tmp)=(u->mu+mur)+u->ms;
    gsl_matrix_complex_set(iprop,1,1,tmp);
    
    tmp=gsl_matrix_complex_get(iprop,11,11);
    GSL_REAL(tmp)=-(u->mu+mur)+u->ms;
    gsl_matrix_complex_set(iprop,11,11,tmp);
    
    mmuu[0]=1.0;
    mmuu[1]=1.0;
    mmuu[10]=-1.0;
    mmuu[11]=-1.0;
    mmu[0]=-1.0;
    mmu[1]=1.0;
    mmu[10]=-1.0;
    mmu[11]=1.0;
    
    tmp=gsl_matrix_complex_get(iprop,6,6);
    GSL_REAL(tmp)=(d->mu+mug)-d->ms;
    gsl_matrix_complex_set(iprop,6,6,tmp);

    tmp=gsl_matrix_complex_get(iprop,4,4);
    GSL_REAL(tmp)=-(d->mu+mug)-d->ms;
    gsl_matrix_complex_set(iprop,4,4,tmp);

    tmp=gsl_matrix_complex_get(iprop,2,2);
    GSL_REAL(tmp)=-(d->mu+mug)+d->ms;
    gsl_matrix_complex_set(iprop,2,2,tmp);

    tmp=gsl_matrix_complex_get(iprop,7,7);
    GSL_REAL(tmp)=(d->mu+mug)+d->ms;
    gsl_matrix_complex_set(iprop,7,7,tmp);

    mmud[2]=-1.0;
    mmud[4]=-1.0;
    mmud[6]=1.0;
    mmud[7]=1.0;
    mmd[2]=1.0;
    mmd[4]=-1.0;
    mmd[6]=-1.0;
    mmd[7]=1.0;
    
    tmp=gsl_matrix_complex_get(iprop,5,5);
    GSL_REAL(tmp)=-(s->mu+mub)-s->ms;
    gsl_matrix_complex_set(iprop,5,5,tmp);

    tmp=gsl_matrix_complex_get(iprop,8,8);
    GSL_REAL(tmp)=(s->mu+mub)-s->ms;
    gsl_matrix_complex_set(iprop,8,8,tmp);

    tmp=gsl_matrix_complex_get(iprop,9,9);
    GSL_REAL(tmp)=(s->mu+mub)+s->ms;
    gsl_matrix_complex_set(iprop,9,9,tmp);

    tmp=gsl_matrix_complex_get(iprop,3,3);
    GSL_REAL(tmp)=-(s->mu+mub)+s->ms;
    gsl_matrix_complex_set(iprop,3,3,tmp);

    mmus[3]=-1.0;
    mmus[5]=-1.0;
    mmus[8]=1.0;
    mmus[9]=1.0;
    mms[3]=1.0;
    mms[5]=-1.0;
    mms[8]=-1.0;
    mms[9]=1.0;
    
    // Set momentum entries
    tmp=gsl_matrix_complex_get(iprop,1,0);
    GSL_REAL(tmp)=lmom;
    gsl_matrix_complex_set(iprop,1,0,tmp);

    tmp=gsl_matrix_complex_get(iprop,4,2);
    GSL_REAL(tmp)=lmom;
    gsl_matrix_complex_set(iprop,4,2,tmp);

    tmp=gsl_matrix_complex_get(iprop,5,3);
    GSL_REAL(tmp)=lmom;
    gsl_matrix_complex_set(iprop,5,3,tmp);

    tmp=gsl_matrix_complex_get(iprop,7,6);
    GSL_REAL(tmp)=lmom;
    gsl_matrix_complex_set(iprop,7,6,tmp);

    tmp=gsl_matrix_complex_get(iprop,9,8);
    GSL_REAL(tmp)=lmom;
    gsl_matrix_complex_set(iprop,9,8,tmp);

    tmp=gsl_matrix_complex_get(iprop,11,10);
    GSL_REAL(tmp)=lmom;
    gsl_matrix_complex_set(iprop,11,10,tmp);
    
    // Set gap entries
    GSL_REAL(tmp)=0.0;

    double z;
    z=u->del;
    tmp=gsl_matrix_complex_get(iprop,7,5);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,7,5,tmp);

    tmp=gsl_matrix_complex_get(iprop,9,4);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,9,4,tmp);

    z=-u->del;
    tmp=gsl_matrix_complex_get(iprop,6,3);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,6,3,tmp);

    tmp=gsl_matrix_complex_get(iprop,8,2);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,8,2,tmp);

    z=d->del;
    tmp=gsl_matrix_complex_get(iprop,3,0);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,3,0,tmp);

    tmp=gsl_matrix_complex_get(iprop,11,8);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,11,8,tmp);

    z=-d->del;
    tmp=gsl_matrix_complex_get(iprop,5,1);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,5,1,tmp);

    tmp=gsl_matrix_complex_get(iprop,10,9);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,10,9,tmp);
    
    z=s->del;
    tmp=gsl_matrix_complex_get(iprop,2,0);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,2,0,tmp);

    tmp=gsl_matrix_complex_get(iprop,11,6);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,11,6,tmp);

    z=-s->del;
    tmp=gsl_matrix_complex_get(iprop,4,1);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,4,1,tmp);

    tmp=gsl_matrix_complex_get(iprop,10,7);
    GSL_IMAG(tmp)=z;
    gsl_matrix_complex_set(iprop,10,7,tmp);

    /*
      cout << "iprop: " << endl;
      for(size_t i=0;i<12;i++) {
      for(size_t j=0;j<12;j++) {
      gsl_complex tmpxx=gsl_matrix_complex_get(iprop,i,j);
      cout << i << " " << j << " " << GSL_REAL(tmpxx) << " "
      << GSL_IMAG(tmpxx) << endl;
      }
      }
      exit(-1);
    */
    
    // Set gap matrices for derivatives (here, we need the full
    // matrix, and not just the lower part)
    
    dipdgapu(7,5)+=littlei2;
    dipdgapu(9,4)+=littlei2;
    dipdgapd(3,0)+=littlei2;
    dipdgapd(11,8)+=littlei2;
    dipdgaps(2,0)+=littlei2;
    dipdgaps(11,6)+=littlei2;
    
    dipdgapu(5,7)-=littlei2;
    dipdgapu(4,9)-=littlei2;
    dipdgapd(0,3)-=littlei2;
    dipdgapd(8,11)-=littlei2;
    dipdgaps(0,2)-=littlei2;
    dipdgaps(6,11)-=littlei2;
    
    dipdgapu(6,3)-=littlei2;
    dipdgapu(8,2)-=littlei2;
    dipdgapd(5,1)-=littlei2;
    dipdgapd(10,9)-=littlei2;
    dipdgaps(4,1)-=littlei2;
    dipdgaps(10,7)-=littlei2;
    
    dipdgapu(3,6)+=littlei2;
    dipdgapu(2,8)+=littlei2;
    dipdgapd(1,5)+=littlei2;
    dipdgapd(9,10)+=littlei2;
    dipdgaps(1,4)+=littlei2;
    dipdgaps(7,10)+=littlei2;
    
    gsl_eigen_hermv(iprop,eval,eivec,w);
    
    for(k=0;k<12;k++) {
      egv[k]=gsl_vector_get(eval,k);

      // Select the eigenvector and compute its conjugate
      ubvector_complex eicol(12), eiconj(12);
      for(size_t ij=0;ij<12;ij++) {
	gsl_complex tmp2;
	tmp2=gsl_matrix_complex_get(eivec,ij,k);
	eicol[ij]=std::complex<double>(GSL_REAL(tmp2),GSL_IMAG(tmp2));
	eiconj[ij]=conj(eicol[ij]);
      }
      
      // Compute v+ M v
      std::complex<double> dp;
      ubvector_complex tmpe;
      
      dp=zero2;
      tmpe=prod(dipdgapu,eicol);
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*tmpe[i];
      }
      dedu[k]=dp.real();

      dp=zero2;
      tmpe=prod(dipdgapd,eicol);
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*tmpe[i];
      }
      dedd[k]=dp.real();

      dp=zero2;
      tmpe=prod(dipdgaps,eicol);
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*tmpe[i];
      }
      deds[k]=dp.real();

      ubvector_complex tx(12);

      dp=zero2;
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*mmuu[i]*eicol[i];
      }
      dedmuu[k]=dp.real();

      dp=zero2;
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*mmud[i]*eicol[i];
      }
      dedmud[k]=dp.real();

      dp=zero2;
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*mmus[i]*eicol[i];
      }
      dedmus[k]=dp.real();

      dp=zero2;
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*mmu[i]*eicol[i];
      }
      dedmu[k]=dp.real();

      dp=zero2;
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*mmd[i]*eicol[i];
      }
      dedmd[k]=dp.real();

      dp=zero2;
      for(size_t i=0;i<12;i++) {
	dp+=eiconj[i]*mms[i]*eicol[i];
      }
      dedms[k]=dp.real();
    } 
    
    add_neg=false;
      
  }

  if (add_neg) {
    for(k=0;k<6;k++) {
      egv[k+6]=-egv[k];
      dedmuu[k+6]=-dedmuu[k];
      dedmud[k+6]=-dedmud[k];
      dedmus[k+6]=-dedmus[k];
      dedmu[k+6]=-dedmu[k];
      dedmd[k+6]=-dedmd[k];
      dedms[k+6]=-dedms[k];
      dedu[k+6]=-dedu[k];
      dedd[k+6]=-dedd[k];
      deds[k+6]=-deds[k];
    }
  }
  
  for(size_t j=0;j<12;j++) {
    dedmu3[j]+=dedmuu[j]*drd3+dedmud[j]*dgd3+dedmus[j]*dbd3;
    dedmu8[j]+=dedmuu[j]*drd8+dedmud[j]*dgd8+dedmus[j]*dbd8;
  }
  
  gapped_eigenvalues(u->ms,d->ms,lmom,u->mu+mug,d->mu+mur,s->del,&egv[12],
		     &dedmuu[12],&dedmud[12],&dedmu[12],&dedmd[12],&deds[12]);
  for(size_t j=12;j<16;j++) {
    dedmu3[j]+=dedmuu[j]*dgd3+dedmud[j]*drd3;
    dedmu8[j]+=dedmuu[j]*dgd8+dedmud[j]*drd8;
  }
  gapped_eigenvalues(d->ms,s->ms,lmom,d->mu+mub,s->mu+mug,u->del,&egv[16],
		     &dedmud[16],&dedmus[16],&dedmd[16],&dedms[16],&dedu[16]);
  for(size_t j=16;j<20;j++) {
    dedmu3[j]+=dedmud[j]*dbd3+dedmus[j]*dgd3;
    dedmu8[j]+=dedmud[j]*dbd8+dedmus[j]*dgd8;
  }
  gapped_eigenvalues(u->ms,s->ms,lmom,u->mu+mub,s->mu+mur,d->del,&egv[20],
		     &dedmuu[20],&dedmus[20],&dedmu[20],&dedms[20],&dedd[20]);
  for(size_t j=20;j<24;j++) {
    dedmu3[j]+=dedmuu[j]*dbd3+dedmus[j]*drd3;
    dedmu8[j]+=dedmuu[j]*dbd8+dedmus[j]*drd8;
  }
  
  for(k=12;k<24;k++) {
    egv[k+12]=-egv[k];
    dedmuu[k+12]=-dedmuu[k];
    dedmud[k+12]=-dedmud[k];
    dedmus[k+12]=-dedmus[k];
    dedmu[k+12]=-dedmu[k];
    dedmd[k+12]=-dedmd[k];
    dedms[k+12]=-dedms[k];
    dedu[k+12]=-dedu[k];
    dedd[k+12]=-dedd[k];
    deds[k+12]=-deds[k];
    dedmu3[k+12]=-dedmu3[k];
    dedmu8[k+12]=-dedmu8[k];
  }
  
  return 0;
}
