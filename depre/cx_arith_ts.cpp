/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/cx_arith.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  
  // -----------------------------------------------
  // Test arithmetic operators for gsl_complex
  // -----------------------------------------------

  gsl_complex g1={{1.0,2.0}}, g2={{0.5,-1.0}};
  std::complex<double> c1(1.0,2.0), c2(0.5,-1.0);
  double val=3.0;

  // -----------------------------------------------
  // Complex, complex binary operators

  gsl_complex gsum=g1+g2;
  std::complex<double> csum=c1+c2;
  gsl_complex gdif=g1-g2;
  std::complex<double> cdif=c1-c2;
  gsl_complex gpro=g1*g2;
  std::complex<double> cpro=c1*c2;
  gsl_complex gquo=g1/g2;
  std::complex<double> cquo=c1/c2;

  t.test_rel(GSL_REAL(gsum),csum.real(),1.0e-15,"sumr");
  t.test_rel(GSL_IMAG(gsum),csum.imag(),1.0e-15,"sumi");
  t.test_rel(GSL_REAL(gdif),cdif.real(),1.0e-15,"difr");
  t.test_rel(GSL_IMAG(gdif),cdif.imag(),1.0e-15,"difi");
  t.test_rel(GSL_REAL(gpro),cpro.real(),1.0e-15,"pror");
  t.test_rel(GSL_IMAG(gpro),cpro.imag(),1.0e-15,"proi");
  t.test_rel(GSL_REAL(gquo),cquo.real(),1.0e-15,"quor");
  t.test_rel(GSL_IMAG(gquo),cquo.imag(),1.0e-15,"quoi");
  
  // -----------------------------------------------
  // Complex, double binary operators

  gsl_complex hsum=g1+val;
  gsl_complex isum=val+g1;
  std::complex<double> dsum=c1+val;
  std::complex<double> esum=val+c1;
  gsl_complex hdif=g1-val;
  gsl_complex idif=val-g1;
  std::complex<double> ddif=c1-val;
  std::complex<double> edif=val-c1;
  gsl_complex hpro=g1*val;
  gsl_complex ipro=val*g1;
  std::complex<double> dpro=c1*val;
  std::complex<double> epro=val*c1;
  gsl_complex hquo=g1/val;
  std::complex<double> dquo=c1/val;

  t.test_rel(GSL_REAL(hsum),dsum.real(),1.0e-15,"sumr");
  t.test_rel(GSL_IMAG(hsum),dsum.imag(),1.0e-15,"sumi");
  t.test_rel(GSL_REAL(isum),esum.real(),1.0e-15,"sumr");
  t.test_rel(GSL_IMAG(isum),esum.imag(),1.0e-15,"sumi");
  t.test_rel(GSL_REAL(hdif),ddif.real(),1.0e-15,"difr");
  t.test_rel(GSL_IMAG(hdif),ddif.imag(),1.0e-15,"difi");
  t.test_rel(GSL_REAL(idif),edif.real(),1.0e-15,"difr");
  t.test_rel(GSL_IMAG(idif),edif.imag(),1.0e-15,"difi");
  t.test_rel(GSL_REAL(hpro),dpro.real(),1.0e-15,"pror");
  t.test_rel(GSL_IMAG(hpro),dpro.imag(),1.0e-15,"proi");
  t.test_rel(GSL_REAL(ipro),epro.real(),1.0e-15,"pror");
  t.test_rel(GSL_IMAG(ipro),epro.imag(),1.0e-15,"proi");
  t.test_rel(GSL_REAL(hquo),dquo.real(),1.0e-15,"quor");
  t.test_rel(GSL_IMAG(hquo),dquo.imag(),1.0e-15,"quoi");

  // -----------------------------------------------
  // Complex, complex increment operators

  gsum+=g1;
  csum+=c1;
  gdif-=g1;
  cdif-=c1;
  gpro*=g1;
  cpro*=c1;
  gquo/=g1;
  cquo/=c1;

  t.test_rel(GSL_REAL(gsum),csum.real(),1.0e-15,"sumr2");
  t.test_rel(GSL_IMAG(gsum),csum.imag(),1.0e-15,"sumi2");
  t.test_rel(GSL_REAL(gdif),cdif.real(),1.0e-15,"difr2");
  t.test_rel(GSL_IMAG(gdif),cdif.imag(),1.0e-15,"difi2");
  t.test_rel(GSL_REAL(gpro),cpro.real(),1.0e-15,"pror2");
  t.test_rel(GSL_IMAG(gpro),cpro.imag(),1.0e-15,"proi2");
  t.test_rel(GSL_REAL(gquo),cquo.real(),1.0e-15,"quor2");
  t.test_rel(GSL_IMAG(gquo),cquo.imag(),1.0e-15,"quoi2");

  // -----------------------------------------------
  // Exponential and trig operators, compare O2scl
  // to the STL

  double z=2.0;
  gsl_complex gc={{-2.0,0.0}};
  gsl_complex gd={{-2.0,1.0}};
  std::complex<double> ge(-2.0,1.0);

  t.test_rel(0.0,sqrt(gc).dat[0],1.0e-10,"sqrt real");
  t.test_rel(sqrt(z),sqrt(gc).dat[1],1.0e-10,"sqrt imag");

  t.test_rel(exp(gd).dat[0],exp(ge).real(),1.0e-10,"exp real");
  t.test_rel(exp(gd).dat[1],exp(ge).imag(),1.0e-10,"exp imag");
  t.test_rel(log(gd).dat[0],log(ge).real(),1.0e-10,"log real");
  t.test_rel(log(gd).dat[1],log(ge).imag(),1.0e-10,"log imag");
  t.test_rel(log10(gd).dat[0],log10(ge).real(),1.0e-10,"log10 real");
  t.test_rel(log10(gd).dat[1],log10(ge).imag(),1.0e-10,"log10 imag");

  t.test_rel(sin(gd).dat[0],sin(ge).real(),1.0e-10,"sin real");
  t.test_rel(sin(gd).dat[1],sin(ge).imag(),1.0e-10,"sin imag");
  t.test_rel(cos(gd).dat[0],cos(ge).real(),1.0e-10,"cos real");
  t.test_rel(cos(gd).dat[1],cos(ge).imag(),1.0e-10,"cos imag");
  t.test_rel(tan(gd).dat[0],tan(ge).real(),1.0e-10,"tan real");
  t.test_rel(tan(gd).dat[1],tan(ge).imag(),1.0e-10,"tan imag");

  t.test_rel(sinh(gd).dat[0],sinh(ge).real(),1.0e-10,"sinh real");
  t.test_rel(sinh(gd).dat[1],sinh(ge).imag(),1.0e-10,"sinh imag");
  t.test_rel(cosh(gd).dat[0],cosh(ge).real(),1.0e-10,"cosh real");
  t.test_rel(cosh(gd).dat[1],cosh(ge).imag(),1.0e-10,"cosh imag");
  t.test_rel(tanh(gd).dat[0],tanh(ge).real(),1.0e-10,"tanh real");
  t.test_rel(tanh(gd).dat[1],tanh(ge).imag(),1.0e-10,"tanh imag");

  // -----------------------------------------------
  // The remaining tests only compare O2scl with GSL
  // -----------------------------------------------

  gsl_complex gf={{-2.0,1.0}};

  t.test_rel(sec(gd).dat[0],gsl_complex_sec(gf).dat[0],1.0e-10,"sec r");
  t.test_rel(sec(gd).dat[1],gsl_complex_sec(gf).dat[1],1.0e-10,"sec i");
  t.test_rel(csc(gd).dat[0],gsl_complex_csc(gf).dat[0],1.0e-10,"csc r");
  t.test_rel(csc(gd).dat[1],gsl_complex_csc(gf).dat[1],1.0e-10,"csc i");
  t.test_rel(cot(gd).dat[0],gsl_complex_cot(gf).dat[0],1.0e-10,"cot r");
  t.test_rel(cot(gd).dat[1],gsl_complex_cot(gf).dat[1],1.0e-10,"cot i");

  t.test_rel(sech(gd).dat[0],gsl_complex_sech(gf).dat[0],1.0e-10,"sech r");
  t.test_rel(sech(gd).dat[1],gsl_complex_sech(gf).dat[1],1.0e-10,"sech i");
  t.test_rel(csch(gd).dat[0],gsl_complex_csch(gf).dat[0],1.0e-10,"csch r");
  t.test_rel(csch(gd).dat[1],gsl_complex_csch(gf).dat[1],1.0e-10,"csch i");
  t.test_rel(coth(gd).dat[0],gsl_complex_coth(gf).dat[0],1.0e-10,"coth r");
  t.test_rel(coth(gd).dat[1],gsl_complex_coth(gf).dat[1],1.0e-10,"coth i");

  // -----------------------------------------------
  
  t.test_rel(asin(gd).dat[0],gsl_complex_arcsin(gf).dat[0],1.0e-10,"asin r");
  t.test_rel(asin(gd).dat[1],gsl_complex_arcsin(gf).dat[1],1.0e-10,"asin i");
  t.test_rel(acos(gd).dat[0],gsl_complex_arccos(gf).dat[0],1.0e-10,"acos r");
  t.test_rel(acos(gd).dat[1],gsl_complex_arccos(gf).dat[1],1.0e-10,"acos i");
  t.test_rel(atan(gd).dat[0],gsl_complex_arctan(gf).dat[0],1.0e-10,"atan r");
  t.test_rel(atan(gd).dat[1],gsl_complex_arctan(gf).dat[1],1.0e-10,"atan i");

  t.test_rel(asinh(gd).dat[0],gsl_complex_arcsinh(gf).dat[0],
	     1.0e-10,"asinh r");
  t.test_rel(asinh(gd).dat[1],gsl_complex_arcsinh(gf).dat[1],
	     1.0e-10,"asinh i");
  t.test_rel(acosh(gd).dat[0],gsl_complex_arccosh(gf).dat[0],
	     1.0e-10,"acosh r");
  t.test_rel(acosh(gd).dat[1],gsl_complex_arccosh(gf).dat[1],
	     1.0e-10,"acosh i");
  t.test_rel(atanh(gd).dat[0],gsl_complex_arctanh(gf).dat[0],
	     1.0e-10,"atanh r");
  t.test_rel(atanh(gd).dat[1],gsl_complex_arctanh(gf).dat[1],
	     1.0e-10,"atanh i");

  // -----------------------------------------------

  t.test_rel(asec(gd).dat[0],gsl_complex_arcsec(gf).dat[0],1.0e-10,"asec r");
  t.test_rel(asec(gd).dat[1],gsl_complex_arcsec(gf).dat[1],1.0e-10,"asec i");
  t.test_rel(acsc(gd).dat[0],gsl_complex_arccsc(gf).dat[0],1.0e-10,"acsc r");
  t.test_rel(acsc(gd).dat[1],gsl_complex_arccsc(gf).dat[1],1.0e-10,"acsc i");
  t.test_rel(acot(gd).dat[0],gsl_complex_arccot(gf).dat[0],1.0e-10,"acot r");
  t.test_rel(acot(gd).dat[1],gsl_complex_arccot(gf).dat[1],1.0e-10,"acot i");

  t.test_rel(asech(gd).dat[0],gsl_complex_arcsech(gf).dat[0],
	     1.0e-10,"asech r");
  t.test_rel(asech(gd).dat[1],gsl_complex_arcsech(gf).dat[1],
	     1.0e-10,"asech i");
  t.test_rel(acsch(gd).dat[0],gsl_complex_arccsch(gf).dat[0],
	     1.0e-10,"acsch r");
  t.test_rel(acsch(gd).dat[1],gsl_complex_arccsch(gf).dat[1],
	     1.0e-10,"acsch i");
  t.test_rel(acoth(gd).dat[0],gsl_complex_arccoth(gf).dat[0],
	     1.0e-10,"acoth r");
  t.test_rel(acoth(gd).dat[1],gsl_complex_arccoth(gf).dat[1],
	     1.0e-10,"acoth i");

  // -----------------------------------------------
  
  t.test_rel(asin_real(val).dat[0],
	     gsl_complex_arcsin_real(val).dat[0],1.0e-10,"asin r");
  t.test_rel(asin_real(val).dat[1],
	     gsl_complex_arcsin_real(val).dat[1],1.0e-10,"asin i");
  t.test_rel(acos_real(val).dat[0],
	     gsl_complex_arccos_real(val).dat[0],1.0e-10,"acos r");
  t.test_rel(acos_real(val).dat[1],
	     gsl_complex_arccos_real(val).dat[1],1.0e-10,"acos i");

  t.test_rel(acosh_real(val).dat[0],
	     gsl_complex_arccosh_real(val).dat[0],1.0e-10,"acosh r");
  t.test_rel(acosh_real(val).dat[1],
	     gsl_complex_arccosh_real(val).dat[1],1.0e-10,"acosh i");
  t.test_rel(atanh_real(val).dat[0],
	     gsl_complex_arctanh_real(val).dat[0],1.0e-10,"atanh r");
  t.test_rel(atanh_real(val).dat[1],
	     gsl_complex_arctanh_real(val).dat[1],1.0e-10,"atanh i");

  // -----------------------------------------------

  t.test_rel(asec_real(val).dat[0],
	     gsl_complex_arcsec_real(val).dat[0],1.0e-10,"asec r");
  t.test_rel(asec_real(val).dat[1],
	     gsl_complex_arcsec_real(val).dat[1],1.0e-10,"asec i");
  t.test_rel(acsc_real(val).dat[0],
	     gsl_complex_arccsc_real(val).dat[0],1.0e-10,"acsc r");
  t.test_rel(acsc_real(val).dat[1],
	     gsl_complex_arccsc_real(val).dat[1],1.0e-10,"acsc i");

  // -----------------------------------------------

  t.test_rel(abs(gd),gsl_complex_abs(gd),1.0e-10,"abs");
  t.test_rel(arg(gd),gsl_complex_arg(gd),1.0e-10,"arg");
  t.test_rel(abs2(gd),gsl_complex_abs2(gd),1.0e-10,"abs2");
  t.test_rel(logabs(gd),gsl_complex_logabs(gd),1.0e-10,"logabs");

  t.test_abs(sqrt_real(-2.0).dat[0],gsl_complex_sqrt_real(-2.0).dat[0],
	     1.0e-10,"sqrt real r");
  t.test_rel(sqrt_real(-2.0).dat[1],gsl_complex_sqrt_real(-2.0).dat[1],
	     1.0e-10,"sqrt real i");

  t.test_rel(conjugate(gd).dat[0],gsl_complex_conjugate(gd).dat[0],
	     1.0e-10,"conjugate r");
  t.test_rel(conjugate(gd).dat[1],gsl_complex_conjugate(gd).dat[1],
	     1.0e-10,"conjugate i");
  
  t.test_rel(pow_real(gd,-2.0).dat[0],gsl_complex_pow_real(gd,-2.0).dat[0],
	     1.0e-10,"pow_real r");
  t.test_rel(pow_real(gd,-2.0).dat[1],gsl_complex_pow_real(gd,-2.0).dat[1],
	     1.0e-10,"pow_real i");

  t.test_rel(log_b(gd,gf).dat[0],gsl_complex_log_b(gd,gf).dat[0],
	     1.0e-10,"log_b r");
  t.test_rel(log_b(gd,gf).dat[1],gsl_complex_log_b(gd,gf).dat[1],
	     1.0e-10,"log_b i");

  t.report();
  
  return 0;
}

