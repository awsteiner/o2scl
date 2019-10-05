/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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

#include <o2scl/funct.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double tx, double &a) {
  return -cos(1.0/(tx+a))/(a+tx)/(a+tx);
}

long double testfun_ld(long double tx, long double &a) {
return -cos(1.0/(tx+a))/(a+tx)/(a+tx);
}

double sin_recip(double x) {
return sin(1.0/(-x+0.01))*pow(-x+0.01,-2.0);
}

long double sin_recip_ld(long double x) {
return sin(1.0/(-x+0.01))*pow(-x+0.01,-2.0);
}

#ifdef O2SCL_LD_TYPES

/*
  __float128 testfun_f128(__float128 tx, __float128 &a) {
  return -cosq(1.0/(tx+a))/(a+tx)/(a+tx);
  }
  
  __float128 sin_recip_f128(__float128 x) {
  return sinq(1.0/(-x+0.01))*pow(-x+0.01,-2.0);
  }
*/

cpp_dec_float_50 testfun_cdf(cpp_dec_float_50 tx, cpp_dec_float_50 &a) {
return -cos(1.0/(tx+a))/(a+tx)/(a+tx);
}

cpp_dec_float_50 sin_recip_cdf(cpp_dec_float_50 x) {
return sin(1.0/(-x+0.01))*pow(-x+0.01,-2.0);
}

typedef std::function<__float128(__float128)> funct_f128;
typedef std::function<cpp_dec_float_50(cpp_dec_float_50)> funct_cdf;

#endif

int main(void) {
  
test_mgr t;
t.set_output_level(2);

cout.setf(ios::scientific);
cout.precision(10);

{
double a, calc, exact, diff, ei;
inte_adapt_cern<funct> cg;
  
a=0.01;
funct tf=std::bind(testfun,std::placeholders::_1,a);

cg.integ_err(tf,0.0,1.0,calc,ei);
exact=sin(1.0/(1.0+a))-sin(1.0/a);
t.test_rel(calc,exact,1.0e-8,"inte_adapt_cern");
diff=fabs(calc-exact);
cout << calc << " " << exact << " " << diff << " " << ei << endl;

// This is a nasty function and takes many subdivisions (68)
cout << cg.get_nsubdivisions() << endl;
size_t n=cg.get_nsubdivisions();
typedef boost::numeric::ublas::vector<double> ubvector;
ubvector xlo(n), xhi(n), val(n), err(n);
cg.get_subdivisions(xlo,xhi,val,err);
for(size_t i=0;i<n;i+=10) {
cout << xlo[i] << " " << xhi[i] << " ";
cout.setf(ios::showpos);
cout << val[i] << " ";
cout.unsetf(ios::showpos);
cout << err[i] << endl;
}

long double a_ld=0.01L, calc_ld, ei_ld, diff_ld;
inte_adapt_cern<funct_ld,100,long double,
		inte_gauss56_coeffs_long_double> cg_ld;
funct_ld tf_ld=std::bind(testfun_ld,std::placeholders::_1,a_ld);
long double exact_ld=sin(1.0/(1.0+a_ld))-sin(1.0/a_ld);
cg_ld.integ_err(tf_ld,0.0L,1.0L,calc_ld,ei_ld);
t.test_rel(calc_ld,exact_ld,1.0e-8L,"inte_adapt_cern_ld");
diff_ld=fabs(calc_ld-exact_ld);
cout << calc_ld << " " << exact_ld << " " << diff_ld << " "
<< ei_ld << endl;
}

{
double calc, ei, diff;
// Test qagil_cern with double precision

inte_qagil_cern<funct> iqc;
double exact=1.0-cos(100.0/101.0);
funct tf2=std::bind(sin_recip,std::placeholders::_1);
iqc.integ_err(tf2,0.0,-1.0,calc,ei);
diff=fabs(calc-exact);
cout << calc << " " << exact << " " << diff << " " << ei << endl;
  
// Test qagil_cern with long double precision
inte_qagil_cern<funct_ld,long double,
		inte_gauss56_coeffs_long_double> iqc_ld;
long double exact_ld=1.0-cos(100.0/101.0);
funct_ld tf2_ld=std::bind(sin_recip_ld,std::placeholders::_1);
long double calc_ld, ei_ld;
iqc_ld.integ_err(tf2_ld,0.0,-1.0,calc_ld,ei_ld);
long double diff_ld=fabs(calc_ld-exact_ld);
cout << calc_ld << " " << exact_ld << " "
<< diff_ld << " " << ei_ld << endl;

#ifdef O2SCL_LD_TYPES

// Test qagil_cern with __float128 precision
/*
  inte_qagil_cern<funct_f128,__float128,
  inte_gauss56_coeffs_float128> iqc_f128;
  __float128 exact_f128=1.0-cos(100.0/101.0);
  funct_f128 tf2_f128=std::bind(sin_recip_f128,std::placeholders::_1);
  __float128 calc_f128, ei_f128;
  iqc_f128.integ_err(tf2_f128,0.0,-1.0,calc_f128,ei_f128);
  __float128 diff_f128=fabsl(calc_f128-exact_f128);
  cout << (long double)calc_f128 << " "
  << (long double)exact_f128 << " "
  << (long double)diff_f128 << " "
  << (long double)ei_f128 << endl;
*/

// Test qagil_cern with cpp_dec_float_50 precision
inte_qagil_cern<funct_cdf,cpp_dec_float_50,
		inte_gauss56_coeffs_cpp_dec_float_50> iqc_cdf;
cpp_dec_float_50 exact_cdf=1.0-cos(100.0/101.0);
funct_cdf tf2_cdf=std::bind(sin_recip_cdf,std::placeholders::_1);
cpp_dec_float_50 calc_cdf, ei_cdf;
iqc_cdf.integ_err(tf2_cdf,0.0,-1.0,calc_cdf,ei_cdf);
cpp_dec_float_50 diff_cdf=fabs(calc_cdf-exact_cdf);
cout << calc_cdf << " " << exact_cdf << " "
<< diff_cdf << " " << ei_cdf << endl;

#endif
}
  
t.report();
return 0;
}

