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

#ifdef O2SCL_LD_TYPES
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

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

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef std::function<cpp_dec_float_50(cpp_dec_float_50)> funct_cdf;

cpp_dec_float_50 testfun_cdf(cpp_dec_float_50 tx, cpp_dec_float_50 &a) {
  cpp_dec_float_50 one=1.0;
  return -cos(one/(tx+a))/(a+tx)/(a+tx);
}

cpp_dec_float_50 sin_recip_cdf(cpp_dec_float_50 x) {
  cpp_dec_float_50 one=1;
  cpp_dec_float_50 hundred=100;
  return sin(one/(-x+one/hundred))*pow(-x+one/hundred,-one-one);
}

#endif

template<class func_t=funct, class fp_t=double,
	 class def_inte_t, size_t nsub>
void test_iac(test_mgr &t, func_t &f, fp_t acc,
	      std::string comment, fp_t &diff, bool output_sub=false) {

  cout << comment << ":\n  ";
  fp_t one=1;
  fp_t ten=10;
  fp_t hundred=100;
  fp_t a=one/hundred, calc, ei;
  fp_t exact=sin(one/(one+a))-sin(one/a);
  
  inte_adapt_cern<func_t,def_inte_t,nsub,fp_t> iac;
  iac.tol_rel=acc;
  iac.tol_abs=acc;
  iac.integ_err(f,0.0,one,calc,ei);
  diff=fabs(calc-exact);
  cout << calc << " " << exact << " " << diff << " "
       << ei << endl;
  cout << "subdivisions: ";
  cout << iac.get_nsubdivisions() << endl;
  cout << endl;

  if (output_sub) {
    size_t n=iac.get_nsubdivisions();
    typedef boost::numeric::ublas::vector<fp_t> ubvector;
    ubvector xlo(n), xhi(n), val(n), err(n);
    iac.get_subdivisions(xlo,xhi,val,err);
    cout << "xlo              xhi               ";
    cout << "val              err              " << endl;;
    for(size_t i=0;i<n;i+=10) {
      cout << xlo[i] << " " << xhi[i] << " ";
      cout.setf(ios::showpos);
      cout << val[i] << " ";
      cout.unsetf(ios::showpos);
      cout << err[i] << endl;
    }
    cout << endl;
  }

  return;
}

int main(void) {
  
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);
  cout.precision(10);

  {
    double a=0.01, diff;
    funct tf=std::bind(testfun,std::placeholders::_1,a);
    test_iac<funct,double,
	     inte_gauss56_cern<funct,double,
			       inte_gauss56_coeffs_double>,100>
      (t,tf,1.0e-8,"iac, double, testfun",diff,true);
    t.test_abs<double>(diff,0.0,1.0e-7,"inte_adapt_cern");

    test_iac<funct,double,
	     inte_newton_cotes_open<funct,double>,1000>
      (t,tf,1.0e-8,"iac, NCO, double, testfun",diff);
    t.test_abs<double>(diff,0.0,1.0e-7,"inte_adapt_cern");

#ifdef O2SCL_LD_TYPES
    
    long double a_ld=0.01L, diff_ld;
    funct_ld tf_ld=std::bind(testfun_ld,std::placeholders::_1,a_ld);
    test_iac<funct_ld,long double,
	     inte_gauss56_cern<funct_ld,long double,
			       inte_gauss56_coeffs_long_double>,1000>
      (t,tf_ld,1.0e-15,"iac, long double, testfun",diff_ld);
    t.test_abs<long double>(diff_ld,0.0,1.0e-14,"inte_adapt_cern_ld");
    
    test_iac<funct_ld,long double,
	     inte_newton_cotes_open<funct_ld,long double>,10000>
      (t,tf_ld,1.0e-15,"iac, NCO, long double, testfun",diff_ld);
    t.test_abs<long double>(diff_ld,0.0,1.0e-14,"inte_adapt_cern_ld");
    
    cout << "iac, cpp_dec_float_50 prec, testfun:\n  ";
    cpp_dec_float_50 one=1.0, diff_cdf;
    cpp_dec_float_50 hundred=100.0;
    cpp_dec_float_50 a_cdf=one/hundred;
    funct_cdf50 tf_cdf=std::bind(testfun_cdf,std::placeholders::_1,a_cdf);
    test_iac<funct_cdf50,cpp_dec_float_50,
	     inte_gauss56_cern<funct_cdf50,cpp_dec_float_50,
			       inte_gauss56_coeffs_cpp_dec_float_50>,10000>
      (t,tf_cdf,1.0e-30,"iac, cpp_dec_float_50, testfun",diff_cdf);

    // AWS 10/8/19: this segfaults I don't know why yet
    //test_iac<funct_cdf50,cpp_dec_float_50,
    //inte_newton_cotes_open<funct_cdf50,cpp_dec_float_50>,100000>
    //(t,tf_cdf,1.0e-30,"iac, NCO, cpp_dec_float_50, testfun",diff_cdf);
    
#endif
    
  }

  if (false) {
    
    double calc, ei, diff;
    // Test qagil_cern with double precision
    
    cout << "iqc, double prec, sin_recip:\n  ";
    inte_qagil_cern<funct> iqc;
    double exact=1.0-cos(100.0/101.0);
    funct tf2=std::bind(sin_recip,std::placeholders::_1);
    iqc.integ_err(tf2,0.0,-1.0,calc,ei);
    diff=fabs(calc-exact);
    cout << calc << " " << exact << " " << diff << " " << ei << endl;
    cout << endl;
  
    // Test qagil_cern with long double precision
    cout << "iqc, long double prec, sin_recip:\n  ";
    inte_qagil_cern<funct_ld,
		    inte_gauss56_cern<funct_ld,long double,
				      inte_gauss56_coeffs_long_double>,
		    100,long double> iqc_ld;
    long double exact_ld=1.0L-cos(100.0L/101.0L);
    funct_ld tf2_ld=std::bind(sin_recip_ld,std::placeholders::_1);
    long double calc_ld, ei_ld;
    iqc_ld.integ_err(tf2_ld,0.0,-1.0,calc_ld,ei_ld);
    long double diff_ld=fabs(calc_ld-exact_ld);
    cout << calc_ld << " " << exact_ld << " "
	 << diff_ld << " " << ei_ld << endl;
    cout << endl;

#ifdef O2SCL_LD_TYPES

    // Test qagil_cern with double precision
    cout << "iqc, double newton-cotes, testfun:\n  ";
    inte_adapt_cern<funct,inte_newton_cotes<funct,double>,
		    100000,double> iqc_d3;
    double a_d3=0.01;
    double exact_d3=sin(1.0/(1.0+a_d3))-sin(1.0/a_d3);
    funct tf2_d3=std::bind(testfun,std::placeholders::_1,a_d3);
    double calc_d3, ei_d3;
    iqc_d3.integ_err(tf2_d3,0.0,1.0,calc_d3,ei_d3);
    double diff_d3=fabs(calc_d3-exact_d3);
    cout << calc_d3 << " " << exact_d3 << " "
	 << diff_d3 << " " << ei_d3 << endl;
    cout << endl;

    // Test qagil_cern with double precision
    cout << "iqc, double newton-cotes2, testfun:\n  ";
    inte_adapt_cern<funct,inte_newton_cotes2<funct,double>,
		    100000,double> iqc_d4;
    double a_d4=0.01;
    double exact_d4=sin(1.0/(1.0+a_d4))-sin(1.0/a_d4);
    funct tf2_d4=std::bind(testfun,std::placeholders::_1,a_d4);
    double calc_d4, ei_d4;
    iqc_d4.integ_err(tf2_d4,0.0,1.0,calc_d4,ei_d4);
    double diff_d4=fabs(calc_d4-exact_d4);
    cout << calc_d4 << " " << exact_d4 << " "
	 << diff_d4 << " " << ei_d4 << endl;
    cout << endl;

    // Test qagil_cern with double precision
    cout << "iqc, double newton-cotes_open, testfun:\n  ";
    inte_adapt_cern<funct,inte_newton_cotes_open<funct,double>,
		    100000,double> iqc_d5;
    double a_d5=0.01;
    double exact_d5=sin(1.0/(1.0+a_d5))-sin(1.0/a_d5);
    funct tf2_d5=std::bind(testfun,std::placeholders::_1,a_d5);
    double calc_d5, ei_d5;
    iqc_d5.integ_err(tf2_d5,0.0,1.0,calc_d5,ei_d5);
    double diff_d5=fabs(calc_d5-exact_d5);
    cout << calc_d5 << " " << exact_d5 << " "
	 << diff_d5 << " " << ei_d5 << endl;
    cout << endl;

    // Test qagil_cern with double precision
    cout << "iqc, double newton-cotes, sin_recip:\n  ";
    inte_qagil_cern<funct,inte_newton_cotes<funct,double>,
		    1000,double> iqc_d2;
    double exact_d2=1.0-cos(100.0/101.0);
    funct tf2_d2=std::bind(sin_recip,std::placeholders::_1);
    double calc_d2, ei_d2;
    iqc_d2.verbose=1;
    iqc_d2.integ_err(tf2_d2,0.0,-1.0,calc_d2,ei_d2);
    double diff_d2=fabs(calc_d2-exact_d2);
    cout << calc_d2 << " " << exact_d2 << " "
	 << diff_d2 << " " << ei_d2 << endl;
    cout << endl;

    // Test qagil_cern with cpp_dec_float_50 precision
    cout << "iqc, cdf 50 prec, sin_recip:\n  ";
    inte_qagil_cern<funct_cdf,
		    inte_newton_cotes<funct_cdf,cpp_dec_float_50>,
		    1000,cpp_dec_float_50> iqc_cdf;
    cpp_dec_float_50 one=1.0;
    cpp_dec_float_50 hundred=100.0;
    cpp_dec_float_50 exact_cdf=one-cos(hundred/(hundred+one));
    funct_cdf tf2_cdf=std::bind(sin_recip_cdf,std::placeholders::_1);
    cpp_dec_float_50 calc_cdf, ei_cdf;
    iqc_cdf.integ_err(tf2_cdf,0.0,-one,calc_cdf,ei_cdf);
    cpp_dec_float_50 diff_cdf=fabs(calc_cdf-exact_cdf);
    cout << calc_cdf << " " << exact_cdf << " "
	 << diff_cdf << " " << ei_cdf << endl;
    cout << endl;

#endif
  }
  
  t.report();
  return 0;
}

