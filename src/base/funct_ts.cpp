/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/funct.h>
#include <o2scl/test_mgr.h>
#include <o2scl/err_hnd.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

class fmc {
public:
  template<class fp_t> fp_t func(fp_t x) {
    return log(1+x);
  }
};

int main(void) {
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  funct_string f("pi*r^2","r");
  f.set_parm("pi",o2scl_const::pi);
  for(double r=1.0;r<=2.0;r+=0.1) {
    t.test_rel(f(r),o2scl_const::pi*r*r,1.0e-12,"funct_string");
  }

#ifdef O2SCL_PYTHON

  o2scl_settings.py_init();
  o2scl_settings.add_python_path("./");
  {
    // We use the brackets to force the funct_python
    // destructor to run before py_final()
    funct_python fp("funct_test","fun");
    t.test_rel(fp(2.0),o2scl_const::pi*2.0,1.0e-12,"funct_python");
  }
  {
    // We use the brackets to force the funct_python
    // destructor to run before py_final()
    funct_python_method fp("funct_test","c","fun2",2);
    t.test_rel(fp(2.0),o2scl_const::pi*4.0,1.0e-12,"funct_python");
  }
  o2scl_settings.py_final();
  
#endif


  fmc f2;
  funct f2_d=std::bind(std::mem_fn<double(double)>(&fmc::func<double>),
                       &f2,std::placeholders::_1);
  funct_ld f2_ld=std::bind(std::mem_fn<long double(long double)>
                           (&fmc::func<long double>),
                           &f2,std::placeholders::_1);
  funct_cdf25 f2_cdf25=std::bind(std::mem_fn<boost::multiprecision::number<
                                 boost::multiprecision::cpp_dec_float<25>>
                                 (boost::multiprecision::number<
                                  boost::multiprecision::cpp_dec_float<25>>)>
                                 (&fmc::func<boost::multiprecision::number<
                                  boost::multiprecision::cpp_dec_float<25>>>),
                                 &f2,std::placeholders::_1);
  funct_cdf35 f2_cdf35=std::bind(std::mem_fn<boost::multiprecision::number<
                                 boost::multiprecision::cpp_dec_float<35>>
                                 (boost::multiprecision::number<
                                  boost::multiprecision::cpp_dec_float<35>>)>
                                 (&fmc::func<boost::multiprecision::number<
                                  boost::multiprecision::cpp_dec_float<35>>>),
                                 &f2,std::placeholders::_1);
  funct_cdf50 f2_cdf50=std::bind(std::mem_fn<boost::multiprecision::number<
                                 boost::multiprecision::cpp_dec_float<50>>
                                 (boost::multiprecision::number<
                                  boost::multiprecision::cpp_dec_float<50>>)>
                                 (&fmc::func<boost::multiprecision::number<
                                  boost::multiprecision::cpp_dec_float<50>>>),
                                 &f2,std::placeholders::_1);
  funct_cdf100 f2_cdf100=std::bind
    (std::mem_fn<boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<100>>
     (boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<100>>)>
     (&fmc::func<boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<100>>>),
     &f2,std::placeholders::_1);
  funct_multip_wrapper fmw(f2_d,f2_ld,f2_cdf25,f2_cdf35,
                           f2_cdf50,f2_cdf100);
  funct_multip fm(fmw);
  fm.verbose=1;
  t.test_rel(log1p(1.0e-4),fm(1.0e-4),1.0e-15,"funct_multip 1");
  t.test_rel<long double>(log1p(1.0e-4L),fm(1.0e-4L),1.0e-18,"funct_multip 2");

  double val, err;
  funct_multip2 fm2;
  fm2.eval_tol_err([f2](auto &&t) mutable
  { return f2.func(std::forward<decltype(t)>(t)); },1.0e-4,val,err);
  cout << fm(1.0e-4) << " " << val << endl;
  
  t.report();
  return 0;
}

