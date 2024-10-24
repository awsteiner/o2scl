/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <o2scl/funct.h>
#include <o2scl/funct_multip.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/test_mgr.h>

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std;
using namespace o2scl;

double gfn(double x) {
  return atan((x-0.2)*4)*(1.0+sin((x-0.2)*50.0)/1.1);
}

double gfn2(double x) {
  return tanh(100.0*(x-0.2));
}

class cl {
public:
  double mfn(double x) {
    return gfn(x);
  }
};

    
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

template<class fp_t> fp_t cbrt_fun(fp_t x) {
  return x*x*x-5;
}

class cl_ld {

public:

  long double mfn(long double x) {
    long double one=1;
    long double five=5;
    long double ten=10;
    return atan((x-one/five)*4)*(one+sin((x-one/five)*five*ten)/
                                 (one+one/ten));
  }
};

class cl_cdf {

public:

  cpp_dec_float_50 mfn(cpp_dec_float_50 x) {
    cpp_dec_float_50 one=1;
    cpp_dec_float_50 five=5;
    cpp_dec_float_50 ten=10;
    return atan((x-one/five)*4)*(one+sin((x-one/five)*five*ten)/
                                 (one+one/ten));
  }
};

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  cl acl;
  double a, b, r;

  funct fmg1=gfn;
  funct fmg2=gfn2;
  funct_gsl fg1(fmg1);
  funct_gsl fg2(fmg2);

  for(size_t k=0;k<2;k++) {

    {
      if (k==0) {
        cout << "Using GSL with the first function:" << endl;
      } else {
        cout << "Using GSL with the second function:" << endl;
      }
      int status;
      int iter=0, max_iter=100;

      a=-1.0;
      b=1.0;
    
      const gsl_root_fsolver_type *T=gsl_root_fsolver_brent;
      gsl_root_fsolver *s=gsl_root_fsolver_alloc(T);
      
      if (k==0) gsl_root_fsolver_set(s,&fg1,a,b);
      else gsl_root_fsolver_set(s,&fg2,a,b);
     
      do {

	iter++;
	status=gsl_root_fsolver_iterate(s);
	r=gsl_root_fsolver_root(s);
	a=gsl_root_fsolver_x_lower(s);
	b=gsl_root_fsolver_x_upper(s);
	status=gsl_root_test_interval(a,b,0,1.0e-6);

	cout << iter << " " << a << " " << r << " " << b << endl;

      } while (status==gsl_continue && iter<max_iter);
     
      gsl_root_fsolver_free(s);

      if (k==0) {
	cout << r << " " << gfn(r) << endl;
      } else {
	cout << r << " " << gfn2(r) << endl;
      }
      t.test_rel(r,0.2,1.0e-10,"1");
      cout << endl;
      
    }

    if (k==0) {
      cout << "Using O2scl set() and iterate() with the first function:"
           << endl;
    } else {
      cout << "Using O2scl set() and iterate() with the second function:"
           << endl;
    }
    
    // Using set() and iterate() with a function pointer
    typedef double (*gfnt)(double);
    root_brent_gsl<gfnt> grb1;
    gfnt gfnv=&gfn;
    if (k==1) gfnv=&gfn2;
    a=-1.0;
    b=1.0;
    grb1.set(gfnv,a,b);
    int status=gsl_continue;
    for(size_t iter=0;iter<100 && status==gsl_continue;iter++) {

      grb1.iterate(gfnv);

      r=grb1.get_root();
      a=grb1.get_lower();
      b=grb1.get_upper();
      status=gsl_root_test_interval(a,b,0,1.0e-6);

      cout << iter << " " << a << " " << r << " " << b << endl;
    }
    
    if (k==0) {
      cout << r << " " << gfn(r) << endl;
    } else {
      cout << r << " " << gfn2(r) << endl;
    }
    t.test_rel(grb1.get_root(),0.2,1.0e-10,"1");
    cout << endl;

  }

  cout << "Using O2scl, a member function, and the simplified interface,\n"
       << "  with the first function:" << endl;
  funct fmf=std::bind(std::mem_fn<double(double)>
                      (&cl::mfn),&acl,std::placeholders::_1);
  root_brent_gsl<> grb2;
  a=-1.0;
  b=1.0;
  grb2.store_funcs=true;
  grb2.solve_bkt(a,b,fmf);
  grb2.store_funcs=false;
  t.test_rel(a,0.2,1.0e-10,"1");
  cout << a << " " << gfn(a) << endl;
  cout << endl;

  cout << "x                      y" << endl;
  for(size_t j=0;j<grb2.tab.get_nlines();j++) {
    cout << dtos(grb2.tab.get("x",j),0) << " ";
    cout << dtos(grb2.tab.get("y",j),0) << endl;
  }
  cout << endl;

  cout << "Now using long double precision:" << endl;
  cl_ld acl_ld;
  funct_ld fmf_ld=std::bind(std::mem_fn<long double(long double)>
                            (&cl_ld::mfn),&acl_ld,std::placeholders::_1);
  root_brent_gsl<funct_ld,long double> grb2_ld;
  long double a_ld=-1, b_ld=1;
  grb2_ld.solve_bkt(a_ld,b_ld,fmf_ld);
  cout << dtos(a_ld,0) << " " << dtos(acl_ld.mfn(a_ld),0) << endl;
  t.test_rel<long double>(10*a_ld,2,
                          static_cast<long double>(grb2_ld.tol_rel),
                          "long double 1");
  cout << endl;

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  cout << "Now using 50-digit precision:" << endl;
  cl_cdf acl_cdf;
  funct_cdf50 fmf_cdf=std::bind(std::mem_fn<cpp_dec_float_50(cpp_dec_float_50)>
                                (&cl_cdf::mfn),&acl_cdf,std::placeholders::_1);
  root_brent_gsl<funct_cdf50,cpp_dec_float_50> grb2_cdf;
  cpp_dec_float_50 a_cdf=-1, b_cdf=1;
  grb2_cdf.solve_bkt(a_cdf,b_cdf,fmf_cdf);
  cout << dtos(a_cdf,0) << " " << dtos(acl_cdf.mfn(a_cdf),0) << endl;
  t.test_rel_boost<cpp_dec_float_50>(10*a_cdf,2,
                          static_cast<cpp_dec_float_50>(grb2_cdf.tol_rel),
                          "50-digit 1");
  cout << endl;
  
  cout << "Using adaptive multiprecision with a simple function and a\n"
       << "  lambda expression:" << endl;
  double am=1.0, bm=5.0, valm, errm;
  grb2.verbose=1;
  int amret=grb2.solve_bkt_multip(am,bm,[](auto &&tx) mutable
  { return cbrt_fun(tx); },errm);
  cout << dtos(am,0) << " " << dtos(cbrt(5.0),0) << endl;
  t.test_rel(am,cbrt(5.0),2.0e-15,"multiprecision");
#endif
  
  t.report();
  return 0;
}

