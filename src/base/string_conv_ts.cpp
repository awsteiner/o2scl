/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/misc.h>
#include <o2scl/test_mgr.h>
#include <climits>

using namespace std;
using namespace o2scl;

#include <boost/multiprecision/cpp_dec_float.hpp>

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  double z1=0.0;
  double z2=-0.0;
  t.test_gen(has_minus_sign(&z1)==false,"hms1");
  t.test_gen(has_minus_sign(&z2)==true,"hms1");

  t.test_gen(stob(" true")==1,"stob1");
  t.test_gen(stob(" false")==0,"stob2");
  t.test_gen(stob("True")==1,"stob3");
  t.test_gen(stob("1")==1,"stob4");
  t.test_gen(stob(" 1")==1,"stob5");
  t.test_gen(stob(" 0")==0,"stob6");
  t.test_gen(stob("-1")==0,"stob1");

  t.test_gen(size_of_exponent(1.0e-111)==3,"soe1");
  t.test_gen(size_of_exponent(1.0e-11)==2,"soe2");
  t.test_gen(size_of_exponent(1.0e-1)==2,"soe3");
  t.test_gen(size_of_exponent(1.0e1)==2,"soe4");
  t.test_gen(size_of_exponent(1.0e11)==2,"soe5");
  t.test_gen(size_of_exponent(1.0e111)==3,"soe6");

  t.test_gen(unbackslashify(backslashify("\n"))=="\n","ub 1");
  t.test_gen(unbackslashify(backslashify("\\n"))=="\\n","ub 2");
  t.test_gen(unbackslashify(backslashify("\\\n"))=="\\\n","ub 3");
  t.test_gen(unbackslashify(backslashify("\\\\n"))=="\\\\n","ub 4");
  t.test_gen(unbackslashify(backslashify("\\ n"))=="\\ n","ub 5");
  t.test_gen(unbackslashify(backslashify("\\ \n"))=="\\ \n","ub 6");
  
  vector<string> ss;
  string stemp;

  t.set_output_level(1);

  // Try with out quotes first
  stemp="this is a  test";
  split_string(stemp,ss);
  cout << "\"" << stemp << "\"" << ": ";
  vector_out_quotes(cout,ss,true);
  t.test_gen(ss.size()==4,"ss1");
  t.test_gen(ss[3]==((string)"test"),"ss2");
  ss.clear();

  // Try various cases of quotes at the beginning and end
  stemp="\"this  is\" a test";
  split_string(stemp,ss);
  cout << "\"" << stemp << "\"" << ": ";
  vector_out_quotes(cout,ss,true);
  t.test_gen(ss.size()==3,"ss3");
  t.test_gen(ss[0]==((string)"this is"),"ss4");
  ss.clear();

  stemp="\"  this is \" a test";
  split_string(stemp,ss);
  cout << "\"" << stemp << "\"" << ": ";
  vector_out_quotes(cout,ss,true);
  t.test_gen(ss.size()==3,"ss5");
  //t.test_gen(ss[0]==((string)" this is "),"ss6");
  ss.clear();

  stemp="\"this  is  \" a test";
  split_string(stemp,ss);
  cout << "\"" << stemp << "\"" << ": ";
  vector_out_quotes(cout,ss,true);
  t.test_gen(ss.size()==3,"ss7");
  t.test_gen(ss[0]==((string)"this is "),"ss8");
  ss.clear();

  stemp="\" this is\"  a test";
  split_string(stemp,ss);
  cout << "\"" << stemp << "\"" << ": ";
  vector_out_quotes(cout,ss,true);
  t.test_gen(ss.size()==3,"ss9");
  //t.test_gen(ss[0]==((string)" this is"),"ss10");
  ss.clear();

  // Quotes in the middle
  stemp="th\"is is a test";
  split_string(stemp,ss);
  cout << "\"" << stemp << "\"" << ": ";
  vector_out_quotes(cout,ss,true);
  t.test_gen(ss.size()==4,"ss11");
  t.test_gen(ss[0]==((string)"th\"is"),"ss12");
  ss.clear();

  stemp="th\"is is a test";
  split_string(stemp,ss);
  cout << "\"" << stemp << "\"" << ": ";
  vector_out_quotes(cout,ss,true);
  t.test_gen(ss.size()==4,"ss11");
  t.test_gen(ss[0]==((string)"th\"is"),"ss12");
  ss.clear();

  stemp="\"th\"is is a\" test";
  split_string(stemp,ss);
  cout << "\"" << stemp << "\"" << ": ";
  vector_out_quotes(cout,ss,true);
  t.test_gen(ss.size()==2,"ss13");
  t.test_gen(ss[0]==((string)"th\"is is a"),"ss14");
  ss.clear();

  t.set_output_level(2);
  
  string longstr=((string)"This is a test of a really long string ")+
    "which occupies several lines in a normal 80 column terminal "+
    "window so that I can test where the rewrap() function will "+
    "split the lines.";
  rewrap(longstr,ss);
  t.test_gen(ss.size()==3,"ss5");
  t.test_gen(ss[0]==((string)"This is a test of a really long ")+
	     "string which occupies several lines in a","ss6");
  ss.clear();

  vector<size_t> list;
  o2scl::string_to_uint_list("1-3,7-9,2-5",list);
  t.test_gen(list.size()==10,"list1");
  o2scl::string_to_uint_list("1-3,2",list);
  t.test_gen(list.size()==4,"list2");
  o2scl::string_to_uint_list("1,3-5",list);
  t.test_gen(list.size()==4,"list3");
  o2scl::string_to_uint_list("4",list);
  t.test_gen(list.size()==1,"list4");
  o2scl::string_to_uint_list("4,10",list);
  t.test_gen(list.size()==2,"list5");
  o2scl::string_to_uint_list("4-11,10",list);
  t.test_gen(list.size()==9,"list6");
  o2scl::string_to_uint_list("4,10-11",list);
  t.test_gen(list.size()==3,"list7");

  typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
  
  {
    double pi=boost::math::constants::pi<double>();
    cout << dtos(pi,-1) << endl;
  }
  {
    cpp_dec_float_50 pi=boost::math::constants::pi<cpp_dec_float_50>();
    cout << dtos(pi,-1) << endl;
  }

  cpp_dec_float_50 x=1.5e-200;
  cpp_dec_float_50 mant;
  int expo;
  float_expo_mant(x,expo,mant);
  cout << expo << " " << dtos(mant,10) << endl;

  std::string us1;
  cout << "\nA: " << endl;
  us1=unc_to_string(1.45245234e-4,2.45e-6,1);
  cout << us1 << endl;
  cout << "\nB: " << endl;
  us1=unc_to_string(1.45245234e+6,2.45e+4,1);
  cout << us1 << endl;
  cout << "\nC: " << endl;
  us1=unc_to_string(1.45245234e+4,2.45e+6,1);
  cout << us1 << endl;
  cout << "\nD: " << endl;
  us1=unc_to_string(1.45245234e+4,0.0,1);
  cout << us1 << endl;
  cout << "\nE: " << endl;
  us1=unc_to_string(1.45245234e-1,2.45e-3,1);
  cout << us1 << endl;
  cout << "\nF: " << endl;
  us1=unc_to_string(1.45245234,2.45e-1,1);
  cout << us1 << endl;
  cout << "\nG: " << endl;
  us1=unc_to_string(1.45245234,2.45e-2,1);
  cout << us1 << endl;

  kwargs kw("dtest=2.0,btest=True");
  cout << kw.get_double("dtest") << endl;
  cout << kw.get_bool("btest") << endl;
  
  t.report();
  return 0;
}

