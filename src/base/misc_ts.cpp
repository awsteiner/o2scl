/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

int main(void) {
  double x=sqrt(2.0);
  int i=3;
  string s1="4", s2="1.7320508";
  test_mgr t;
  t.set_output_level(2);
  
  string s[10]={"test1","test_of_string2","test_of_string3",
		"test_of_string4","test5","test_of_string6",
		"test_of_string7","test_of_string8","test_of_string9",
		"test_of_string10"};
  vector<string> sout;
  screenify(10,s,sout);
  for(size_t ij=0;ij<sout.size();ij++) {
    cout << sout[ij] << endl;
  }
  t.test_gen(sout.size()==3,"screenify");

  t.test_gen(o2scl::count_words(" This is a hard test ")==5,"count1");
  t.test_gen(o2scl::count_words("This is an easy test")==5,"count2");
  t.test_gen(o2scl::count_words("This is  \t a\t test  \twith\ttabs")==6,
	     "count3");

  string word_test="x x";
  cout << "Whitespace characters: " << flush;
  for(i=0;i<256;i++) {
    word_test[1]=((char)i);
    if (count_words(word_test)==2) {
      cout << i << " " << flush;
    }
  }
  cout << endl;
  
  cout << fermi_function(39.9,0.0,1.0) << " " 
       << fermi_function(40.1,0.0,1.0) << " " 
       << 1.0/(1.0+exp(40.1)) << endl;

  // This generates a list of numbers for the gen_test_number
  // class documentation
  cout.setf(ios::scientific);
  gen_test_number<10> gn;
  for(size_t ij=0;ij<15;ij++) {
    cout.width(2);
    cout.setf(ios::left);
    cout << ij << " ";
    cout.setf(ios::showpos);
    cout << gn.gen() << endl;
    cout.unsetf(ios::showpos);
  }

  string ss=" This is a\n test ";
  remove_whitespace(ss);
  t.test_gen(ss=="Thisisatest","remove_whitespace()");
  ss="\tThis is a\n test\n";
  remove_whitespace(ss);
  t.test_gen(ss=="Thisisatest","remove_whitespace()");

#ifndef O2SCL_OLDER_COMPILER

  // vec_index doesn't work in older compilers, so
  // we comment this out for now.
  
  vec_index vi;
  vi.append("a1");
  vi.append("c1");
  vi.append("b1");
  t.test_gen(vi(1)==((string)"c1"),"vec_index 1");
  t.test_gen(vi[2]==((string)"b1"),"vec_index 2");
  t.test_gen(vi("a1")==0,"vec_index 3");
  t.test_gen(vi("b1")==2,"vec_index 4");
  vec_index vi2={"a1","c1","b1"};
  t.test_gen(vi2(1)==((string)"c1"),"vec_index 5");
  t.test_gen(vi2[2]==((string)"b1"),"vec_index 6");
  t.test_gen(vi2("a1")==0,"vec_index 7");
  t.test_gen(vi2("b1")==2,"vec_index 8");
  
#endif

#ifdef O2SCL_NEVER_DEFINED

  // This appears not to be sufficiently platform independent,
  // see 
  // vec_index doesn't work in older compilers, so
  // we comment this out for now.
  
  vec_index vi;
  vi.append("a1");
  vi.append("c1");
  vi.append("b1");
  t.test_gen(vi(1)==((string)"c1"),"vec_index 1");
  t.test_gen(vi[2]==((string)"b1"),"vec_index 2");
  t.test_gen(vi("a1")==0,"vec_index 3");
  t.test_gen(vi("b1")==2,"vec_index 4");
  vec_index vi2={"a1","c1","b1"};
  t.test_gen(vi2(1)==((string)"c1"),"vec_index 5");
  t.test_gen(vi2[2]==((string)"b1"),"vec_index 6");
  t.test_gen(vi2("a1")==0,"vec_index 7");
  t.test_gen(vi2("b1")==2,"vec_index 8");
  
#endif

#ifdef O2SCL_NEVER_DEFINED

  // This appears not to be sufficiently platform independent, see
  // https://github.com/awsteiner/o2scl/issues/8 for problem with
  // openSUSE.
  vector<std::string> matches;
  glob_wrapper("../anneal/*.h",matches);
  t.test_gen(matches.size()==4,"glob test");
  
#endif
  
  t.report();
  return 0;
}

