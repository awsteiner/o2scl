/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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

#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

bool test_mgr::report() const {
  if (output_level>0) {
    if (success) {
      cout << ntests << " tests performed." << endl;
      cout << "All tests passed." << endl;
    } else {
      cout << ntests << " tests performed." << endl;
      cout << "At least one test failed." << endl;
      cout << "Last failed test: " << last_fail << endl;
    }
  }
  return success;
}

void test_mgr::process_test(bool ret, string d2, string description) {
  ntests++;
  if (!ret) {
    last_fail=description;
    success=false;
  }
  if (output_level>1 || (output_level>0 && ret==false)) {
    if (ret==true) {
      cout << "PASS: " << d2 << " " << description << endl;
    } else {
      cout << "FAIL: " << d2 << " " << description << endl;
    }
  }
  return;
}

bool test_mgr::test_gen(bool value, std::string description) {
  
  process_test(value,"general",description);

  return value;
}

bool test_mgr::test_str(std::string result, std::string expected, 
			std::string description) {
  bool ret;

  ret=(result==expected);

  description=result+" vs. "+expected+"\n "+
    description;
  process_test(ret,"string",description);

  return ret;
}

const test_mgr operator+(const test_mgr &left,
			 const test_mgr &right) {

  bool success=(left.get_success() && right.get_success());
  string lf;
  if (left.get_success()==false) {
    lf="operator+() "+left.get_last_fail();
  } 
  if (right.get_success()==false) {
    if (lf.length()==0) lf="operator+()";
    lf=" "+right.get_last_fail();
  }
  int ntests=left.get_ntests()+right.get_ntests();
  int output_level=left.get_output_level();
  if (right.get_output_level()>output_level) {
    output_level=right.get_output_level();
  }

  return test_mgr(success,lf,ntests,output_level);
}

