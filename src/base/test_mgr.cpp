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

#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

bool test_mgr::report() {
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

bool test_mgr::test_rel(double result, double expected, double rel_error,
			std::string description) {
  bool ret;
  if (std::isnan(expected)) {
    ret=(std::isnan(expected)==std::isnan(result));
  } else if (std::isinf(expected)) {
    ret=(std::isinf(expected)==std::isinf(result));
  } else if (expected==0.0) {
    ret=test_abs(result,expected,rel_error,description);
    return ret;
  } else {
    ret=((fabs(expected-result))/fabs(expected)<rel_error);
  }
  
  description=dtos(result)+" vs. "+dtos(expected)+
    " is "+dtos(fabs(expected-result)/fabs(expected))+"\n "+
    description;
  process_test(ret,"relative",description);

  return ret;
}

bool test_mgr::test_abs(double result, double expected, double abs_error,
			std::string description) {
  bool ret;
  if (std::isnan(expected)) {
    ret=(std::isnan(expected)==std::isnan(result));
    description=dtos(result)+" vs. "+ dtos(expected)+
      "\n "+description;
  } else if (std::isinf(expected)) {
    ret=(std::isinf(expected)==std::isinf(result));
    description=dtos(result)+" vs. "+ dtos(expected)+
      "\n "+description;
  } else {
    ret=(fabs(expected-result)<abs_error);
    description=dtos(result)+" vs. "+ dtos(expected)+" is "
      +dtos(fabs(expected-result))+"\n "+description;
  }
  
  process_test(ret,"absolute",description);

  return ret;
}

bool test_mgr::test_fact(double result, double expected, double factor,
			 std::string description) {
  bool ret;
  double ratio;
  if (std::isnan(expected)) {
    ret=(std::isnan(expected)==std::isnan(result));
  } else if (std::isinf(expected)) {
    ret=(std::isinf(expected)==std::isinf(result));
  } else {
    ratio=expected/result;
    ret=(ratio<factor && ratio>1.0/factor);
  }

  description= dtos(result)+" vs. "+ dtos(expected)+"\n "+
    description;
  process_test(ret,"factor",description);

  return ret;
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

const test_mgr operator+(const test_mgr& left,
			 const test_mgr& right) {
  test_mgr ret;
  ret.success=(left.success && right.success);
  if (left.success==false) {
    ret.last_fail=left.last_fail;
  } else if (right.success==false) {
    ret.last_fail=right.last_fail;
  } else {
    ret.last_fail="";
  }
  return ret;
}

