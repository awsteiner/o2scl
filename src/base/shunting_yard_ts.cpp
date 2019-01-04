/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2019, Andrew W. Steiner
  
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
#include <o2scl/shunting_yard.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  calculator calc;
  calc.compile("1.0e-100<2.0e-100",0);
  t.test_gen(calc.eval(0)==1,"calc1");
  calc.compile("-(10)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-10.0,"calc2");
  calc.compile("((10+1)*2+3)*2+1",0);
  cout << calc.RPN_to_string() << endl;

  calc.compile("(3 && true) == true",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==true,"calc3");
  calc.compile("(3 && 0) == true",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==false,"calc4");
  calc.compile("(3 || 0) == true",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==true,"calc5");
  calc.compile("(false || 0) == true",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==false,"calc6");

  // Test new functions
  calc.compile("sin(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==sin(0.2),"calc7");
  calc.compile("-sin(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-sin(0.2),"calc8");
  calc.compile("1-sin(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==1-sin(0.2),"calc9");
  calc.compile("1+sin(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==1+sin(0.2),"calc10");
  calc.compile("-(-sin(0.2))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==sin(0.2),"calc11");
  calc.compile("-sin(0.2+3*(4+5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-sin(0.2+3*(4+5)),"calc12");

  // Test new functions
  calc.compile("sqrt(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==sqrt(0.2),"calc13");
  calc.compile("-sqrt(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-sqrt(0.2),"calc14");
  calc.compile("1-sqrt(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==1-sqrt(0.2),"calc15");
  calc.compile("1+sqrt(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==1+sqrt(0.2),"calc16");
  calc.compile("-(-sqrt(0.2))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==sqrt(0.2),"calc17");
  calc.compile("-sqrt(0.2+3*(4+5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-sqrt(0.2+3*(4+5)),"calc18");

  // Test new functions
  calc.compile("exp(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==exp(0.2),"calc19");
  calc.compile("-exp(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-exp(0.2),"calc20");
  calc.compile("1-exp(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==1-exp(0.2),"calc21");
  calc.compile("1+exp(0.2)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==1+exp(0.2),"calc22");
  calc.compile("-(-exp(0.2))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==exp(0.2),"calc23");
  calc.compile("-exp(0.2+3*(4+5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-exp(0.2+3*(4+5)),"calc24");

  // Mixing functions
  calc.compile("-exp(0.2+sin(4+5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-exp(0.2+sin(4+5)),"calc25");
  calc.compile("-exp(0.2)+sin(4)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==-exp(0.2)+sin(4),"calc26");
  calc.compile("exp(0.2)+sin(4)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==exp(0.2)+sin(4),"calc27");
  calc.compile("exp(0.2)^2*sin(4)",0);
  cout << calc.RPN_to_string() << endl;
  t.test_gen(calc.eval(0)==pow(exp(0.2),2.0)*sin(4),"calc28");

  // Test new functions
  calc.compile("asin(sin(0.5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_rel(calc.eval(0),0.5,1.0e-14,"calc29");
  calc.compile("acos(cos(0.5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_rel(calc.eval(0),0.5,1.0e-14,"calc30");
  calc.compile("atan(tan(0.5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_rel(calc.eval(0),0.5,1.0e-14,"calc31");
  calc.compile("asinh(sinh(0.5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_rel(calc.eval(0),0.5,1.0e-14,"calc32");
  calc.compile("acosh(cosh(0.5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_rel(calc.eval(0),0.5,1.0e-14,"calc33");
  calc.compile("atanh(tanh(0.5))",0);
  cout << calc.RPN_to_string() << endl;
  t.test_rel(calc.eval(0),0.5,1.0e-14,"calc34");

  t.report();
  return 0;
}

