/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2022, Andrew W. Steiner

  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------
*/
#include "format_float.h"
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(int argc, char *argv[]) {
  test_mgr t;
  t.set_output_level(2);

  format_float ff;

  ff.latex_mode();
  t.test_gen(ff.convert(-sqrt(2.0)*1.0e-5)=="$-$1.4142 $\\times 10^{-5}$",
	     "test 1");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-4)=="1.4142 $\\times 10^{-4}$",
	     "test 2");
  t.test_gen(ff.convert(-sqrt(2.0)*1.0e-3)=="$-$1.4142 $\\times 10^{-3}$",
	     "test 3");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-2)=="0.014142",
	     "test 4");
  t.test_gen(ff.convert(-sqrt(2.0)*1.0e-1)=="$-$0.14142",
	     "test 5");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e0)=="1.4142",
	     "test 6");
  t.test_gen(ff.convert(-sqrt(2.0)*1.0e1)=="$-$14.142",
	     "test 7");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e2)=="141.42",
	     "test 8");
  t.test_gen(ff.convert(-sqrt(2.0)*1.0e3)=="$-$1414.2",
	     "test 9");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e4)=="1.4142 $\\times 10^{4}$",
	     "test 10");
  t.test_gen(ff.convert(-sqrt(2.0)*1.0e5)=="$-$1.4142 $\\times 10^{5}$",
	     "test 11");

  ff.set_sig_figs(3);
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-5)=="1.41 $\\times 10^{-5}$",
	     "test 12");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-2)=="0.0141","test 13");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e0)=="1.41","test 14");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e2)=="141","test 15");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e5)=="1.41 $\\times 10^{5}$","test 16");
 
  ff.html_mode();
  ff.set_sig_figs(5);
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-5)=="1.4142 &times; 10<sup>&mdash;5</sup>",
	     "test 17");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-4)=="1.4142 &times; 10<sup>&mdash;4</sup>",
	     "test 18");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-3)=="1.4142 &times; 10<sup>&mdash;3</sup>",
	     "test 19");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-2)=="0.014142",
	     "test 20");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e-1)=="0.14142",
	     "test 21");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e0)=="1.4142",
	     "test 22");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e1)=="14.142",
	     "test 23");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e2)=="141.42",
	     "test 24");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e3)=="1414.2",
	     "test 25");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e4)=="1.4142 &times; 10<sup>4</sup>",
	     "test 26");
  t.test_gen(ff.convert(sqrt(2.0)*1.0e5)=="1.4142 &times; 10<sup>5</sup>",
	     "test 27");

  ff.set_sig_figs(7);
  ff.set_pad_zeros(true);
  double num=1.414;
  t.test_gen(ff.convert(num*1.0e-5)=="1.414000 &times; 10<sup>&mdash;5</sup>",
	     "test 28");
  t.test_gen(ff.convert(num*1.0e-4)=="1.414000 &times; 10<sup>&mdash;4</sup>",
	     "test 29");
  t.test_gen(ff.convert(num*1.0e-3)=="1.414000 &times; 10<sup>&mdash;3</sup>",
	     "test 30");
  t.test_gen(ff.convert(num*1.0e-2)=="0.01414000",
	     "test 31");
  t.test_gen(ff.convert(num*1.0e-1)=="0.1414000",
	     "test 32");
  t.test_gen(ff.convert(num*1.0e0)=="1.414000",
	     "test 33");
  t.test_gen(ff.convert(num*1.0e1)=="14.14000",
	     "test 34");
  t.test_gen(ff.convert(num*1.0e2)=="141.4000",
	     "test 35");
  t.test_gen(ff.convert(num*1.0e3)=="1414.000",
	     "test 36");
  t.test_gen(ff.convert(num*1.0e4)=="1.414000 &times; 10<sup>4</sup>",
	     "test 37");
  t.test_gen(ff.convert(num*1.0e5)=="1.414000 &times; 10<sup>5</sup>",
	     "test 38");

  ff.set_sig_figs(3);
  ff.set_pad_zeros(true);
  t.test_gen(ff.convert(num*1.0e-5)=="1.41 &times; 10<sup>&mdash;5</sup>",
	     "test 39");
  t.test_gen(ff.convert(num*1.0e-4)=="1.41 &times; 10<sup>&mdash;4</sup>",
	     "test 40");
  t.test_gen(ff.convert(num*1.0e-3)=="1.41 &times; 10<sup>&mdash;3</sup>",
	     "test 41");
  t.test_gen(ff.convert(num*1.0e-2)=="0.0141",
	     "test 42");
  t.test_gen(ff.convert(num*1.0e-1)=="0.141",
	     "test 43");
  t.test_gen(ff.convert(num*1.0e0)=="1.41",
	     "test 44");
  t.test_gen(ff.convert(num*1.0e1)=="14.1",
	     "test 45");
  t.test_gen(ff.convert(num*1.0e2)=="141",
	     "test 46");
  t.test_gen(ff.convert(num*1.0e3)=="1410",
	     "test 47");
  t.test_gen(ff.convert(num*1.0e4)=="1.41 &times; 10<sup>4</sup>",
	     "test 48");
  t.test_gen(ff.convert(num*1.0e5)=="1.41 &times; 10<sup>5</sup>",
	     "test 49");

  num=1.4;
  t.test_gen(ff.convert(num*1.0e-5)=="1.40 &times; 10<sup>&mdash;5</sup>",
	     "test 50");
  t.test_gen(ff.convert(num*1.0e-4)=="1.40 &times; 10<sup>&mdash;4</sup>",
	     "test 51");
  t.test_gen(ff.convert(num*1.0e-3)=="1.40 &times; 10<sup>&mdash;3</sup>",
	     "test 52");
  t.test_gen(ff.convert(num*1.0e-2)=="0.0140",
	     "test 53");
  t.test_gen(ff.convert(num*1.0e-1)=="0.140",
	     "test 54");
  t.test_gen(ff.convert(num*1.0e0)=="1.40",
	     "test 55");
  t.test_gen(ff.convert(num*1.0e1)=="14.0",
	     "test 56");
  t.test_gen(ff.convert(num*1.0e2)=="140",
	     "test 57");
  t.test_gen(ff.convert(num*1.0e3)=="1400",
	     "test 58");
  t.test_gen(ff.convert(num*1.0e4)=="1.40 &times; 10<sup>4</sup>",
	     "test 59");
  t.test_gen(ff.convert(num*1.0e5)=="1.40 &times; 10<sup>5</sup>",
	     "test 60");

  ff.latex_mode();

  ff.set_pad_zeros(true);

  ff.set_sig_figs(3);
  t.test_gen(ff.convert(0.0)=="0","test 61");
  t.test_gen(ff.convert(0.001)=="1.00 $\\times 10^{-3}$","test 62");
  t.test_gen(ff.convert(0.002)=="2.00 $\\times 10^{-3}$","test 63");
  t.test_gen(ff.convert(0.003)=="3.00 $\\times 10^{-3}$","test 64");
  t.test_gen(ff.convert(0.004)=="4.00 $\\times 10^{-3}$","test 65");
  t.test_gen(ff.convert(0.005)=="5.00 $\\times 10^{-3}$","test 66");
  t.test_gen(ff.convert(0.006)=="6.00 $\\times 10^{-3}$","test 67");
  ff.set_sig_figs(3);
  t.test_gen(ff.convert(1.0)=="1.00","test 68");
  t.test_gen(ff.convert(1.001)=="1.00","test 69");
  t.test_gen(ff.convert(1.002)=="1.00","test 70");
  t.test_gen(ff.convert(1.003)=="1.00","test 71");
  t.test_gen(ff.convert(1.004)=="1.00","test 72");
  t.test_gen(ff.convert(1.005)=="1.00","test 73");
  t.test_gen(ff.convert(1.006)=="1.01","test 74");
  ff.set_sig_figs(2);
  t.test_gen(ff.convert(1.0)=="1.0","test 75");
  t.test_gen(ff.convert(1.01)=="1.0","test 76");
  t.test_gen(ff.convert(1.02)=="1.0","test 77");
  t.test_gen(ff.convert(1.03)=="1.0","test 78");
  t.test_gen(ff.convert(1.04)=="1.0","test 79");
  t.test_gen(ff.convert(1.05)=="1.1","test 80");
  t.test_gen(ff.convert(1.06)=="1.1","test 81");
  ff.set_sig_figs(1);
  t.test_gen(ff.convert(1.0)=="1","test 82");
  t.test_gen(ff.convert(1.0)=="1","test 83");
  t.test_gen(ff.convert(1.0)=="1","test 84");
  t.test_gen(ff.convert(1.1)=="1","test 85");
  t.test_gen(ff.convert(1.2)=="1","test 86");
  t.test_gen(ff.convert(1.3)=="1","test 87");
  t.test_gen(ff.convert(1.4)=="1","test 88");
  t.test_gen(ff.convert(1.5)=="2","test 89");
  t.test_gen(ff.convert(1.6)=="2","test 90");
  t.test_gen(ff.convert(1.0)=="1","test 91");
  t.test_gen(ff.convert(10.0)=="10","test 92");
  t.test_gen(ff.convert(100.0)=="100","test 93");
  t.test_gen(ff.convert(1000.0)=="1000","test 94");
  t.test_gen(ff.convert(10000.0)=="1 $\\times 10^{4}$","test 95");
  t.test_gen(ff.convert(0.1)=="0.1","test 96");
  t.test_gen(ff.convert(0.01)=="0.01","test 97");
  t.test_gen(ff.convert(0.001)=="1 $\\times 10^{-3}$","test 98");
  t.test_gen(ff.convert(0.0001)=="1 $\\times 10^{-4}$","test 99");
  ff.set_sig_figs(4);
  t.test_gen(ff.convert(1.23)=="1.230","test 100");
  t.test_gen(ff.convert(12.3)=="12.30","test 101");
  t.test_gen(ff.convert(123.0)=="123.0","test 102");
  t.test_gen(ff.convert(1230.0)=="1230","test 103");
  t.test_gen(ff.convert(12300.0)=="1.230 $\\times 10^{4}$","test 104");
  t.test_gen(ff.convert(0.123)=="0.1230","test 105");
  t.test_gen(ff.convert(0.0123)=="0.01230","test 106");
  t.test_gen(ff.convert(0.00123)=="1.230 $\\times 10^{-3}$","test 107");
  t.test_gen(ff.convert(0.000123)=="1.230 $\\times 10^{-4}$","test 108");

  ff.set_pad_zeros(false);

  ff.set_sig_figs(3);
  t.test_gen(ff.convert(0.0)=="0","test 61b");
  t.test_gen(ff.convert(0.001)=="1 $\\times 10^{-3}$","test 62b");
  t.test_gen(ff.convert(0.002)=="2 $\\times 10^{-3}$","test 63b");
  t.test_gen(ff.convert(0.003)=="3 $\\times 10^{-3}$","test 64b");
  t.test_gen(ff.convert(0.004)=="4 $\\times 10^{-3}$","test 65b");
  t.test_gen(ff.convert(0.005)=="5 $\\times 10^{-3}$","test 66b");
  t.test_gen(ff.convert(0.006)=="6 $\\times 10^{-3}$","test 67b");
  ff.set_sig_figs(3);
  t.test_gen(ff.convert(1.0)=="1","test 68b");
  t.test_gen(ff.convert(1.001)=="1","test 69b");
  t.test_gen(ff.convert(1.002)=="1","test 70b");
  t.test_gen(ff.convert(1.003)=="1","test 71b");
  t.test_gen(ff.convert(1.004)=="1","test 72b");
  t.test_gen(ff.convert(1.005)=="1","test 73b");
  t.test_gen(ff.convert(1.006)=="1.01","test 74b");
  ff.set_sig_figs(2);
  t.test_gen(ff.convert(1.0)=="1","test 75b");
  t.test_gen(ff.convert(1.01)=="1","test 76b");
  t.test_gen(ff.convert(1.02)=="1","test 77b");
  t.test_gen(ff.convert(1.03)=="1","test 78b");
  t.test_gen(ff.convert(1.04)=="1","test 79b");
  t.test_gen(ff.convert(1.05)=="1.1","test 80b");
  t.test_gen(ff.convert(1.06)=="1.1","test 81b");
  ff.set_sig_figs(1);
  t.test_gen(ff.convert(1.0)=="1","test 82b");
  t.test_gen(ff.convert(1.0)=="1","test 83b");
  t.test_gen(ff.convert(1.0)=="1","test 84b");
  t.test_gen(ff.convert(1.1)=="1","test 85b");
  t.test_gen(ff.convert(1.2)=="1","test 86b");
  t.test_gen(ff.convert(1.3)=="1","test 87b");
  t.test_gen(ff.convert(1.4)=="1","test 88b");
  t.test_gen(ff.convert(1.5)=="2","test 89b");
  t.test_gen(ff.convert(1.6)=="2","test 90b");
  t.test_gen(ff.convert(1.0)=="1","test 91b");
  t.test_gen(ff.convert(10.0)=="10","test 92b");
  t.test_gen(ff.convert(100.0)=="100","test 93b");
  t.test_gen(ff.convert(1000.0)=="1000","test 94b");
  t.test_gen(ff.convert(10000.0)=="1 $\\times 10^{4}$","test 95b");
  t.test_gen(ff.convert(0.1)=="0.1","test 96b");
  t.test_gen(ff.convert(0.01)=="0.01","test 97b");
  t.test_gen(ff.convert(0.001)=="1 $\\times 10^{-3}$","test 98b");
  t.test_gen(ff.convert(0.0001)=="1 $\\times 10^{-4}$","test 99b");
  ff.set_sig_figs(4);
  t.test_gen(ff.convert(1.23)=="1.23","test 100b");
  t.test_gen(ff.convert(12.3)=="12.3","test 101b");
  t.test_gen(ff.convert(123.0)=="123","test 102b");
  t.test_gen(ff.convert(1230.0)=="1230","test 103b");
  t.test_gen(ff.convert(12300.0)=="1.23 $\\times 10^{4}$","test 104b");
  t.test_gen(ff.convert(0.123)=="0.123","test 105b");
  t.test_gen(ff.convert(0.0123)=="0.0123","test 106b");
  t.test_gen(ff.convert(0.00123)=="1.23 $\\times 10^{-3}$","test 107b");
  t.test_gen(ff.convert(0.000123)=="1.23 $\\times 10^{-4}$","test 108b");

  // Test the algorithm by comparing the results in c_mode()
  // to the output from ostringstream

  ff.c_mode();
  ff.set_sig_figs(6);
  static const size_t N=20;
  gen_test_number<> gtn;
  for(size_t i=0;i<N;i++) {

    double x=gtn.gen();
    {
      ostringstream s1;
      ostringstream s2;
      s1 << x;
      s2 << ff.convert(x);
      t.test_gen(s1.str()==s2.str(),"test 1c");
    }
    double x2=x*1.0e3;
    {
      ostringstream s1;
      ostringstream s2;
      s1 << x2;
      s2 << ff.convert(x2);
      t.test_gen(s1.str()==s2.str(),"test 2c");
    }
    double x3=x*1.0e6;
    {
      ostringstream s1;
      ostringstream s2;
      s1 << x3;
      s2 << ff.convert(x3);
      t.test_gen(s1.str()==s2.str(),"test 3c");
    }
    double x4=x*1.0e-3;
    {
      ostringstream s1;
      ostringstream s2;
      s1 << x4;
      s2 << ff.convert(x4);
      t.test_gen(s1.str()==s2.str(),"test 4c");
    }
    double x5=x*1.0e-6;
    {
      ostringstream s1;
      ostringstream s2;
      s1 << x5;
      s2 << ff.convert(x5);
      t.test_gen(s1.str()==s2.str(),"test 5c");
    }
  }

  cout << "Unicode: " << endl;
  ff.unicode_mode();
  cout << ff.convert(sqrt(2.0)*1.0e-10) << endl;
  cout << ff.convert(sqrt(2.0)*1.0e-5) << endl;
  cout << ff.convert(sqrt(2.0)*1.0e7) << endl;
  cout << ff.convert(sqrt(2.0)*1.0e14) << endl;
  
  t.report();
  return 0;
}

