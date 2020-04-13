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
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/convert_units.h>

using namespace o2scl;
using namespace std;

int main(int argc, char *argv[]) {
  cout.setf(ios::scientific);
  test_mgr t;
  t.set_output_level(1);

  // Test all the permutations of mixing two conversions
  // to generate a third

  //1 
  {
    convert_units<double> cu4;
    cu4.insert_cache("in","ft",1.0/12.0);
    cu4.insert_cache("in","cm",2.54);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //2 
  {
    convert_units<double> cu4;
    cu4.insert_cache("ft","in",12.0);
    cu4.insert_cache("in","cm",2.54);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //3 
  {
    convert_units<double> cu4;
    cu4.insert_cache("in","ft",1.0/12.0);
    cu4.insert_cache("cm","in",1.0/2.54);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //4
  {
    convert_units<double> cu4;
    cu4.insert_cache("ft","in",12.0);
    cu4.insert_cache("cm","in",1.0/2.54);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //5
  {
    convert_units<double> cu4;
    cu4.insert_cache("in","cm",2.54);
    cu4.insert_cache("in","ft",1.0/12.0);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //6
  {
    convert_units<double> cu4;
    cu4.insert_cache("in","cm",2.54);
    cu4.insert_cache("ft","in",12.0);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //7
  {
    convert_units<double> cu4;
    cu4.insert_cache("cm","in",1.0/2.54);
    cu4.insert_cache("in","ft",1.0/12.0);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //8
  {
    convert_units<double> cu4;
    cu4.insert_cache("cm","in",1.0/2.54);
    cu4.insert_cache("ft","in",12.0);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  if (argc>=2 && ((string)argv[1])==((string)"test-cache")) {
    // An exhaustive check not intended for the end-user
    convert_units<double> &cu=o2scl_settings.get_convert_units();
    //cu.units_cmd_string=((std::string)"units -f /home/awsteiner")+
    ///wcs/int4/misc/units_hck.dat ";
    cu.test_cache();
  }

  if (argc>=2 && ((string)argv[1])==((string)"make-units")) {
    convert_units<double> cu;
    cu.make_units_dat("units.dat");
  }
  
  if (argc>=2 && ((string)argv[1])==((string)"make-units-hck")) {
    convert_units<double> cu;
    cu.make_units_dat("units_hck.dat",true,true,true);
  }

  {
    convert_units<double> cux;
    double d1, d2;
    int ix;
    cout << endl;
    cout << "Normal conversions:" << endl;
    cout << "nJ kg*m/^2/s^2: ";
    ix=cux.convert_calc("nJ","kg*m^2/s^2",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << endl;
    } else {
      cout << "conversion failed." << endl;
    }
    cout << "kg MeV/c^2: ";
    ix=cux.convert_calc("kg","MeV/c^2",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("kg","MeV/c^2",2.0) << endl;
    cout << "s fm/c: ";
    ix=cux.convert_calc("s","fm/c",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("s","fm/c",2.0) << endl;
    cout << endl;

    // Check that c=1 works
    cout << "Checking c=1:" << endl;
    cout << "kg MeV: ";
    ix=cux.convert_calc("kg","MeV",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << endl;
    } else {
      cout << "conversion failed." << endl;
    }
    cux.set_natural_units(1,0,0);
    
    cout << "kg MeV: ";
    ix=cux.convert_calc("kg","MeV",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("kg","MeV",2.0) << endl;
    
    cout << "kg*m/s MeV: ";
    ix=cux.convert_calc("kg*m/s","MeV",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("kg*m/s","MeV/c",2.0) << endl;

    cout << "MeV/fm^3 Msun/km^3: ";
    ix=cux.convert_calc("MeV/fm^3","Msun/km^3",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("MeV/fm^3","Msun/km^3",2.0) << endl;

    cux.set_natural_units(0,0,0);
    cout << endl;

    // Check that hbar=c=1 works
    cout << "Checking hbar=c=1:" << endl;
    cout << "kg 1/s: ";
    ix=cux.convert_calc("kg","1/s",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << endl;
    } else {
      cout << "conversion failed." << endl;
    }
    
    cux.set_natural_units(1,1,0);
    
    // This requires a factor of hbar/c^2
    cout << "kg 1/s: ";
    ix=cux.convert_calc("kg","1/s",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("kg","hbar/s/c/c",2.0) << endl;
    
    // This requires a factor of hbar/c
    cout << "kg 1/m: ";
    ix=cux.convert_calc("kg","1/m",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("kg","hbar/m/c",2.0) << endl;

    // This requires a factor of hbar 
    cout << "kg s/m^2: ";
    ix=cux.convert_calc("kg","s/m^2",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("kg","hbar*s/m^2",2.0) << endl;

    // This requires a factor of hbar/c^3
    cout << "kg m/s^2: ";
    ix=cux.convert_calc("kg","m/s^2",2.0,d1,d2);
    if (ix==0) {
      cout << "ret,converted: " << ix << " " << d1 << " ";
    } else {
      cout << "conversion failed." << endl;
    }
    cout << o2scl_settings.get_convert_units().convert
      ("kg","hbar*m/s^2/c^3",2.0) << endl;

    cout << endl;
  }
  
  t.report();
  return 0;
}
