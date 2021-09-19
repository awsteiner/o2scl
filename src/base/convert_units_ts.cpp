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

  if (argc>=2 && ((string)argv[1])==((string)"test-cache")) {
    // An exhaustive check not intended for the end-user
    convert_units<double> &cu=o2scl_settings.get_convert_units();
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

  if (true) {
    convert_units<double> cux;
    double d1, d2;
    int ix;

    cux.set_natural_units(1,1,1);
    cux.verbose=2;
    ix=cux.convert_calc("meV","eV",1.0,d1,d2);
    cout << ix << " " << d1 << " " << d2 << endl;
    cux.verbose=0;
    //exit(-1);
    
    // hbar is in kg m^2/s, hbar*c is in kg m^3/s^2

    for(int inu=7;inu>=0;inu--) {
      
      bool c_is_1=(inu&4)>0;
      bool hbar_is_1=(inu&2)>0;
      bool kb_is_1=(inu&1)>0;
      
      cux.set_natural_units(c_is_1,hbar_is_1,kb_is_1);

      cout << "-----------------------------------"
           << "-----------------------------------" << endl;
      cout << "c: " << c_is_1 << " hbar: " << hbar_is_1 << " kb: "
           << kb_is_1 << endl;
      cout << endl;
    
      cout << "J K 1/kB" << endl;
      ix=cux.convert_calc("J","K",2.0,d1,d2);
      if (ix==0) {
        cout << "ix,factor: " << ix << " " << d2 << endl;
        system("acol -get-conv J K");
        t.test_gen(kb_is_1,"kb_is_1 conv 1a");
      } else {
        cout << "ix: " << ix << endl;
        t.test_gen(!kb_is_1,"kb_is_1 conv 1b");
      }
      cout << endl;
    
      cout << "J kg 1/c^2" << endl;
      ix=cux.convert_calc("J","kg",2.0,d1,d2);
      if (ix==0) {
        cout << "ix,factor: " << ix << " " << d2 << endl;
        system("acol -get-conv kg 1/c^2");
        t.test_gen(c_is_1,"c_is_1 conv 2a");
      } else {
        cout << "ix: " << ix << endl;
        t.test_gen(!c_is_1,"c_is_1 conv 2b");
      }
      cout << endl;

      cout << "J^2*s^2 J*s 1/hbar" << endl;
      ix=cux.convert_calc("J^2*s^2","J*s",2.0,d1,d2);
      if (ix==0) {
        cout << "ix,factor: " << ix << " " << d2 << endl;
        system("acol -get-conv J^2*s^2 J*s");
        t.test_gen(hbar_is_1,"hbar_is_1 conv 3a");
      } else {
        cout << "ix: " << ix << endl;
        t.test_gen(!hbar_is_1,"hbar_is_1 conv 3b");
      }
      cout << endl;
    
      cout << "J^2 K^2 1/kB^2" << endl;
      ix=cux.convert_calc("J^2","K^2",2.0,d1,d2);
      if (ix==0) {
        cout << "ix,factor: " << ix << " " << d2 << endl;
        system("acol -get-conv J^2 K^2");
        t.test_gen(kb_is_1,"kb_is_1 conv 4a");
      } else {
        cout << "ix: " << ix << endl;
        t.test_gen(!kb_is_1,"kb_is_1 conv 4b");
      }
      cout << endl;
    
      cout << "erg^2 K^2 1/kB^2" << endl;
      ix=cux.convert_calc("erg^2","K^2",2.0,d1,d2);
      if (ix==0) {
        cout << "ix,factor: " << ix << " " << d2 << endl;
        system("acol -get-conv erg^2 K^2");
        t.test_gen(kb_is_1,"kb_is_1 conv 5a");
      } else {
        cout << "ix: " << ix << endl;
        t.test_gen(!kb_is_1,"kb_is_1 conv 5b");
      }
      cout << endl;
    
      cout << "1/fm^2 K^2 hbar^2*c^2/kB^2" << endl;
      ix=cux.convert_calc("1/fm^2","K^2",2.0,d1,d2);
      if (ix==0) {
        cout << "ix,factor: " << ix << " " << d2 << endl;
        system("acol -get-conv 1/fm^2 K^2");
        t.test_gen(kb_is_1,"kb_is_1 conv 6a");
        t.test_gen(hbar_is_1,"hbar_is_1 conv 6b");
        t.test_gen(c_is_1,"c_is_1 conv 6c");
      } else {
        cout << "ix: " << ix << endl;
        t.test_gen(!c_is_1 || !kb_is_1 || !hbar_is_1,"conv 6d");
      }
      cout << endl;
    
      cout << "1/s^2 nK^2 hbar^2/kB^2" << endl;
      ix=cux.convert_calc("1/s^2","nK^2",2.0,d1,d2);
      if (ix==0) {
        cout << "ix,factor: " << ix << " " << d2 << endl;
        system("acol -get-conv 1/s^2 nK^2");
        t.test_gen(kb_is_1,"kb_is_1 conv 7a");
        t.test_gen(hbar_is_1,"hbar_is_1 conv 7b");
      } else {
        cout << "ix: " << ix << endl;
        t.test_gen(!kb_is_1 || !hbar_is_1,"conv 7c");
      }
      cout << endl;
    
      cout << "m^2 nK^2" << endl;
      ix=cux.convert_calc("m^2","nK^2",2.0,d1,d2);
      if (ix==0) {
        t.test_gen(true,"conv 8a");
      } else {
        cout << "ix: " << ix << endl;
      }
      cout << endl;
    }
    
    cout << "-----------------------------------"
         << "-----------------------------------" << endl;
    cout << "Print units: " << endl;
    cout << endl;
    
    cux.print_units(cout);
    cout << endl;
    
    cout << "-----------------------------------"
         << "-----------------------------------" << endl;
    cout << "Test unique: " << endl;
    cux.test_unique();
    cout << endl;
    
    cout << "-----------------------------------"
         << "-----------------------------------" << endl;
    cout << "Test cache_calc(): " << endl;
    cout << endl;
    
    cux.set_natural_units(1,1,1);
    cux.default_conversions();
    cux.test_cache_calc<test_mgr>(t);
    
  }
  
  t.report();
  return 0;
}
