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
#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/convert_units.h>

using namespace o2scl;
using namespace std;

int main(int argc, char *argv[]) {
  
  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(1);

  convert_units<double> cux;
  double d1, d2;
  int ix;

  cux.set_natural_units(1,1,1);
  cux.verbose=2;
  ix=cux.convert_ret("meV","eV",1.0,d1);
  cout << ix << " " << d1 << endl;
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

    int sret=0;

    // AWS, 10/6/21: this testing code is nice, but it uses
    // acol, and acol isn't available via system() until
    // LD_LIBRARY_FLAGS is set, so I need to rewrite these
    // tests. FIXME. Additionally, convert_calc() is now
    // protected.

#ifdef O2SCL_NEVER_DEFINED
    if (false) {
      
      cout << "J K 1/kB" << endl;
      ix=cux.convert_calc("J","K",2.0,d1,d2);
      if (ix==0) {
        cout << "ix,factor: " << ix << " " << d2 << endl;
        sret=system("acol -convert J K");
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
        sret=system("acol -convert kg 1/c^2");
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
        sret=system("acol -convert J^2*s^2 J*s");
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
        sret=system("acol -convert J^2 K^2");
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
        sret=system("acol -convert erg^2 K^2");
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
        sret=system("acol -convert 1/fm^2 K^2");
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
        sret=system("acol -convert 1/s^2 nK^2");
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
    
#endif
    
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
  cout << endl;
  
  cux.verbose=1;
  int cret=cux.convert_ret("g","solarmass",1.0,d1);
  t.test_rel(1.0/o2scl_const::solar_mass_f<double>
             (o2scl_const::o2scl_cgs),d1,1.0e-6,"calc2");

  cux.err_on_fail=false;
  int iret=cux.convert_ret("α","N/K",3.0,d1);
  t.test_gen(iret!=0,"convert with new unit 0");
  
  // With these values, alpha is basically 3 Newtons per Kelvin
  convert_units<>::der_unit d;
  d.label="α";
  d.val=3.0;
  d.name="alpha unit";
  d.m=1;
  d.k=1;
  d.s=-2;
  d.K=-1;
  d.A=0;
  d.mol=0;
  d.cd=0;
  cux.add_unit(d);

  // Convert from 3 alpha to N/K to get 9. FIXME. This doesn't
  // work yet.
  iret=cux.convert_ret("α","N/K",3.0,d1);
  t.test_gen(iret==0,"convert with new unit 1");
  t.test_rel(d1,9.0,1.0e-15,"convert with new unit 2");

  t.test_gen(cux.is_in_cache("α","N/K")==1,"in_cache() 1");
  t.test_gen(cux.is_in_cache("N/K","α")==2,"in_cache() 2");
  
  cux.del_unit("α");

  // This succeeds because the unit conversion is still in the
  // cache
  iret=cux.convert_ret("α","N/K",3.0,d1);
  t.test_gen(iret==0,"convert with new unit 3");

  // Now additionally remove it from the cache
  cux.remove_cache("α","N/K");

  // This now fails
  iret=cux.convert_ret("α","N/K",3.0,d1);
  t.test_gen(iret!=0,"convert with new unit 4");
  
  t.report();
  return 0;
}
