/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

int main(void) {
  cout.setf(ios::scientific);
  test_mgr t;
  t.set_output_level(1);

  // Test all the permutations of mixing two conversions
  // to generate a third

  //1 
  {
    convert_units cu4;
    cu4.insert_cache("in","ft",1.0/12.0);
    cu4.insert_cache("in","cm",2.54);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //2 
  {
    convert_units cu4;
    cu4.insert_cache("ft","in",12.0);
    cu4.insert_cache("in","cm",2.54);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //3 
  {
    convert_units cu4;
    cu4.insert_cache("in","ft",1.0/12.0);
    cu4.insert_cache("cm","in",1.0/2.54);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //4
  {
    convert_units cu4;
    cu4.insert_cache("ft","in",12.0);
    cu4.insert_cache("cm","in",1.0/2.54);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //5
  {
    convert_units cu4;
    cu4.insert_cache("in","cm",2.54);
    cu4.insert_cache("in","ft",1.0/12.0);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //6
  {
    convert_units cu4;
    cu4.insert_cache("in","cm",2.54);
    cu4.insert_cache("ft","in",12.0);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //7
  {
    convert_units cu4;
    cu4.insert_cache("cm","in",1.0/2.54);
    cu4.insert_cache("in","ft",1.0/12.0);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  //8
  {
    convert_units cu4;
    cu4.insert_cache("cm","in",1.0/2.54);
    cu4.insert_cache("ft","in",12.0);

    double res=cu4.convert("ft","cm",1.0);
    t.test_rel(res,12.0*2.540,1.0e-10,"1");
    res=cu4.convert("cm","ft",1.0);
    t.test_rel(res,1.0/(12.0*2.540),1.0e-10,"1");
  }

  // This is an exhaustive check not for the end-user
  // and is thus commented out
  if (false) {
    convert_units &cu=o2scl_settings.get_convert_units();
    //cu.units_cmd_string=((std::string)"units -f /home/awsteiner")+
    ///wcs/int4/misc/units_hck.dat ";
    cu.test_cache();
  }
  
  t.report();
  return 0;
}
