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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/funct_to_fp.h>
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

int o2scl::function_to_double_nothrow(std::string s, double &result,
                                      int verbose, rng<> *r) {
  convert_units<double> &cu=o2scl_settings.get_convert_units();
  return function_to_fp_nothrow<double>
    (s,result,cu,verbose,r);
}

double o2scl::function_to_double(std::string s, int verbose) {
  double res;
  int ret=function_to_double_nothrow(s,res,verbose);
  if (ret!=0) {
    O2SCL_ERR("Function function_to_double() failed.",ret);
  }
  return res;
}
