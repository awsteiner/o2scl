/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/nucdist.h>

using namespace std;
using namespace o2scl;

void o2scl::nucdist_set(vector<nucleus> &dist, nucmass &nm, 
			std::string expr, int maxA,
			bool include_neutron) {
  
  nucleus n;

  if (dist.size()>0) dist.clear();

  /// The function parser
  calc_utf8<> calc;
  std::map<std::string,double> vars;
  calc.compile(expr.c_str(),&vars);
  
  size_t dist_size=0;
  
  // Now fill the vector with the nuclei
  size_t ix=0;
  for(int A=1;A<=maxA;A++) {
    for(int Z=0;Z<=A;Z++) {
      if (A==1 && Z==0) {
	if (include_neutron && nm.is_included(0,1)) {
	  dist.push_back(n);
	  nm.get_nucleus(Z,A-Z,dist[ix]);
	  ix++;
	}
      } else {
	vars["Z"]=Z;
	vars["A"]=A;
	vars["N"]=A-Z;
	if (nm.is_included(Z,A-Z) && calc.eval(&vars)>0.5) {
	  dist.push_back(n);
	  nm.get_nucleus(Z,A-Z,dist[ix]);
	  ix++;
	}
      }
    }
  }

  return;
}

