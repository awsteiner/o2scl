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
#include <o2scl/nucdist.h>
#include <o2scl/fparser.h>

using namespace std;
using namespace o2scl;

void o2scl::nucdist_set(vector<nucleus> &dist, nucmass &nm, 
			std::string expr, int maxA,
			bool include_neutron) {
  
  nucleus n;

  if (dist.size()>0) dist.clear();

  /// The function parser
  FunctionParser fp;
  double vals[2];
  
  // Parse the formula
  int ret=fp.Parse(expr,"Z,N");
  if (ret!=-1) {
    O2SCL_ERR("Failed to parse in nucdist_set().",exc_einval);
  }
  
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
	vals[0]=Z;
	vals[1]=A-Z;
	if (nm.is_included(Z,A-Z) && fp.Eval(vals)) {
	  dist.push_back(n);
	  nm.get_nucleus(Z,A-Z,dist[ix]);
	  ix++;
	}
      }
    }
  }

  return;
}

bool o2scl::operator==(const nucdist::iterator &i1,
		       const nucdist::iterator &i2) {
  return i1.np==i2.np;
}

bool o2scl::operator!=(const nucdist::iterator &i1,
		       const nucdist::iterator &i2) {
  return i1.np!=i2.np;
}

nucdist_full::nucdist_full(nucmass &nm, int maxA, bool include_neutron) {
  list_size=0;
  set_dist(nm,maxA,include_neutron);
}

void nucdist_full::set_dist(nucmass &nm, int maxA, bool include_neutron) {
  
  // Delete previous list if necessary
  if (list_size>0) delete[] list;

  nucleus n;

  // First pass, count number of nuclei
  list_size=0;
  // Count neutron
  if (include_neutron && nm.is_included(0,1)) {
    list_size++;
  }
  // Count remaining nuclei
  for(int A=1;A<=maxA;A++) {
    for(int Z=0;Z<=A;Z++) {
      if (A==1 && Z==0) {
	if (include_neutron && nm.is_included(0,1)) {
	  list_size++;
	}
      } else {
	if (nm.is_included(Z,A-Z)) {
	  list_size++;
	}
      }
    }
  }
  // Allocate the vector
  list=new nucleus[list_size];

  // Now fill the vector with the nuclei
  size_t ix=0;
  for(int A=1;A<=maxA;A++) {
    for(int Z=0;Z<=A;Z++) {
      if (A==1 && Z==0) {
	if (include_neutron && nm.is_included(0,1)) {
	  nm.get_nucleus(Z,A-Z,list[ix]);
	  ix++;
	}
      } else {
	if (nm.is_included(Z,A-Z)) {
	  nm.get_nucleus(Z,A-Z,list[ix]);
	  ix++;
	}
      }
    }
  }

  if (ix!=list_size) {
    O2SCL_ERR("Sanity check in nucdist_full::set_dist().",exc_esanity);
  }
  
  return;
}

