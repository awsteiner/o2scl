/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/nuclear_dist.h>

using namespace std;
using namespace o2scl;

bool o2scl::operator==(const nuclear_dist::iterator &i1,
		      const nuclear_dist::iterator &i2) {
  return i1.np==i2.np;
}

bool o2scl::operator!=(const nuclear_dist::iterator &i1,
	       const nuclear_dist::iterator &i2) {
  return i1.np!=i2.np;
}

full_dist::full_dist(nuclear_mass &nm, int maxA, 
		     bool include_neutron) {
  nucleus n;
  list_size=0;
  if (include_neutron && nm.is_included(0,1)) {
    list_size++;
  }
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
  list=new nucleus[list_size];
  int ix=0;
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
}

int full_dist::set_dist(nuclear_mass &nm, int maxA, 
			bool include_neutron) {
  
  if (list_size>0) delete[] list;

  nucleus n;
  list_size=0;
  if (include_neutron && nm.is_included(0,1)) {
    list_size++;
  }
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
  list=new nucleus[list_size];
  int ix=0;
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
  
  return 0;
}

