/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#include <o2scl/nucmass.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_ame::nucmass_ame() {
  n=0;
  reference="";
  mass=0;
  last=0;
}

nucmass_ame::~nucmass_ame() {
  if (n>0) {
    delete[] mass;
  }
}

double nucmass_ame::mass_excess(int Z, int N) {
  entry ret;
  ret=get_ZN(Z,N);
  if (ret.Z==0 && ret.N==0) return 0.0;
  return ret.mass/1.0e3;
}

bool nucmass_ame::is_included(int l_Z, int l_N) {

  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::is_included().",
		  exc_einval);
  }

  if (false) {

    int l_A=l_Z+l_N, mid=last, lo=0, hi=0;

    // Binary search for the correct A first
    if (mass[mid].Z+mass[mid].N!=l_A) {
      if (mass[mid].Z+mass[mid].N>l_A) {
	lo=0;
	hi=mid;
      } else {
	lo=mid;
	hi=n-1;
      }
      while (hi>lo+1) {
	int mp=(lo+hi)/2;
	if (mass[mid].Z+mass[mid].N<l_A) {
	  lo=mp;
	} else {
	  hi=mp;
	}
      }
      mid=lo;
      if (mass[mid].Z+mass[mid].N!=l_A) {
	mid=hi;
      }
      if (mass[mid].Z+mass[mid].N!=l_A) {
	return false;
      }
    }

    // The cached point was the right one, so we're done
    if (mass[mid].N==l_N) {
      if (fabs(mass[mid].mass)>1.0e-20 &&
	  fabs(mass[mid].mass)<1.0e90) {
	return true;
      } else {
	return false;
      }
    }

    // Now look for the right N among all the A's
    while (mass[mid].Z+mass[mid].N==l_A) {
      if (mass[mid].N==l_N) {
	if (fabs(mass[mid].mass)>1.0e-20 &&
	    fabs(mass[mid].mass)<1.0e90) {
	  return true;
	} else {
	  return false;
	}
      } else if (mass[mid].N>l_N) {
	if (mid==0) return false;
	mid--;
      } else {
	if (mid==n-1) return false;
	mid++;
      }
    }

    return false;
  }

  for(int i=0;i<n;i++) {
    if (mass[i].Z==l_Z && mass[i].N==l_N) {
      return true;
    }
  }
  return false;
}

bool nucmass_ame_exp::is_included(int l_Z, int l_N) {
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame_exp::is_included().",
		  exc_einval);
  }
  for(int i=0;i<n;i++) {
    if (mass[i].Z==l_Z && mass[i].N==l_N && mass[i].mass_acc==0) {
      return true;
    }
  }
  return false;
}

nucmass_ame::entry nucmass_ame::get_ZN(int l_Z, int l_N) {
  nucmass_ame::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::get_ZN().",
	      exc_einval);
    return ret;
  }
  for(int i=0;i<n;i++) {
    if (mass[i].Z==l_Z && mass[i].N==l_N) {
      ret=mass[i];
    }
  }
  return ret;
}

nucmass_ame::entry nucmass_ame::get_ZA(int l_Z, int l_A) {
  nucmass_ame::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::get_ZA().",
	      exc_einval);
    return ret;
  }
  for(int i=0;i<n;i++) {
    if (mass[i].Z==l_Z && mass[i].A==l_A) {
      ret=mass[i];
    }
  }
  return ret;
}

nucmass_ame::entry nucmass_ame::get_elA(string l_el, int l_A) {
  nucmass_ame::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::get_elA().",
	      exc_einval);
    return ret;
  }
  for(int i=0;i<n;i++) {
    if (mass[i].el==l_el && mass[i].A==l_A) {
      ret=mass[i];
    }
  }
  return ret;
}

nucmass_ame::entry nucmass_ame::get(string nucleus) {
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::get().",
	      exc_einval);
    nucmass_ame::entry ret;
    ret.Z=0;
    ret.A=0;
    ret.N=0;
    return ret;
  }
  if (isalpha(nucleus[1])) {
    string el=nucleus.substr(0,2);
    int A=o2scl::stoi(nucleus.substr(2,nucleus.length()-2));
    return get_elA(el,A);
  } 
  string el=nucleus.substr(0,1);
  int A=o2scl::stoi(nucleus.substr(1,nucleus.length()-1));
  return get_elA(el,A);
}

