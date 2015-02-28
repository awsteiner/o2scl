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

#include <o2scl/nucmass_hfb.h>

using namespace std;
using namespace o2scl;

nucmass_hfb::nucmass_hfb() {
  n=0;
}

nucmass_hfb::~nucmass_hfb() {
  if (n>0) {
    delete[] mass;
  }
}

double nucmass_hfb::mass_excess(int Z, int N) {
  nucmass_hfb::entry ret;
  ret=get_ZN(Z,N);
  if (ret.Z==0 && ret.N==0) return 0.0;
  return ret.Mcal;
}

int nucmass_hfb::set_data(int n_mass, nucmass_hfb::entry *m, std::string ref) {
  n=n_mass;
  mass=m;
  reference=ref;
  last=n/2;
  return 0;
}

bool nucmass_hfb::is_included(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;

  // binary search for the correct Z first
  if (mass[mid].Z!=l_Z) {
    if (mass[mid].Z>l_Z) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (mass[mp].Z<l_Z) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (mass[mid].Z!=l_Z) mid=hi;
    if (mass[mid].Z!=l_Z) {
      return false;
    }
  }

  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    return true;
  }

  int it=0;

  // Now look for the right N among all the Z's
  while (mass[mid].Z==l_Z) {

    // This hack is necessary because some nuclei are missing from
    // some of the HFB tables
    if (it>14 && mid!=0 && mass[mid].Z==mass[mid-1].Z && 
	mass[mid].N!=mass[mid-1].N+1) {
      return false;
    }
    if (mass[mid].N==l_N) {
      return true;
    } else if (mass[mid].N>l_N) {
      if (mid==0) return false;
      mid--;
    } else {
      if (mid==n-1) return false;
      mid++;
    }

    it++;
  }
  
  return false;
}

nucmass_hfb::entry nucmass_hfb::get_ZN(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;

  nucmass_hfb::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  
  // binary search for the correct Z first
  if (mass[mid].Z!=l_Z) {
    if (mass[mid].Z>l_Z) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (mass[mp].Z<l_Z) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (mass[mid].Z!=l_Z) mid=hi;
    if (mass[mid].Z!=l_Z) {
      O2SCL_ERR((((string)"Nuclei with Z=")+itos(l_Z) 
		 +" not found in nucmass_hfb::get_ZN().").c_str(),
		exc_enotfound);
    }
  }

  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    ret=mass[mid];
    last=mid;
    return ret;
  }

  // Now look for the right N among all the Z's
  while (mass[mid].Z==l_Z) {
    if (mass[mid].N==l_N) {
      ret=mass[mid];
      last=mid;
      return ret;
    } else if (mass[mid].N>l_N) {
      mid--;
    } else {
      mid++;
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)
	     +" not found in nucmass_hfb::get_ZN().").c_str(),exc_enotfound);
  return ret;
}

nucmass_hfb_sp::nucmass_hfb_sp() {
  n=0;
}

nucmass_hfb_sp::~nucmass_hfb_sp() {
  if (n>0) {
    delete[] mass;
  }
}

double nucmass_hfb_sp::mass_excess(int Z, int N) {
  nucmass_hfb_sp::entry ret;
  ret=get_ZN(Z,N);
  if (ret.Z==0 && ret.N==0) return 0.0;
  return ret.Mcal;
}

int nucmass_hfb_sp::set_data(int n_mass, nucmass_hfb_sp::entry *m, std::string ref) {
  n=n_mass;
  mass=m;
  reference=ref;
  last=n/2;
  return 0;
}

bool nucmass_hfb_sp::is_included(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;

  // binary search for the correct Z first
  if (mass[mid].Z!=l_Z) {
    if (mass[mid].Z>l_Z) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (mass[mp].Z<l_Z) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (mass[mid].Z!=l_Z) mid=hi;
    if (mass[mid].Z!=l_Z) {
      return false;
    }
  }

  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    return true;
  }

  int it=0;

  // Now look for the right N among all the Z's
  while (mass[mid].Z==l_Z) {

    // This hack is necessary because some nuclei are missing from
    // some of the HFB tables
    if (it>14 && mid!=0 && mass[mid].Z==mass[mid-1].Z && 
	mass[mid].N!=mass[mid-1].N+1) {
      return false;
    }
    if (mass[mid].N==l_N) {
      return true;
    } else if (mass[mid].N>l_N) {
      if (mid==0) return false;
      mid--;
    } else {
      if (mid==n-1) return false;
      mid++;
    }

    it++;
  }
  
  return false;
}

nucmass_hfb_sp::entry nucmass_hfb_sp::get_ZN(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;

  nucmass_hfb_sp::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  
  // binary search for the correct Z first
  if (mass[mid].Z!=l_Z) {
    if (mass[mid].Z>l_Z) {
      lo=0;
      hi=mid;
    } else {
      lo=mid;
      hi=n-1;
    }
    while (hi>lo+1) {
      int mp=(lo+hi)/2;
      if (mass[mp].Z<l_Z) {
	lo=mp;
      } else {
	hi=mp;
      }
    }
    mid=lo;
    if (mass[mid].Z!=l_Z) mid=hi;
    if (mass[mid].Z!=l_Z) {
      O2SCL_ERR((((string)"Nuclei with Z=")+itos(l_Z) 
		 +" not found in nucmass_hfb::get_ZN().").c_str(),
		exc_enotfound);
    }
  }

  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    ret=mass[mid];
    last=mid;
    return ret;
  }

  // Now look for the right N among all the Z's
  while (mass[mid].Z==l_Z) {
    if (mass[mid].N==l_N) {
      ret=mass[mid];
      last=mid;
      return ret;
    } else if (mass[mid].N>l_N) {
      mid--;
    } else {
      mid++;
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)
	     +" not found in nucmass_hfb::get_ZN().").c_str(),exc_enotfound);
  return ret;
}
