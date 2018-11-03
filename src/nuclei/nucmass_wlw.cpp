/*
  -------------------------------------------------------------------
  
  Copyright (C) 2014, Andrew W. Steiner
  
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
#include <o2scl/nucmass_wlw.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_wlw::nucmass_wlw() {
  n=0;
}
  
int nucmass_wlw::load(std::string model, bool external) {
  
  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    if (model=="WS3.2") {
      fname=dir+"/nucmass/wlw10.o2";
    } else if (model=="WS3.3") {
      fname=dir+"/nucmass/wllw10.o2";
    } else if (model=="WS3.6") {
      fname=dir+"/nucmass/lwdw11.o2";
    } else if (model=="WS3_RBF") {
      fname=dir+"/nucmass/wl11.o2";
    } else if (model=="WS4_RBF") {
      fname=dir+"/nucmass/wlwm14.o2";
    } else {
      O2SCL_ERR("Invalid model in nucmass_sdnp().",exc_einval);
    }
  }
  
  table<> data;
  o2scl_hdf::hdf_file hf;
  hf.open(fname);
  string name;
#ifndef O2SCL_NO_HDF_INPUT  
  hdf_input(hf,data,name);
#endif
  hf.close();
  
  n=data.get_nlines();
  
  mass=new nucmass_wlw::entry[n];
  if (model=="WS3_RBF") {
    for(size_t i=0;i<n;i++) {
      nucmass_wlw::entry nde={((int)(data.get("Z",i)+1.0e-6)),
			      ((int)(data.get("A",i)+1.0e-6))-
			      ((int)(data.get("Z",i))),
			      data.get("WS3_RBF",i)};
      mass[i]=nde;
    }
  } else if (model=="WS4_RBF") {
    for(size_t i=0;i<n;i++) {
      nucmass_wlw::entry nde={((int)(data.get("Z",i)+1.0e-6)),
			      ((int)(data.get("A",i)+1.0e-6))-
			      ((int)(data.get("Z",i))),
			      data.get("WS4_RBF",i)};
      mass[i]=nde;
    }
  } else {
    for(size_t i=0;i<n;i++) {
      nucmass_wlw::entry nde={((int)(data.get("Z",i)+1.0e-6)),
			      ((int)(data.get("A",i)+1.0e-6))-
			      ((int)(data.get("Z",i))),
			      data.get("Mth",i)};
      mass[i]=nde;
    }
  }

  last=n/2;
  return 0;
}

nucmass_wlw::~nucmass_wlw() {
}

bool nucmass_wlw::is_included(int l_Z, int l_N) {
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

  // Now look for the right N among all the N's
  while (mass[mid].Z==l_Z) {

    if (mass[mid].N==l_N) {
      return true;
    } else if (mass[mid].N>l_N) {
      if (mid==0) return false;
      mid--;
    } else {
      if (mid==((int)n-1)) return false;
      mid++;
    }
  }
  
  return false;
}

double nucmass_wlw::mass_excess(int l_Z, int l_N) {
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
      O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		 " not found in nucmass_wlw::mass_excess().").c_str(),
		exc_enotfound);
    }
  }
  
  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    int A=l_Z+l_N;
    return mass[mid].Mth;
  }
  
  // Now look for the right N among all the N's
  while (mass[mid].Z==l_Z) {
    
    if (mass[mid].N==l_N) {
      int A=l_Z+l_N;
      return mass[mid].Mth;
    } else if (mass[mid].N>l_N) {
      if (mid==0) {
	O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		   " not found in nucmass_wlw::mass_excess().").c_str(),
		  exc_enotfound);
      }
      mid--;
    } else {
      if (mid==((int)n-1)) {
	O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		   " not found in nucmass_wlw::mass_excess().").c_str(),
		  exc_enotfound);
      }
      mid++;
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
	     " not found in nucmass_wlw::mass_excess().").c_str(),
	    exc_enotfound);
  
  return 0.0;
}
