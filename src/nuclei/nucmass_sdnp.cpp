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
#include <o2scl/nucmass_sdnp.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_sdnp::nucmass_sdnp() {
  n=0;
}
  
int nucmass_sdnp::load(std::string model, bool external) {
  
  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    if (model=="sdnp03") {
      fname=dir+"/nucmass/sdnp03.o2";
    } else if (model=="sd_skp_04") {
      fname=dir+"/nucmass/sd_skp_04.o2";
    } else if (model=="sd_sly4_04") {
      fname=dir+"/nucmass/sd_sly4_04.o2";
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
  
  mass=new nucmass_sdnp::entry[n];
  for(size_t i=0;i<n;i++) {
    nucmass_sdnp::entry nde={((int)(data.get("Z",i)+1.0e-6)),
			     ((int)(data.get("N",i)+1.0e-6)),
			     data.get("ENERGY",i)};
    mass[i]=nde;
  }

  last=n/2;
  return 0;
}

nucmass_sdnp::~nucmass_sdnp() {
}

bool nucmass_sdnp::is_included(int l_Z, int l_N) {
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

double nucmass_sdnp::mass_excess(int l_Z, int l_N) {
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
		 " not found in nucmass_sdnp::mass_excess().").c_str(),
		exc_enotfound);
    }
  }
  
  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    int A=l_Z+l_N;
    return mass[mid].ENERGY-A*m_amu+l_Z*(m_prot+m_elec)+l_N*m_neut;
  }
  
  // Now look for the right N among all the N's
  while (mass[mid].Z==l_Z) {
    
    if (mass[mid].N==l_N) {
      int A=l_Z+l_N;
      return mass[mid].ENERGY-A*m_amu+l_Z*(m_prot+m_elec)+l_N*m_neut;
    } else if (mass[mid].N>l_N) {
      if (mid==0) {
	O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		   " not found in nucmass_sdnp::mass_excess().").c_str(),
		  exc_enotfound);
      }
      mid--;
    } else {
      if (mid==((int)n-1)) {
	O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		   " not found in nucmass_sdnp::mass_excess().").c_str(),
		  exc_enotfound);
      }
      mid++;
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
	     " not found in nucmass_sdnp::mass_excess().").c_str(),
	    exc_enotfound);
  
  return 0.0;
}
