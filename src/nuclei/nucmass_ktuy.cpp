/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/nucmass_ktuy.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_ktuy::nucmass_ktuy() {
  n=0;
}

int nucmass_ktuy::load(std::string model, bool external) {
  
  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    if (model=="04") {
      fname=dir+"/nucmass/ktuy04.o2";
    } else {
      fname=dir+"/nucmass/ktuy05.o2";
    }
  }
  
  table<> data;
  o2scl_hdf::hdf_file hf;
  hf.open(fname);
  string name;
#ifndef O2SCL_NO_HDF_INPUT  
  hdf_input_n(hf,data,name);
#endif
  hf.close();
  
  n=data.get_nlines();

  mass=new nucmass_ktuy::entry[n];
  for(size_t i=0;i<n;i++) {
    nucmass_ktuy::entry kme={((int)(data.get("NN",i)+1.0e-6)),
			 ((int)(data.get("ZZ",i)+1.0e-6)),
			 ((int)(data.get("NN",i)+data.get("ZZ",i)+1.0e-6)),
			 data.get("Mcal",i),data.get("Esh",i),
			 data.get("alpha2",i),data.get("alpha4",i),
			 data.get("alpha6",i)};
    mass[i]=kme;
  }

  last=n/2;

  return 0;
}

nucmass_ktuy::~nucmass_ktuy() {
  if (n>0) {
    delete[] mass;
  }
}

bool nucmass_ktuy::is_included(int l_Z, int l_N) {
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
    // some of the KTUY tables
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
      if (mid==((int)n-1)) return false;
      mid++;
    }

    it++;
  }
  
  return false;
}

nucmass_ktuy::entry nucmass_ktuy::get_ZN(int l_Z, int l_N) {
  int lo=0, hi=0, mid=last;

  nucmass_ktuy::entry ret;
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
		 +" not found in nucmass_ktuy::get_ZN().").c_str(),
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
	     +" not found in nucmass_ktuy::get_ZN().").c_str(),exc_enotfound);
  return ret;
}

double nucmass_ktuy::mass_excess(int Z, int N) {
  nucmass_ktuy::entry ret;
  ret=get_ZN(Z,N);
  if (ret.Z==0 && ret.N==0) return 0.0;
  return ret.Mcal;
}

