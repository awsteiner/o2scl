/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2020, Andrew W. Steiner
  
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
#include <o2scl/nucmass_dglg.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_dglg::nucmass_dglg(std::string model, bool external) {
  n=0;
  
  std::string fname;
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (external) {
    fname=model;
  } else {
    fname=dir+"/nucmass/dglg10.o2";
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
  
  mass=new nucmass_dglg::entry[n];
  for(size_t i=0;i<n;i++) {
    nucmass_dglg::entry nde={((int)(data.get("Z",i)+1.0e-6)),
			     ((int)(data.get("N",i)+1.0e-6)),
			     data.get("EHFB",i),data.get("BMIN",i),
			     data.get("GMIN",i),data.get("RCHFB",i),
			     data.get("RPHFB",i),data.get("RNHFB",i),
			     data.get("EABS",i),data.get("ECORR",i),
			     data.get("BET01",i),data.get("GAM01",i),
			     data.get("DELB01",i),data.get("DELG01",i),
			     data.get("E21",i),data.get("E41",i),
			     data.get("E61",i),data.get("E02",i),
			     data.get("E22",i),data.get("E23",i),
			     data.get("PK0_2_1",i),data.get("PK2_2_2",i),
			     data.get("PK2_2_3",i),data.get("BE2_2_1_0_1",i),
			     data.get("BE2_2_3_0_1",i),
			     data.get("BE2_2_1_0_2",i),
			     data.get("BE2_4_1_2_1",i),
			     data.get("BE2_2_3_2_1",i),
			     data.get("BE2_2_3_0_2",i),
			     data.get("RC5DCH",i),data.get("RP5DCH",i),
			     data.get("RN5DCH",i),data.get("ROE0TH",i),
			     ((int)(data.get("NMIN",i)+1.0e-6)),
			     ((int)(data.get("NMAX",i)+1.0e-6))};
    mass[i]=nde;
  }

  last=n/2;
}

nucmass_dglg::~nucmass_dglg() {
}

bool nucmass_dglg::is_included(int l_Z, int l_N) {
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

double nucmass_dglg::mass_excess(int l_Z, int l_N) {
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
		 " not found in nucmass_dglg::mass_excess().").c_str(),
		exc_enotfound);
    }
  }
  
  // The cached point was the right one, so we're done
  if (mass[mid].N==l_N) {
    int A=l_Z+l_N;
    return mass[mid].EHFB-A*m_amu+l_Z*(m_prot+m_elec)+l_N*m_neut;
  }
  
  // Now look for the right N among all the N's
  while (mass[mid].Z==l_Z) {
    
    if (mass[mid].N==l_N) {
      int A=l_Z+l_N;
      return mass[mid].EHFB-A*m_amu+l_Z*(m_prot+m_elec)+l_N*m_neut;
    } else if (mass[mid].N>l_N) {
      if (mid==0) {
	O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		   " not found in nucmass_dglg::mass_excess().").c_str(),
		  exc_enotfound);
      }
      mid--;
    } else {
      if (mid==((int)n-1)) {
	O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
		   " not found in nucmass_dglg::mass_excess().").c_str(),
		  exc_enotfound);
      }
      mid++;
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
	     " not found in nucmass_dglg::mass_excess().").c_str(),
	    exc_enotfound);
  
  return 0.0;
}
