/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2012-2025, Andrew W. Steiner
  
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
#include <o2scl/nucmass_gen.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_gen::nucmass_gen() {
}

nucmass_gen::~nucmass_gen() {
}

int nucmass_gen::load_be(std::string fname, std::string be_col,
			 double be_units, bool external) {

  if (n>0) {
    data.clear();
    n=0;
  }
  
  std::string dir=o2scl::o2scl_settings.get_data_dir();
  if (!external) {
    if (fname=="ddme2" || fname=="ddmed" || fname=="ddpc1" || fname=="nl3s") {
      this->reference=((std::string)"S. E. Agbemava et al., ")+
        "Phys. Rev. C 89, 054320 (2014); A. V. Afanasjev et al., "+
        "Phys. Rev. C 91, 014324 (2015); A. V. Afanasjev and S. E. "+
        "Agbemava, Phys. Rev. C 93, 054310 (2016).";
    } else if (fname=="skms_all" || fname=="skms_even") {
      reference="Bartel et al., Nucl. Phys. A 386, 79 (1982)";
    } else if (fname=="skp_all" || fname=="skp_even") {
      reference="Dobaczewski et al., Nucl. Phys. A 422, 103 (1984).";
    } else if (fname=="sly4_all" || fname=="sly4_even") {
      reference="Chabanat et al., Nucl. Phys. A 635, 231 (1998).";
    } else if (fname=="sv-min_all" || fname=="sv-min_even") {
      reference="Klupfel et al., Phys. Rev. C 79, 034310 (2009).";
    } else if (fname=="unedf0_all" || fname=="unedf0_even") {
      reference="Kortelainen et al., Phys. Rev. C 82, 024313 (2010).";
    } else if (fname=="unedf1_all" || fname=="unedf1_even") {
      reference="Kortelainen et al., Phys. Rev. C 85, 024304 (2012).";
    }
    if (fname.substr(fname.length()-3,3)!=".o2") {
      fname=fname+".o2";
    }
    fname=dir+"/nucmass/frib_mex/"+fname;
  }
  
  o2scl_hdf::hdf_file hf;
  hf.open(fname);
  string name;
  hdf_input_n(hf,data,name);
  hf.close();

  if (data.is_column("N")==false) {
    data.function_column("A-Z","N");
  }
  data.sort_table("Z");

  n=data.get_nlines();

  data.new_column("mex");
  for(size_t i=0;i<n;i++) {
    double be=data.get(be_col,i);
    double dZ=data.get("Z",i);
    double dN=data.get("N",i);
    double mex=be-(dZ+dN)*m_amu+dZ*m_elec+dN*m_neut+dZ*m_prot;
    data.set("mex",i,mex);
  }
  mex_col_ix=data.lookup_column("mex");
  
  last=n/2;

  return 0;
}

bool nucmass_gen::is_included(int l_Z, int l_N) {

  for(size_t i=0;i<n;i++) {
    if (fabs(data.get("Z",i)-l_Z)<1.0e-12 &&
	fabs(data.get("N",i)-l_N)<1.0e-12) {
      return true;
    }
  }

  return false;
}

double nucmass_gen::mass_excess(int l_Z, int l_N) {
  for(size_t i=0;i<n;i++) {
    if (fabs(data.get("Z",i)-l_Z)<1.0e-12 &&
	fabs(data.get("N",i)-l_N)<1.0e-12) {
      return data.get(mex_col_ix,i);
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
	     " not found in nucmass_gen::mass_excess().").c_str(),
	    exc_enotfound);
  
  return 0.0;
}

double nucmass_gen::get_string(int l_Z, int l_N, std::string column) {
  for(size_t i=0;i<n;i++) {
    if (fabs(data.get("Z",i)-l_Z)<1.0e-12 &&
	fabs(data.get("N",i)-l_N)<1.0e-12) {
      return data.get(column,i);
    }
  }
  
  O2SCL_ERR((((string)"Nucleus with Z=")+itos(l_Z)+" and N="+itos(l_N)+
	     " not found in nucmass_gen::mass_excess().").c_str(),
	    exc_enotfound);
  
  return 0.0;
}

