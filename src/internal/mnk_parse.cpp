/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2018, Andrew W. Steiner

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
/*
  This code generates the O2scl HDF files for the Moller et al. mass
  formula from the previous O2scl formatted files
*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// for exit()
#include <cstdlib>

#include <o2scl/string_conv.h>
#include <o2scl/nuclear_mass.h>
#include <o2scl/hdf_file.h>

#include <hdf5_hl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  ofstream fout;

  if (argc<2) {
    cout << "Usage: mnk_parse <dir>, where <dir> is the directory" << endl;
    exit(-1);
  }
  string outname="mnmsk.o2";

  vector<mnmsk_mass_entry> list;
  
  string dir=argv[1];
  string fname=dir+"/mnmsk.o2";
  ifstream fin(fname.c_str());

  string stemp;

  for(size_t j=0;j<7;j++) {
    getline(fin,stemp);
  }

  mnmsk_mass_entry me;
  string ssp, ssn;
  int N2, A2;
  for(size_t j=0;j<8979;j++) {
    fin >> me.N;
    fin >> me.A >> me.eps2 >> me.eps3 >> me.eps4;
    fin >> me.eps6 >> me.eps6sym >> me.beta2 >> me.beta3;
    fin >> me.beta4 >> me.beta6 >> me.Emic >> me.Mth >> me.Mexp;
    fin >> me.sigmaexp >> me.EmicFL >> me.MthFL >> N2;
    fin >> A2 >> ssp >> ssn >> me.gapp >> me.gapn >> me.be;
    fin >> me.S1n >> me.S2n >> me.PA >> me.PAm1 >> me.PAm2;
    fin >> me.Qbeta >> me.Tbeta >> me.S1p >> me.S2p;
    fin >> me.Qalpha >> me.Talpha;
    me.Z=me.A-me.N;
    for(size_t k=0;k<6;k++) {
      me.spinp[k]='\0';
      me.spinn[k]='\0';
    }
    for(size_t k=0;k<ssp.size();k++) me.spinp[k]=ssp[k];
    for(size_t k=0;k<ssn.size();k++) me.spinn[k]=ssn[k];
    list.push_back(me);
    if (false) {
      cout << me.N << " " << me.A << " " << me.eps2 << " " 
	   << me.eps3 << " " << me.eps4;
      cout << me.eps6 << " " << me.eps6sym << " " << me.beta2 
	   << " " << me.beta3;
      cout << me.beta4 << " " << me.beta6 << " " << me.Emic 
	   << " " << me.Mth << " " << me.Mexp;
      cout << me.sigmaexp << " " << me.EmicFL << " " << me.MthFL << " " << N2;
      cout << A2 << " " << ssp << " " << ssn << " " << me.gapp 
	   << " " << me.gapn << " " << me.be;
      cout << me.S1n << " " << me.S2n << " " << me.PA << " " 
	   << me.PAm1 << " " << me.PAm2;
      cout << me.Qbeta << " " << me.Tbeta << " " << me.S1p << " " << me.S2p;
      cout << me.Qalpha << " " << me.Talpha << endl;
      char ch;
      cin >> ch;
    }
  }

  // Make HDF table
  {
    size_t offset[34]={HOFFSET(mnmsk_mass_entry,N),
		       HOFFSET(mnmsk_mass_entry,Z),
		       HOFFSET(mnmsk_mass_entry,A),
		       HOFFSET(mnmsk_mass_entry,eps2),
		       HOFFSET(mnmsk_mass_entry,eps3),
		       HOFFSET(mnmsk_mass_entry,eps4),
		       HOFFSET(mnmsk_mass_entry,eps6),
		       HOFFSET(mnmsk_mass_entry,eps6sym),
		       HOFFSET(mnmsk_mass_entry,beta2),
		       HOFFSET(mnmsk_mass_entry,beta3),
		       HOFFSET(mnmsk_mass_entry,beta4),
		       HOFFSET(mnmsk_mass_entry,beta6),
		       HOFFSET(mnmsk_mass_entry,Emic),
		       HOFFSET(mnmsk_mass_entry,Mth),
		       HOFFSET(mnmsk_mass_entry,Mexp),
		       HOFFSET(mnmsk_mass_entry,sigmaexp),
		       HOFFSET(mnmsk_mass_entry,EmicFL),
		       HOFFSET(mnmsk_mass_entry,MthFL),
		       HOFFSET(mnmsk_mass_entry,spinp),
		       HOFFSET(mnmsk_mass_entry,spinn),
		       HOFFSET(mnmsk_mass_entry,gapp),
		       HOFFSET(mnmsk_mass_entry,gapn),
		       HOFFSET(mnmsk_mass_entry,be),
		       HOFFSET(mnmsk_mass_entry,S1n),
		       HOFFSET(mnmsk_mass_entry,S2n),
		       HOFFSET(mnmsk_mass_entry,PA),
		       HOFFSET(mnmsk_mass_entry,PAm1),
		       HOFFSET(mnmsk_mass_entry,PAm2),
		       HOFFSET(mnmsk_mass_entry,Qbeta),
		       HOFFSET(mnmsk_mass_entry,Tbeta),
		       HOFFSET(mnmsk_mass_entry,S1p),
		       HOFFSET(mnmsk_mass_entry,S2p),
		       HOFFSET(mnmsk_mass_entry,Qalpha),
		       HOFFSET(mnmsk_mass_entry,Talpha)};
    
    size_t sizes[34]={sizeof(me.N),
		      sizeof(me.Z),
		      sizeof(me.A),
		      sizeof(me.eps2),
		      sizeof(me.eps3),
		      sizeof(me.eps4),
		      sizeof(me.eps6),
		      sizeof(me.eps6sym),
		      sizeof(me.beta2),
		      sizeof(me.beta3),
		      sizeof(me.beta4),
		      sizeof(me.beta6),
		      sizeof(me.Emic),
		      sizeof(me.Mth),
		      sizeof(me.Mexp),
		      sizeof(me.sigmaexp),
		      sizeof(me.EmicFL),
		      sizeof(me.MthFL),
		      sizeof(me.spinp),
		      sizeof(me.spinn),
		      sizeof(me.gapp),
		      sizeof(me.gapn),
		      sizeof(me.be),
		      sizeof(me.S1n),
		      sizeof(me.S2n),
		      sizeof(me.PA),
		      sizeof(me.PAm1),
		      sizeof(me.PAm2),
		      sizeof(me.Qbeta),
		      sizeof(me.Tbeta),
		      sizeof(me.S1p),
		      sizeof(me.S2p),
		      sizeof(me.Qalpha),
		      sizeof(me.Talpha)
    };
    
    const char *names[34]={
      "Neutron number",
      "Proton number",
      "Atomic number",
      "Quadrupole deformation",
      "Octupole deformation",
      "Hexadecapole deformation",
      "Hexacontatetrapole deformation",
      "Hexacontatetrapole deformation without mass asymmetry",
      "Quadrupole deformation (SH)",
      "Octupole deformation (SH)",
      "Hexadecapole deformation (SH)",
      "Hexacontatetrapole deformation (SH)",
      "The ground-state microscopic energy",
      "The theoretical mass excess (in MeV)",
      "The experimental mass excess (in MeV)",
      "Experimental mass excess error",
      "The ground-state microscopic energy in the FRLDM",
      "The theoretical mass excess in the FRLDM",
      "Spin and pairity of odd proton ",
      "Spin and pairity of odd neutron",
      "Lipkin-Nogami proton gap",
      "Lipkin-Nogami neutron gap",
      "Total binding energy",
      "One neutron separation energy",
      "Two neutron separation energy",
      "Percentage of daughters in beta decay before first neutron",
      "Percentage of daughters in beta decay before second neutron",
      "Percentage of daughters in beta decay before third neutron",
      "Energy released in beta-decay",
      "Half-life w.r.t. GT beta-decay",
      "One proton separation energy",
      "Two proton separation energy",
      "Energy released in alpha-decay",
      "Half-life w.r.t. alpha-decay"};
      
    // string type
    hid_t string_type6=H5Tcopy(H5T_C_S1);
    H5Tset_size(string_type6,6);

    hid_t field_type[34]={H5T_NATIVE_INT,H5T_NATIVE_INT,H5T_NATIVE_INT,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,string_type6,string_type6,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE};

    hid_t file=H5Fcreate(outname.c_str(),H5F_ACC_TRUNC,
			 H5P_DEFAULT,H5P_DEFAULT);
      
    hdf_file hf;
    hf.set_current_id(file);
    hf.seti("nrecords",list.size());
    hf.sets_fixed("comment",
		  ((string)"HDF5 version of M\"{o}ller, et al. ")+
		  "mass data created for O2scl. "
		  "See http://o2scl.sourceforge.net for details.");
      
    herr_t status=H5TBmake_table("mnk97.dat",file,outname.c_str(),
				 34,list.size(),sizeof(mnmsk_mass_entry),
				 names,offset,field_type,100,0,0,&list[0]);

    hf.sets("orig_file","mnk97.dat");
    hf.sets("reference",
	    ((string)"P. M\"{o}ller, J. R. Nix, W. D. Myers, and W. J. ")+
	    "Swiatecki, Atomic Data Nucl. Data Tables 59 (1995) 185.");
      
    H5Tclose(string_type6);
    H5Fclose(file);
  }

  return 0;
}
