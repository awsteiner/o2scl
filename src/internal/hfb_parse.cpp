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
/*
  This code generates the O2scl HDF files for the HFB mass formula
  from the previous O2scl formatted files.
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

  if (argc<3) {
    cout << "Usage: hfb_parse <dir> <file>." << endl;
    exit(-1);
  }
  string outname=argv[2];

  vector<hfb_mass_entry> list;
  vector<hfb_sp_mass_entry> list_sp;
  
  string dir=argv[1];
  string fname=dir+argv[2];
  ifstream fin(fname.c_str());
  size_t jmax;
  string stemp, orig_file;
  bool inc_spin_parity=false;
  
  if ((string)(argv[2])=="hfb2.o2") {
    cout << fname << endl;
    //for(size_t j=0;j<6;j++) {
    //getline(fin,stemp);
    //cout << stemp << endl;
    //}
    jmax=9203;
    orig_file="hfb2-plain";
  } else if ((string)(argv[2])=="hfb8.o2") {
    cout << fname << endl;
    //for(size_t j=0;j<6;j++) {
    //getline(fin,stemp);
    //cout << stemp << endl;
    //}
    jmax=9197;
    orig_file="hfb8-plain";
  } else if ((string)(argv[2])=="hfb14.o2") {
    cout << fname << endl;
    //for(size_t j=0;j<9;j++) {
    //getline(fin,stemp);
    //cout << stemp << endl;
    //}
    jmax=8388;
    orig_file="hfb14-plain";
  } else if ((string)(argv[2])=="hfb17.o2") {
    cout << fname << endl;
    //for(size_t j=0;j<9;j++) {
    //getline(fin,stemp);
    //cout << stemp << endl;
    //}
    jmax=8389;
    orig_file="hfb17-plain";
    inc_spin_parity=true;
  } else {
    cout << fname << endl;
    //for(size_t j=0;j<9;j++) {
    //getline(fin,stemp);
    //cout << stemp << endl;
    //}
    jmax=8389;
    orig_file="hfb21-dat";
    inc_spin_parity=true;
  }
  
  if (inc_spin_parity) {

  hfb_sp_mass_entry he;
  string ssp, ssn;
  int N2, A2;
  for(size_t j=0;j<jmax;j++) {
    double dtemp;
    fin >> dtemp;
    he.Z=((int)(dtemp+1.0e-4));
    fin >> dtemp;
    he.A=((int)(dtemp+1.0e-4));
    he.N=he.A-he.Z;
    fin >> he.bet2;
    fin >> he.bet4;
    fin >> he.Rch;
    fin >> he.def_wig;
    fin >> he.Sn;
    fin >> he.Sp;
    fin >> he.Qbet;
    fin >> he.Mcal;
    fin >> he.Err;
    fin >> he.Jexp;
    fin >> he.Jth;
    fin >> he.Pexp;
    fin >> he.Pth;
    if (j==0 || j==jmax-1) {
      cout << j << " " << he.Z << " " << he.A << " " << he.N << endl;
    }
    list.push_back(he);
  }


    // Make HDF table
    size_t offset[16]={HOFFSET(hfb_mass_entry,N),
		       HOFFSET(hfb_mass_entry,Z),
		       HOFFSET(hfb_mass_entry,A),
		       HOFFSET(hfb_mass_entry,bet2),
		       HOFFSET(hfb_mass_entry,bet4),
		       HOFFSET(hfb_mass_entry,Rch),
		       HOFFSET(hfb_mass_entry,def_wig),
		       HOFFSET(hfb_mass_entry,Sn),
		       HOFFSET(hfb_mass_entry,Sp),
		       HOFFSET(hfb_mass_entry,Qbet),
		       HOFFSET(hfb_mass_entry,Mcal),
		       HOFFSET(hfb_mass_entry,Err),
		       HOFFSET(hfb_mass_entry,Jexp),
		       HOFFSET(hfb_mass_entry,Jth),
		       HOFFSET(hfb_mass_entry,Pexp),
		       HOFFSET(hfb_mass_entry,Pth)};
    
    size_t sizes[16]={sizeof(he.N),
		      sizeof(he.Z),
		      sizeof(he.A),
		      sizeof(he.bet2),
		      sizeof(he.bet4),
		      sizeof(he.Rch),
		      sizeof(he.def_wig),
		      sizeof(he.Sn),
		      sizeof(he.Sp),
		      sizeof(he.Qbet),
		      sizeof(he.Mcal),
		      sizeof(he.Err),
		      sizeof(he.Jexp),
		      sizeof(he.Jth),
		      sizeof(he.Pexp),
		      sizeof(he.Pth)};
    
    const char *names[16]={
      "Neutron number",
      "Proton number",
      "Atomic number",
      "Beta 2 deformation",
      "Beta 4 deformation",
      "RMS charge radius",
      "Deformation and Wigner energies",
      "Neutron separation energy",
      "Proton separation energy",
      "Beta-decay energy",
      "Calculated mass excess",
      "Error between experimental and calculated mass excess",
      "Experimental spin","Theoretical spin",
      "Experimental parity","Theoretical parity"};
      
    hid_t field_type[16]={H5T_NATIVE_INT,H5T_NATIVE_INT,H5T_NATIVE_INT,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_INT,H5T_NATIVE_INT};

    hid_t file=H5Fcreate(outname.c_str(),H5F_ACC_TRUNC,
			 H5P_DEFAULT,H5P_DEFAULT);
      
    hdf_file hf;
    hf.set_current_id(file);
    hf.seti("nrecords",list.size());
    cout << "nrecords: " << list.size() << endl;
    hf.sets_fixed("comment",
		  ((string)"HDF5 version of HFB ")+
		  "mass data created for O2scl. "
		  "See http://o2scl.sourceforge.net for details.");
    
    herr_t status=H5TBmake_table(orig_file.c_str(),file,outname.c_str(),
				 16,list.size(),sizeof(hfb_mass_entry),
				 names,offset,field_type,100,0,0,&list[0]);

    hf.sets("orig_file",orig_file);
    if ((string)(argv[2])=="hfb17.o2") {
      hf.sets("reference",
	      ((string)"S. Goriely, M. Samyn, P.-H. Heenen, J.M. Pearson, ")+
	      "and F. Tondeur, Phys. Rev. C 66 (2002) 024326.");
    } else {
      hf.sets("reference",
	      ((string)"S. Goriely, M. Samyn, and J. M. Pearson, ")+
	      "Phys. Rev. C 75 (2007) 064312.");
    }
    H5Fclose(file);

  } else {
    
    hfb_mass_entry he;
    string ssp, ssn;
    int N2, A2;
    for(size_t j=0;j<jmax;j++) {
      double dtemp;
      fin >> dtemp;
      he.Z=((int)(dtemp+1.0e-4));
      fin >> dtemp;
      he.A=((int)(dtemp+1.0e-4));
      he.N=he.A-he.Z;
      fin >> he.bet2;
      fin >> he.bet4;
      fin >> he.Rch;
      fin >> he.def_wig;
      fin >> he.Sn;
      fin >> he.Sp;
      fin >> he.Qbet;
      fin >> he.Mcal;
      fin >> he.Err;
      if (j==0 || j==jmax-1) {
	cout << j << " " << he.Z << " " << he.A << " " << he.N << endl;
      }
      list.push_back(he);
    }

    // Make HDF table
    size_t offset[12]={HOFFSET(hfb_mass_entry,N),
		       HOFFSET(hfb_mass_entry,Z),
		       HOFFSET(hfb_mass_entry,A),
		       HOFFSET(hfb_mass_entry,bet2),
		       HOFFSET(hfb_mass_entry,bet4),
		       HOFFSET(hfb_mass_entry,Rch),
		       HOFFSET(hfb_mass_entry,def_wig),
		       HOFFSET(hfb_mass_entry,Sn),
		       HOFFSET(hfb_mass_entry,Sp),
		       HOFFSET(hfb_mass_entry,Qbet),
		       HOFFSET(hfb_mass_entry,Mcal),
		       HOFFSET(hfb_mass_entry,Err)};
    
    size_t sizes[12]={sizeof(he.N),
		      sizeof(he.Z),
		      sizeof(he.A),
		      sizeof(he.bet2),
		      sizeof(he.bet4),
		      sizeof(he.Rch),
		      sizeof(he.def_wig),
		      sizeof(he.Sn),
		      sizeof(he.Sp),
		      sizeof(he.Qbet),
		      sizeof(he.Mcal),
		      sizeof(he.Err)};
    
    const char *names[12]={
      "Neutron number",
      "Proton number",
      "Atomic number",
      "Beta 2 deformation",
      "Beta 4 deformation",
      "RMS charge radius",
      "Deformation and Wigner energies",
      "Neutron separation energy",
      "Proton separation energy",
      "Beta-decay energy",
      "Calculated mass excess",
      "Error between experimental and calculated mass excess"};
      
    hid_t field_type[12]={H5T_NATIVE_INT,H5T_NATIVE_INT,H5T_NATIVE_INT,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			  H5T_NATIVE_DOUBLE};

    hid_t file=H5Fcreate(outname.c_str(),H5F_ACC_TRUNC,
			 H5P_DEFAULT,H5P_DEFAULT);
      
    hdf_file hf;
    hf.set_current_id(file);
    hf.seti("nrecords",list.size());
    cout << "nrecords: " << list.size() << endl;
    hf.sets_fixed("comment",
		  ((string)"HDF5 version of HFB ")+
		  "mass data created for O2scl. "
		  "See http://o2scl.sourceforge.net for details.");
    
    herr_t status=H5TBmake_table(orig_file.c_str(),file,outname.c_str(),
				 12,list.size(),sizeof(hfb_mass_entry),
				 names,offset,field_type,100,0,0,&list[0]);

    hf.sets("orig_file",orig_file);
    if ((string)(argv[2])=="hfb2.o2") {
      hf.sets("reference",
	      ((string)"S. Goriely, M. Samyn, P.-H. Heenen, J.M. Pearson, ")+
	      "and F. Tondeur, Phys. Rev. C 66 (2002) 024326.");
    } else if ((string)(argv[2])=="hfb8.o2") {
      hf.sets("reference",
	      ((string)"M. Samyn, S. Goriely, M. Bender and J.M. Pearson, ")+
	      "Phys. Rev. C 70 (2004) 044309.");
    } else {
      hf.sets("reference",
	      ((string)"S. Goriely, M. Samyn, and J. M. Pearson, ")+
	      "Phys. Rev. C 75 (2007) 064312.");
    }
    H5Fclose(file);
  }

  return 0;
}
