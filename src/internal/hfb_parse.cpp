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
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_hfb.h>
#include <o2scl/hdf_file.h>

#include <hdf5_hl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  ofstream fout;

  if (argc<3) {
    cout << "Usage: hfb_parse <dir> <out file>." << endl;
    exit(-1);
  }
  string outname=argv[2];

  vector<nucmass_hfb_entry> list;
  vector<nucmass_hfb_sp_entry> list_sp;
  
  string dir=argv[1];
  string out_fname=argv[2];
  string orig_file;
  size_t jmax;
  string stemp;
  bool inc_spin_parity=false;
  
  ifstream fin;

  if ((string)(argv[2])=="hfb2.o2") {
    orig_file="hfb2-plain";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<4;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=9203;
  } else if ((string)(argv[2])=="hfb8.o2") {
    orig_file="hfb8-plain";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=9197;
  } else if ((string)(argv[2])=="hfb14.o2") {
    orig_file="hfb14-plain";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=8388;
  } else if ((string)(argv[2])=="hfb14_v0.o2") {
    orig_file="hfb14-plain_v0";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=8382;
  } else if ((string)(argv[2])=="hfb17.o2") {
    orig_file="hfb17-plain";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=8389;
    inc_spin_parity=true;
  } else if ((string)(argv[2])=="hfb21.o2") {
    orig_file="hfb21-dat";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=8387;
    inc_spin_parity=true;
  } else if ((string)(argv[2])=="hfb22.o2") {
    orig_file="hfb22-dat";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=8392;
    inc_spin_parity=true;
  } else if ((string)(argv[2])=="hfb23.o2") {
    orig_file="hfb23-dat";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=8392;
    inc_spin_parity=true;
  } else if ((string)(argv[2])=="hfb24.o2") {
    orig_file="hfb24-dat";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=8388;
    inc_spin_parity=true;
  } else if ((string)(argv[2])=="hfb25.o2") {
    orig_file="hfb25-dat";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=9484;
    inc_spin_parity=true;
  } else if ((string)(argv[2])=="hfb26.o2") {
    orig_file="hfb26-dat";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=9511;
    inc_spin_parity=true;
  } else if ((string)(argv[2])=="hfb27.o2") {
    orig_file="hfb27-dat";
    string in_fname=dir+"/"+orig_file;
    cout << "Opening file '" << in_fname << "'." << endl;
    fin.open(in_fname.c_str());
    for(size_t j=0;j<3;j++) {
      getline(fin,stemp);
      cout << stemp << endl;
    }
    jmax=9482;
    inc_spin_parity=true;
  } else {
    O2SCL_ERR("Bad argument 2.",exc_efailed);
  }

  if (!inc_spin_parity) {
    
    nucmass_hfb_entry he;
    string ssp, ssn;
    int N2, A2;
    for(size_t j=0;j<jmax;j++) {
      getline(fin,stemp);
      he.Z=o2scl::stoi(stemp.substr(0,4));
      he.A=o2scl::stoi(stemp.substr(4,4));
      he.N=he.A-he.Z;
      if (orig_file!=((string)"hfb14-plain_v0")) {
	if (o2scl::count_words(stemp.substr(8,6))==0) {
	  he.bet2=1.0e99;
	} else {
	  he.bet2=o2scl::stod(stemp.substr(8,6));
	}
	if (o2scl::count_words(stemp.substr(14,6))==0) {
	  he.bet4=1.0e99;
	} else {
	  he.bet4=o2scl::stod(stemp.substr(14,6));
	}
	if (o2scl::count_words(stemp.substr(20,6))==0) {
	  he.Rch=1.0e99;
	} else {
	  he.Rch=o2scl::stod(stemp.substr(20,6));
	}
	if (o2scl::count_words(stemp.substr(26,6))==0) {
	  he.def_wig=1.0e99;
	} else {
	  he.def_wig=o2scl::stod(stemp.substr(26,6));
	}
	if (o2scl::count_words(stemp.substr(32,9))==0) {
	  he.Sn=1.0e99;
	} else {
	  he.Sn=o2scl::stod(stemp.substr(32,9));
	}
	if (o2scl::count_words(stemp.substr(41,9))==0) {
	  he.Sp=1.0e99;
	} else {
	  he.Sp=o2scl::stod(stemp.substr(41,9));
	}
	if (o2scl::count_words(stemp.substr(50,9))==0) {
	  he.Qbet=1.0e99;
	} else {
	  he.Qbet=o2scl::stod(stemp.substr(50,9));
	}
	if (o2scl::count_words(stemp.substr(59,9))==0) {
	  he.Mcal=1.0e99;
	} else {
	  he.Mcal=o2scl::stod(stemp.substr(59,9));
	}
	if (o2scl::count_words(stemp.substr(68,9))==0) {
	  he.Err=1.0e99;
	} else {
	  he.Err=o2scl::stod(stemp.substr(68,9));
	}
      } else {
	if (o2scl::count_words(stemp.substr(8,8))==0) {
	  he.bet2=1.0e99;
	} else {
	  he.bet2=o2scl::stod(stemp.substr(8,8));
	}
	if (o2scl::count_words(stemp.substr(16,8))==0) {
	  he.bet4=1.0e99;
	} else {
	  he.bet4=o2scl::stod(stemp.substr(16,8));
	}
	if (o2scl::count_words(stemp.substr(24,8))==0) {
	  he.Rch=1.0e99;
	} else {
	  he.Rch=o2scl::stod(stemp.substr(24,8));
	}
	if (o2scl::count_words(stemp.substr(32,8))==0) {
	  he.def_wig=1.0e99;
	} else {
	  he.def_wig=o2scl::stod(stemp.substr(32,8));
	}
	if (o2scl::count_words(stemp.substr(40,8))==0) {
	  he.Sn=1.0e99;
	} else {
	  he.Sn=o2scl::stod(stemp.substr(40,8));
	}
	if (o2scl::count_words(stemp.substr(48,8))==0) {
	  he.Sp=1.0e99;
	} else {
	  he.Sp=o2scl::stod(stemp.substr(48,8));
	}
	if (o2scl::count_words(stemp.substr(56,8))==0) {
	  he.Qbet=1.0e99;
	} else {
	  he.Qbet=o2scl::stod(stemp.substr(56,8));
	}
	if (o2scl::count_words(stemp.substr(64,8))==0) {
	  he.Mcal=1.0e99;
	} else {
	  he.Mcal=o2scl::stod(stemp.substr(64,8));
	}
	if (o2scl::count_words(stemp.substr(72,8))==0) {
	  he.Err=1.0e99;
	} else {
	  he.Err=o2scl::stod(stemp.substr(72,8));
	}
      }
      if (j==0 || j==jmax-1) {
	if (j==0) {
	  cout << "First line: " << endl;
	} else {
	  cout << "Last line: " << endl;
	}
	cout << j << " " << he.Z << " " << he.A << " " << he.N << endl;
	cout << "\t" << he.bet2 << " " << he.bet4 << " " << he.Rch << endl;
	cout << "\t" << he.def_wig << " " << he.Sn << " " << he.Sp << endl;
	cout << "\t" << he.Qbet << " " << he.Mcal << " " << he.Err << endl;
      }
      list.push_back(he);
    }
    
    // Make HDF table
    size_t offset[12]={HOFFSET(nucmass_hfb_entry,N),
		       HOFFSET(nucmass_hfb_entry,Z),
		       HOFFSET(nucmass_hfb_entry,A),
		       HOFFSET(nucmass_hfb_entry,bet2),
		       HOFFSET(nucmass_hfb_entry,bet4),
		       HOFFSET(nucmass_hfb_entry,Rch),
		       HOFFSET(nucmass_hfb_entry,def_wig),
		       HOFFSET(nucmass_hfb_entry,Sn),
		       HOFFSET(nucmass_hfb_entry,Sp),
		       HOFFSET(nucmass_hfb_entry,Qbet),
		       HOFFSET(nucmass_hfb_entry,Mcal),
		       HOFFSET(nucmass_hfb_entry,Err)};
    
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
				 12,list.size(),sizeof(nucmass_hfb_entry),
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

  } else {

    nucmass_hfb_sp_entry he;
    string ssp, ssn;
    int N2, A2;
    for(size_t j=0;j<jmax;j++) {
      getline(fin,stemp);
      he.Z=o2scl::stoi(stemp.substr(0,4));
      he.A=o2scl::stoi(stemp.substr(4,4));
      he.N=he.A-he.Z;
      if (o2scl::count_words(stemp.substr(8,6))==0) {
	he.bet2=1.0e99;
      } else {
	he.bet2=o2scl::stod(stemp.substr(8,6));
      }
      if (o2scl::count_words(stemp.substr(14,6))==0) {
	he.bet4=1.0e99;
      } else {
	he.bet4=o2scl::stod(stemp.substr(14,6));
      }
      if (o2scl::count_words(stemp.substr(20,6))==0) {
	he.Rch=1.0e99;
      } else {
	he.Rch=o2scl::stod(stemp.substr(20,6));
      }
      if (o2scl::count_words(stemp.substr(26,6))==0) {
	he.def_wig=1.0e99;
      } else {
	he.def_wig=o2scl::stod(stemp.substr(26,6));
      }
      if (o2scl::count_words(stemp.substr(32,9))==0) {
	he.Sn=1.0e99;
      } else {
	he.Sn=o2scl::stod(stemp.substr(32,9));
	if (fabs(he.Sn-999.99)<1.0e-6) {
	  he.Sn=1.0e99;
	}
      }
      if (o2scl::count_words(stemp.substr(41,9))==0) {
	he.Sp=1.0e99;
      } else {
	he.Sp=o2scl::stod(stemp.substr(41,9));
	if (fabs(he.Sp-999.99)<1.0e-6) {
	  he.Sp=1.0e99;
	}
      }
      if (o2scl::count_words(stemp.substr(50,9))==0) {
	he.Qbet=1.0e99;
      } else {
	he.Qbet=o2scl::stod(stemp.substr(50,9));
	if (fabs(he.Qbet-999.99)<1.0e-6) {
	  he.Qbet=1.0e99;
	}
      }
      if (o2scl::count_words(stemp.substr(59,9))==0) {
	he.Mcal=1.0e99;
      } else {
	he.Mcal=o2scl::stod(stemp.substr(59,9));
      }
      if (o2scl::count_words(stemp.substr(68,9))==0) {
	he.Err=1.0e99;
      } else {
	he.Err=o2scl::stod(stemp.substr(68,9));
	if (fabs(he.Err-999.99)<1.0e-6) {
	  he.Err=1.0e99;
	}
      }
      if (o2scl::count_words(stemp.substr(80,5))==0) {
	he.Jexp=1.0e99;
      } else {
	he.Jexp=o2scl::stod(stemp.substr(80,5));
	if (fabs(he.Jexp-99.9)<1.0e-6) {
	  he.Jexp=1.0e99;
	}
      }
      if (o2scl::count_words(stemp.substr(85,5))==0) {
	he.Jth=1.0e99;
      } else {
	he.Jth=o2scl::stod(stemp.substr(85,5));
      }
      if (o2scl::count_words(stemp.substr(90,5))==0) {
	he.Pexp=99;
      } else {
	he.Pexp=o2scl::stoi(stemp.substr(90,5));
	if (he.Pexp==9) he.Pexp=99;
      }
      if (o2scl::count_words(stemp.substr(95,5))==0) {
	he.Pth=99;
      } else {
	he.Pth=o2scl::stoi(stemp.substr(95,5));
      }
      if (j==0 || j==jmax-1) {
	if (j==0) {
	  cout << "First line: " << endl;
	} else {
	  cout << "Last line: " << endl;
	}
	cout << j << " " << he.Z << " " << he.A << " " << he.N << endl;
	cout << "\t" << he.bet2 << " " << he.bet4 << " " << he.Rch << endl;
	cout << "\t" << he.def_wig << " " << he.Sn << " " << he.Sp << endl;
	cout << "\t" << he.Qbet << " " << he.Mcal << " " << he.Err << endl;
	cout << "\t" << he.Jexp << " " << he.Jth << " " 
	     << he.Pexp << " " << he.Pth << endl;
      }
      list_sp.push_back(he);
    }
    
    cout << list_sp.size() << endl;
    
    // Make HDF table
    size_t offset[16]={HOFFSET(nucmass_hfb_sp_entry,N),
		       HOFFSET(nucmass_hfb_sp_entry,Z),
		       HOFFSET(nucmass_hfb_sp_entry,A),
		       HOFFSET(nucmass_hfb_sp_entry,bet2),
		       HOFFSET(nucmass_hfb_sp_entry,bet4),
		       HOFFSET(nucmass_hfb_sp_entry,Rch),
		       HOFFSET(nucmass_hfb_sp_entry,def_wig),
		       HOFFSET(nucmass_hfb_sp_entry,Sn),
		       HOFFSET(nucmass_hfb_sp_entry,Sp),
		       HOFFSET(nucmass_hfb_sp_entry,Qbet),
		       HOFFSET(nucmass_hfb_sp_entry,Mcal),
		       HOFFSET(nucmass_hfb_sp_entry,Err),
		       HOFFSET(nucmass_hfb_sp_entry,Jexp),
		       HOFFSET(nucmass_hfb_sp_entry,Jth),
		       HOFFSET(nucmass_hfb_sp_entry,Pexp),
		       HOFFSET(nucmass_hfb_sp_entry,Pth)};
    
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
    hf.seti("nrecords",list_sp.size());
    cout << "nrecords: " << list_sp.size() << endl;
    hf.sets_fixed("comment",
		  ((string)"HDF5 version of HFB ")+
		  "mass data created for O2scl. "
		  "See http://o2scl.sourceforge.net for details.");
    
    herr_t status=H5TBmake_table(orig_file.c_str(),file,outname.c_str(),
				 16,list_sp.size(),sizeof(nucmass_hfb_sp_entry),
				 names,offset,field_type,100,0,0,&list_sp[0]);

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

  }


  return 0;
}
