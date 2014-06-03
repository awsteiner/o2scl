/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2014, Andrew W. Steiner

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
  This code generates the O2scl HDF files original Audi et al. data
  tables.

  Columns in AME '03 and '12 (not indexed by zero)
  1 unused
  2-4 N-Z (3)
  5-9 N (5)
  10-14 Z (5)
  15-19 A (5)
  20 unused
  21-23 element (3)
  24-27 origin (4)
  28 unused
  29-41 mass excess (13)
  42-52 mass excess unc. (11)
  53-63 binding energy (11)
  64-72 binding energy unc. (9)
  73 unused
  74-75 beta decay mode (2)
  76-86 beta-decay energy (11)
  87-95 beta-decay energy unc. (9)
  96 unused
  97-99 A (3)
  100 unused
  101-112 atomic mass (12)
  113-123 atomic mass unc. (11)
*/
#include <iostream>
#include <fstream>
#include <string>

// for exit()
#include <cstdlib>

#include <o2scl/string_conv.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucmass.h>
#include <o2scl/hdf_file.h>

#include <hdf5_hl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

nucmass_ame::entry ae;

int parse(string s1, string s2, double &d1, double &d2, int &acc) {
  if (s1.find('*')!=string::npos) {
    d1=0.0;
    d2=0.0;
    acc=nucmass_ame::not_calculable;
    return 0;
  } 
  if (s1.find('#')!=string::npos) {
    acc=nucmass_ame::estimated;
  } else {
    acc=nucmass_ame::measured;
  }
  d1=stod(s1,false);
  d2=stod(s2,false);
  return 0;
}

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  ofstream fout;

  if (argc<2) {
    cout << "Usage: ame_parse <dir>, where <dir> is the directory\n"
	 << "containing the original Audi et al. data files." << endl;
    exit(-1);
  }
  string dir=argv[1];
  string fnames[5]={"mass_exp.mas95","mass_rmd.mas95","mass.mas03",
		    "mass.mas03round","mass.mas12"};
  string outnames[5]={"ame95exp.o2","ame95rmd.o2","ame03.o2",
		      "ame03round.o2","ame12.o2"};

  nucmass_info nmi;
  
  int count=0;
  const size_t output=1000;

  for (size_t ik=0;ik<5;ik++) {
    
    cout << "--------------------------------------------------------" << endl;

    string fname=dir+"/"+fnames[ik], tmp, tmp2;
    ifstream fin(fname.c_str());

    vector<nucmass_ame::entry> list;
  
    for(size_t i=0;i<39;i++) getline(fin,tmp);
    cout << "Filename: " << fname << endl;
    cout << endl;
    
    while (getline(fin,tmp)) {

      ae.el[0]='\0';
      ae.orig[0]='\0';
      ae.bdmode[0]='\0';

      ae.N=o2scl::stoi(tmp.substr(4,5));
      ae.Z=o2scl::stoi(tmp.substr(9,5));
      if (ik!=3) {
	ae.NMZ=o2scl::stoi(tmp.substr(1,3));
	ae.A=o2scl::stoi(tmp.substr(14,5));
	if (ae.NMZ!=ae.N-ae.Z) {
	  cout << "N-Z not correct. N=" << ae.N << " Z=" << ae.Z << endl;
	  exit(-1);
	}
	if (ae.A!=ae.N+ae.Z) {
	  cout << "N+Z not correct. N=" << ae.N << " Z=" << ae.Z << endl;
	  exit(-1);
	}
      } else {
	ae.NMZ=ae.N-ae.Z;
	ae.A=ae.N+ae.Z;
      }

      if (count%output==0) {
	cout << "Z, N, A, N-Z: ";
	cout.width(3);
	cout << ae.Z << " ";
	cout.width(3);
	cout << ae.N << " ";
	cout.width(3);
	cout << ae.A << " ";
	cout.width(3);
	cout << ae.NMZ << endl;
      }
      
      tmp2=tmp.substr(20,3);
      remove_whitespace(tmp2);
      if (tmp2.length()>0) { ae.el[0]=tmp2[0]; ae.el[1]='\0'; }
      if (tmp2.length()>1) { ae.el[1]=tmp2[1]; ae.el[2]='\0'; }
      if (tmp2.length()>2) { ae.el[1]=tmp2[2]; ae.el[3]='\0'; }
      
      if ((ik<2 && ae.Z<103) || (ik>=2 && ik<4 && ae.Z<110) || 
	  (ik==4 && ae.Z<113)) {
	if (((string)ae.el)!=nmi.Ztoel(ae.Z)) {
	  cout << "Element name incorrect: " << ae.el << " " << ae.Z << endl;
	  exit(-1);
	}
      }

      tmp2=tmp.substr(23,4);
      remove_whitespace(tmp2);
      if (tmp2.length()>0) { ae.orig[0]=tmp2[0]; ae.orig[1]='\0'; }
      if (tmp2.length()>1) { ae.orig[1]=tmp2[1]; ae.orig[2]='\0'; }
      if (tmp2.length()>2) { ae.orig[2]=tmp2[2]; ae.orig[3]='\0'; }
      if (tmp2.length()>3) { ae.orig[3]=tmp2[3]; ae.orig[4]='\0'; }

      if (count%output==0) {
	cout << "el,orig: '" << ae.el << "' '" << ae.orig << "'" << endl;
      }
      
      if (ik<2) {
	
	// 95 format
	parse(tmp.substr(28,11),tmp.substr(39,9),ae.mass,ae.dmass,ae.mass_acc);
	if (count%output==0) {
	  cout << "mass: '" << tmp.substr(28,11) << "' '" << tmp.substr(39,9) 
	       << "' " << ae.mass << " " << ae.dmass << " " 
	       << ae.mass_acc << endl;
	}
	parse(tmp.substr(48,11),tmp.substr(59,9),ae.be,ae.dbe,ae.be_acc);
	if (count%output==0) {
	  cout << "binding: '" << tmp.substr(48,11) << "' '" 
	       << tmp.substr(59,9) 
	       << "' " << ae.beoa << " " << ae.dbeoa << " " 
	       << ae.beoa_acc << endl;
	}
	
	tmp2=tmp.substr(72,2);
	remove_whitespace(tmp2);
	if (tmp2.length()>0) { ae.bdmode[0]=tmp2[0]; ae.bdmode[1]='\0'; }
	if (tmp2.length()>1) { ae.bdmode[1]=tmp2[1]; ae.bdmode[2]='\0'; }
	
	if (count%output==0) {
	  cout << "bdmode: '" << ae.bdmode << "'" << endl;
	}

	parse(tmp.substr(74,11),tmp.substr(85,9),ae.bde,ae.dbde,ae.bde_acc);
	if (count%output==0) {
	  cout << "beta: '" << tmp.substr(74,11) << "' '" 
	       << tmp.substr(85,9) 
	       << "' " << ae.bde << " " << ae.dbde << " " 
	       << ae.bde_acc << endl;
	}

	ae.A2=o2scl::stoi(tmp.substr(96,3));
	parse(tmp.substr(100,10),tmp.substr(110,9),ae.amass,ae.damass,
	      ae.amass_acc);
	if (count%output==0) {
	  cout << "amass: '" << tmp.substr(100,10) << "' '" 
	       << tmp.substr(110,9) 
	       << "' " << ae.amass << " " << ae.damass << " " 
	       << ae.amass_acc << endl;
	  cout << endl;
	}

	ae.beoa=ae.be/ae.A;
	ae.dbeoa=ae.dbe/ae.A;
	
      } else {

	// 03 and 2012 format
	parse(tmp.substr(28,13),tmp.substr(41,11),ae.mass,ae.dmass,
	      ae.mass_acc);
	if (count%output==0) {
	  cout << "mass: '" << tmp.substr(28,13) << "' '" << tmp.substr(41,11) 
	       << "' " << ae.mass << " " << ae.dmass << " " 
	       << ae.mass_acc << endl;
	}
	parse(tmp.substr(52,11),tmp.substr(63,9),ae.beoa,ae.dbeoa,ae.beoa_acc);
	if (count%output==0) {
	  cout << "binding: '" << tmp.substr(52,11) << "' '" 
	       << tmp.substr(63,9) 
	       << "' " << ae.beoa << " " << ae.dbeoa << " " 
	       << ae.beoa_acc << endl;
	}
	
	tmp2=tmp.substr(73,2);
	remove_whitespace(tmp2);
	if (tmp2.length()>0) { ae.bdmode[0]=tmp2[0]; ae.bdmode[1]='\0'; }
	if (tmp2.length()>1) { ae.bdmode[1]=tmp2[1]; ae.bdmode[2]='\0'; }
	
	if (count%output==0) {
	  cout << "bdmode: '" << ae.bdmode << "'" << endl;
	}

	parse(tmp.substr(75,11),tmp.substr(86,9),ae.bde,ae.dbde,ae.bde_acc);
	if (count%output==0) {
	  cout << "beta: '" << tmp.substr(75,11) << "' '" 
	       << tmp.substr(86,9) 
	       << "' " << ae.bde << " " << ae.dbde << " " 
	       << ae.bde_acc << endl;
	}
	
	ae.A2=o2scl::stoi(tmp.substr(96,3));
	parse(tmp.substr(100,12),tmp.substr(112,11),ae.amass,ae.damass,
	      ae.amass_acc);
	if (count%output==0) {
	  cout << "amass: '" << tmp.substr(100,12) << "' '" 
	       << tmp.substr(112,11) 
	       << "' " << ae.amass << " " << ae.damass << " " 
	       << ae.amass_acc << endl;
	  cout << endl;
	}

	ae.be=ae.beoa*ae.A;
	ae.dbe=ae.dbeoa*ae.A;

      }

      list.push_back(ae);

      if (false) {
	cout << ae.NMZ << ", ," << ae.N << ", ," << ae.Z << ", ," 
	     << ae.A << ", ," << ae.el << ", ," << ae.orig << ", ," 
	     << ae.mass << ", ," << ae.dmass << ", ," 
	     << ae.mass_acc << ", ," << ae.beoa << ", ,"
	     << ae.dbeoa << ", ," << ae.beoa_acc << ", ," 
	     << ae.bdmode << ", ," << ae.bde << ", ," 
	     << ae.dbde << ", ," << ae.bde_acc << ", ,"
	     << ae.amass << ", ," << ae.damass << ", ," 
	     << ae.amass_acc << endl;
	cout << endl;
	//char ch;
	//cin >> ch;
      }

      count++;
    }
    cout << "count: " << count << endl;

    if (true) {
      cout << "Checking: " << endl;
      nucmass_ame ame12;
      string stmp2="12";
      if (ik==0) stmp2="95rmd";
      else if (ik==1) stmp2="95exp";
      else if (ik==2) stmp2="03";
      else if (ik==3) stmp2="03round";
      ame_load(ame12,stmp2);
      for(size_t i=0;i<list.size();i++) {
	if (i%100==0) {
	  cout << "Checking " << i << endl;
	}
	int Z=list[i].Z;
	int N=list[i].N;
	nucmass_ame::entry e=ame12.get_ZN(Z,N);
	if (e.NMZ!=list[i].NMZ || e.N!=list[i].N ||
	    e.Z!=list[i].Z || e.A!=list[i].A ||
	    ((string)e.el)!=((string)list[i].el) ||
	    ((string)e.orig)!=((string)list[i].orig) ||
	    fabs(e.mass-list[i].mass)/fabs(e.mass)>1.0e-14 ||
	    fabs(e.dmass-list[i].dmass)/fabs(e.dmass)>1.0e-14 ||
	    e.mass_acc!=list[i].mass_acc ||
	    fabs(e.be-list[i].be)/fabs(e.be)>1.0e-14 ||
	    fabs(e.dbe-list[i].dbe)/fabs(e.dbe)>1.0e-14 ||
	    e.be_acc!=list[i].be_acc ||
	    fabs(e.beoa-list[i].beoa)/fabs(e.beoa)>1.0e-14 ||
	    fabs(e.dbeoa-list[i].dbeoa)/fabs(e.dbeoa)>1.0e-14 ||
	    e.beoa_acc!=list[i].beoa_acc ||
	    ((string)e.bdmode)!=((string)list[i].bdmode) ||
	    fabs(e.bde-list[i].bde)/fabs(e.bde)>1.0e-14 ||
	    fabs(e.dbde-list[i].dbde)/fabs(e.dbde)>1.0e-14 ||
	    e.bde_acc!=list[i].bde_acc ||
	    fabs(e.A2-list[i].A2)/fabs(e.A2)>1.0e-14 ||
	    fabs(e.amass-list[i].amass)/fabs(e.amass)>1.0e-14 ||
	    fabs(e.damass-list[i].damass)/fabs(e.damass)>1.0e-14 ||
	    e.amass_acc!=list[i].amass_acc) {
	  cout << "Problem: " << i << " " << Z << " " << N << endl;
	  cout << e.NMZ << " "
	       << e.N << " "
	       << e.Z << " "
	       << e.A << " "
	       << e.el << " "
	       << e.orig << " "
	       << e.mass << " "
	       << e.dmass << " "
	       << e.mass_acc << " "
	       << e.be << " "
	       << e.dbe << " "
	       << e.be_acc << " "
	       << e.beoa << " "
	       << e.dbeoa << " "
	       << e.beoa_acc << " "
	       << e.bdmode << " "
	       << e.bde << " "
	       << e.dbde << " "
	       << e.bde_acc << " "
	       << e.A2 << " "
	       << e.amass << " "
	       << e.damass << " "
	       << e.amass_acc << endl;
	  cout << list[i].NMZ << " "
	       << list[i].N << " "
	       << list[i].Z << " "
	       << list[i].A << " "
	       << list[i].el << " "
	       << list[i].orig << " "
	       << list[i].mass << " "
	       << list[i].dmass << " "
	       << list[i].mass_acc << " "
	       << list[i].be << " "
	       << list[i].dbe << " "
	       << list[i].be_acc << " "
	       << list[i].beoa << " "
	       << list[i].dbeoa << " "
	       << list[i].beoa_acc << " "
	       << list[i].bdmode << " "
	       << list[i].bde << " "
	       << list[i].dbde << " "
	       << list[i].bde_acc << " "
	       << list[i].A2 << " "
	       << list[i].amass << " "
	       << list[i].damass << " "
	       << list[i].amass_acc << endl;
	}
      }
      
    }

    // Make HDF table
    {
      size_t offset[23]={HOFFSET(nucmass_ame::entry,NMZ),
			 HOFFSET(nucmass_ame::entry,N),
			 HOFFSET(nucmass_ame::entry,Z),
			 HOFFSET(nucmass_ame::entry,A),
			 HOFFSET(nucmass_ame::entry,el),
			 HOFFSET(nucmass_ame::entry,orig),
			 HOFFSET(nucmass_ame::entry,mass),
			 HOFFSET(nucmass_ame::entry,dmass),
			 HOFFSET(nucmass_ame::entry,mass_acc),
			 HOFFSET(nucmass_ame::entry,be),
			 HOFFSET(nucmass_ame::entry,dbe),
			 HOFFSET(nucmass_ame::entry,be_acc),
			 HOFFSET(nucmass_ame::entry,beoa),
			 HOFFSET(nucmass_ame::entry,dbeoa),
			 HOFFSET(nucmass_ame::entry,beoa_acc),
			 HOFFSET(nucmass_ame::entry,bdmode),
			 HOFFSET(nucmass_ame::entry,bde),
			 HOFFSET(nucmass_ame::entry,dbde),
			 HOFFSET(nucmass_ame::entry,bde_acc),
			 HOFFSET(nucmass_ame::entry,A2),
			 HOFFSET(nucmass_ame::entry,amass),
			 HOFFSET(nucmass_ame::entry,damass),
			 HOFFSET(nucmass_ame::entry,amass_acc)};
      size_t sizes[23]={sizeof(ae.NMZ),
			sizeof(ae.N),
			sizeof(ae.Z),
			sizeof(ae.A),
			sizeof(ae.el),
			sizeof(ae.orig),
			sizeof(ae.mass),
			sizeof(ae.dmass),
			sizeof(ae.mass_acc),
			sizeof(ae.be),
			sizeof(ae.dbe),
			sizeof(ae.be_acc),
			sizeof(ae.beoa),
			sizeof(ae.dbeoa),
			sizeof(ae.beoa_acc),
			sizeof(ae.bdmode),
			sizeof(ae.bde),
			sizeof(ae.dbde),
			sizeof(ae.bde_acc),
			sizeof(ae.A2),
			sizeof(ae.amass),
			sizeof(ae.damass),
			sizeof(ae.amass_acc)};
      
      const char *names[23]={
	"N-Z",
	"Neutron number",
	"Proton number",
	"Mass number",
	"Element name",
	"Data origin",
	"Mass excess (in keV)",
	"Mass excess uncertainty (in keV)",
	"Mass excess accuracy flag",
	"Binding energy (in keV)",
	"Binding energy uncertainty (in keV)",
	"Binding energy accuracy flag",
	"Binding energy / A (in keV)",
	"Binding energy / A uncertainty (in keV)",
	"Binding energy / A accuracy flag",
	"Beta decay mode",
	"Beta-decay energy (in keV)",
	"Beta-decay energy uncertainty (in keV)",
	"Beta-decay energy accuracy flag",
	"A2",
	"Atomic mass (in keV)",
	"Atomic mass uncertainty (in keV)",
	"Atomic mass accuracy flag"};
      
      // Three character string type
      hid_t string_type3=H5Tcopy(H5T_C_S1);
      H5Tset_size(string_type3,3);
      // Four character string type
      hid_t string_type4=H5Tcopy(H5T_C_S1);
      H5Tset_size(string_type4,4);
      // Five character string type
      hid_t string_type5=H5Tcopy(H5T_C_S1);
      H5Tset_size(string_type5,5);

      hid_t field_type[23]={H5T_NATIVE_INT,H5T_NATIVE_INT,H5T_NATIVE_INT,
			    H5T_NATIVE_INT,string_type4,string_type5,
			    H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			    H5T_NATIVE_INT,H5T_NATIVE_DOUBLE,
			    H5T_NATIVE_DOUBLE,H5T_NATIVE_INT,
			    H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
			    H5T_NATIVE_INT,string_type3,H5T_NATIVE_DOUBLE,
			    H5T_NATIVE_DOUBLE,H5T_NATIVE_INT,
			    H5T_NATIVE_INT,H5T_NATIVE_DOUBLE,
			    H5T_NATIVE_DOUBLE,H5T_NATIVE_INT};

      hid_t file=H5Fcreate(outnames[ik].c_str(),H5F_ACC_TRUNC,
			   H5P_DEFAULT,H5P_DEFAULT);
      
      hdf_file hf;
      hf.set_current_id(file);
      hf.seti("nrecords",list.size());
      hf.sets_fixed("comment",
		    ((string)"HDF5 version of Audi, et al. ")+
		    "mass data created for O2scl. "
		    "See http://o2scl.sourceforge.net for details.");

      herr_t status;
      status=H5TBmake_table(fnames[ik].c_str(),file,
			    outnames[ik].c_str(),
			    23,list.size(),sizeof(nucmass_ame::entry),
			    names,offset,field_type,100,0,1,&list[0]);

      if (ik==0) {
	hf.sets("orig_file","mass_exp.mas95");
	hf.sets_fixed("reference",
		((string)"G. Audi and A. H. Wapstra, ")+
		"Nucl. Phys. A, 595 (1995) 409.");
      } else if (ik==1) {
	hf.sets("orig_file","mass_rmd.mas95");
	hf.sets_fixed("reference",
		((string)"G. Audi and A. H. Wapstra, ")+
		"Nucl. Phys. A, 595 (1995) 409.");
      } else if (ik==2) {
	hf.sets("orig_file","mass.mas03");
	hf.sets_fixed("reference",
		((string)"G. Audi, A. H. Wapstra and C. Thibault, ")+
		"Nucl. Phys. A, 729 (2003) 337.");
      } else if (ik==3) {
	hf.sets("orig_file","mass.mas03round");
	hf.sets_fixed("reference",
		      ((string)"G. Audi, A. H. Wapstra and C. Thibault, ")+
		      "Nucl. Phys. A, 729 (2003) 337.");
      } else {
	hf.sets("orig_file","mass.mas12");
	hf.sets_fixed
	  ("reference",((string)"G. Audi, M. Wang, A. H. Wapstra, ")+
	   "F. G. Kondev, M. MacCormick, X. Xu, and B. Pfeiffer, "+
	   "Chin. Phys. C, 36 (2012) 1287; "+
	   "M. Wang, G. Audi, A. H. Wapstra, "+
	   "F. G. Kondev, M. MacCormick, X. Xu, and B. Pfeiffer, "+
	   "Chin. Phys. C, 36 (2012) 1603.");
      }
      
      H5Tclose(string_type3);
      H5Tclose(string_type4);
      H5Tclose(string_type5);
      H5Fclose(file);
    }
    
  }
  
  return 0;
}
