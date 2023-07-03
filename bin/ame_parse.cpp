/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2023, Andrew W. Steiner

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
/*
  This code generates the O2scl HDF files original AME data tables

  6/3/14 - This code produces the HDF .o2 tables and compares
  the results to the internal .o2 files and found no differences.

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

  Columns in Nubase '16
  1-3     A
  4       Space
  5-9     Z
  10-11   Space
  12-17   Element name 
  18      Space
  19-29   Mass excess
  30-38   Mass excess uncertainty
  39-55   Excitation energy and uncertainty (column 48 is mixed,
  containing excitation energy or uncertainty)
  57-59   Excitation energy origin code
  60      Space
  61-69   Half life
  70-71   Half life unit
  72-77   Half life uncertainty
  80-92   Spin, parity, and isospin multiplet
  93      Space
  94-95   Year of ENSDF archive
  96      Space
  97-105  Reference
  106-109 Year of discovery
  110     Space
  111-186 Decay mode and intensity
*/
#include <iostream>
#include <fstream>
#include <string>

// for exit()
#include <cstdlib>

#include <o2scl/string_conv.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucmass.h>
#include <o2scl/hdf_io.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_nucmass_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

/** \brief Parse strings \c s1 and \c s2 from the AME into a value,
    \c d1, an uncertainty, \c d2, and an accuracy flag, \c acc
    
    - If string \c s1 has an asterisk, then \c d1 and \c d2 are
    set to zero and \c acc is set to \ref nucmass_ame::not_calculable.
    - If string \c s2 contains the letter 'a', then \c d2 is set to
    zero and \c ass is set to \ref nucmass_ame::unc_less_than_half_eV.
    The value of d1 is computed from <tt>stod_nothrow()</tt>.
    - Otherwise, if string \c s1 has a pound sign, then \c acc is set
    to \ref nucmass_ame::estimated, otherwise, \c acc is set to \ref
    nucmass_ame::measured. The values of \c d1 and \c d2 are computed
    from \c s1 and \c s2 using <tt>stod_nothrow()</tt>.
    - If either of the stod_nothrow() calls returns a non-zero value,
    then the error handler is called.
*/
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
  int ret1=o2scl::stod_nothrow(s1,d1);
  if (ret1!=0) {
    cerr << "Failed to convert: '" << s1 << "'." << endl;
    O2SCL_ERR("Failed to convert first string in parse().",
              o2scl::exc_einval);
  }
  if (s2.find('a')!=string::npos) {
    d2=0.0;
    acc=nucmass_ame::unc_less_than_half_eV;
  } else {
    int ret2=o2scl::stod_nothrow(s2,d2);
    if (ret2!=0) {
      cerr << "Failed to convert: '" << s2 << "'." << endl;
      O2SCL_ERR("Failed to convert second string in parse().",
                o2scl::exc_einval);
    }
  }
  return 0;
}

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  static const size_t n_files=11;
  string fnames[n_files]={"ame95/mass_exp.mas95",
                          "ame95/mass_rmd.mas95",
                          "ame03/mass.mas03",
                          "ame03/mass.mas03round",
                          "ame12/mass.mas12",
                          "ame16/mass16.txt",
                          "ame16/mass16round.txt",
                          "ame16/nubase2016.txt",
                          "ame20/mass.mas20.txt",
                          "ame20/massround.mas20.txt",
                          "ame20/nubase2020.txt"};
  string outnames[n_files]={"ame95exp.o2",
                            "ame95rmd.o2",
                            "ame03.o2",
                            "ame03round.o2",
                            "ame12.o2",
                            "ame16.o2",
                            "ame16round.o2",
                            "nubase16.o2",
                            "ame20.o2",
                            "ame20round.o2",
                            "nubase20.o2"};
		      
  if (argc<3) {
    cout << "Usage: ame_parse <dir> <index>, where <dir> is the "
	 << "directory containing the\noriginal data files (inside "
         << "their respective subdirectories and <index> is \n"
         << "the index from 0 to 10.\n" << endl;
    cout << "ix ";
    cout.width(27);
    cout << "input filename" << " ";
    cout.width(15);
    cout << "output file" << endl;
    cout << "-- ";
    cout.width(27);
    cout << "--------------" << " ";
    cout.width(15);
    cout << "-----------" << endl;
    for(size_t i=0;i<n_files;i++) {
      cout.width(2);
      cout << i << " ";
      cout.width(27);
      cout << fnames[i] << " ";
      cout.width(15);
      cout << outnames[i] << endl;
    }
    exit(-1);
  }
  string dir=argv[1];
  int ik=o2scl::stoi(argv[2]);
  
  nucmass_info nmi;
  
  int count=0;

  // Output every 1000 lines
  //const size_t output=1000;
  const size_t output=1;

  /*
    if (true) {

    cout << "--------------------------------------------------------" << endl;

    string fname=dir+"/"+fnames[ik], tmp, tmp2;
    ifstream fin(fname.c_str());
    
    vector<nucmass_nubase::entry> list;

    while (getline(fin,tmp)) {
      
    }
    
    }
  */
  
  //for (size_t ik=0;ik<n_files;ik++) {
    
  cout << "--------------------------------------------------------" << endl;
  
  string fname=dir+"/"+fnames[ik], tmp, tmp2;
  ifstream fin(fname.c_str());
  
  int sret=system((((string)"rm -f ")+outnames[ik]).c_str());
  
  vector<nucmass_ame::entry> list;

  if (ik==8) {
    for(size_t i=0;i<36;i++) getline(fin,tmp);
  } else if (ik==9) {
    for(size_t i=0;i<34;i++) getline(fin,tmp);
  } else {
    for(size_t i=0;i<39;i++) getline(fin,tmp);
  }
  cout << "Filename: " << fname << endl;
  cout << endl;
    
  nucmass_ame::entry ae;
    
  while (getline(fin,tmp)) {
    
    vector<string> entries;
    
    if (ik==0 || ik==1) {
      
      // 1995 format
      parse_fortran_format(tmp,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f11.3,f9.3,f11.3,f9.3,4x,a2,"+
                           "f11.3,f9.3,2x,i3,1x,f10.3,f9.3",entries);
        
    } else if (ik==2 || ik==3) {
        
      // 2003 format
      parse_fortran_format(tmp,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f13.5,f11.5,f11.3,f9.3,1x,a2,"+
                           "f11.3,f9.3,1x,i3,1x,f12.5,f11.3,1x",entries);
        
    } else if (ik==4 || ik==5 || ik==6) {
        
      // 2012 and 2016 format
      parse_fortran_format(tmp,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f13.5,f11.5,f11.3,f9.3,1x,a2,"+
                           "f11.3,f9.3,1x,i3,1x,f12.5,f11.5",entries);
        
        
    } else if (ik==7) {
        
    } else if (ik==8) {
        
      // 2020 experimental format
      parse_fortran_format(tmp,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f14.6,f12.6,f15.5,f11.5,1x,a2,"+
                           "f13.5,f11.5,1x,i3,1x,f13.6,f12.6",entries);
        
    } else if (ik==9) {
        
      // 2020 recommended format
      parse_fortran_format(tmp,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f13.5,f11.5,f11.3,f9.3,1x,a2,"+
                           "f11.3,f9.3,1x,i3,1x,f12.5,f11.5",entries);
        
    } else if (ik==10) {

      // 2020 Nubase format
      parse_fortran_format(tmp,((string)"a3,1x,a4,3x,a5,a1,1x,f13.6,")+
                           "f11.6,f12.6,f11.6,a2,a1,a1,f9.4,"+
                           "a2,1x,a7,a14,a2,10x,a4,1x,a90",entries);
        
    }
      
    if (ik<7 || ik==8 || ik==9) {
        
      // The line feed character is not read

      // Read N and Z first to compute NMZ and A
      ae.N=o2scl::stoi(entries[2]);
      ae.Z=o2scl::stoi(entries[3]);
      if (ik==3 || ik==6 || ik==9) {
        ae.NMZ=ae.N-ae.Z;
        ae.A=ae.N+ae.Z;
      } else {
        ae.NMZ=o2scl::stoi(entries[1]);
        ae.A=o2scl::stoi(entries[4]);
      }
      string_to_char_array(entries[5],ae.el,4);
      string_to_char_array(entries[6],ae.orig,5);
      parse(entries[7],entries[8],ae.mass,ae.dmass,
            ae.mass_acc);
        
      if (ik<2) {
        parse(entries[9],entries[10],ae.be,ae.dbe,
              ae.be_acc);
        // The 1995 files tabulate the binding energy
        // rather than the binding energy per nucleon
        ae.beoa=ae.be/ae.A;
        ae.dbeoa=ae.dbe/ae.A;
        ae.beoa_acc=nucmass_ame::intl_computed;
      } else {
        parse(entries[9],entries[10],ae.beoa,ae.dbeoa,
              ae.beoa_acc);
        ae.be=ae.beoa*ae.A;
        ae.dbe=ae.dbeoa*ae.A;
        ae.be_acc=nucmass_ame::intl_computed;
      }
        
      string_to_char_array(entries[11],ae.bdmode,3);
      parse(entries[12],entries[13],ae.bde,ae.dbde,ae.bde_acc);
	
      ae.A2=o2scl::stoi(entries[14]);
      parse(entries[15],entries[16],ae.amass,ae.damass,
            ae.amass_acc);
          
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
        //char ch;
        //cin >> ch;
      }
      
    } else if (ik==10) {
      
    }

    count++;
  }
  
  cout << "count: " << count << endl;

  // Directly compare the new file to the o2scl version
  if (true) {
    
    cout << "Checking: " << endl;
    nucmass_ame ame_o2scl;
    string stmp2="12";
    if (ik==0) stmp2="95exp";
    else if (ik==1) stmp2="95rmd";
    else if (ik==2) stmp2="03";
    else if (ik==3) stmp2="03round";
    else if (ik==4) stmp2="12";
    else if (ik==5) stmp2="16";
    else if (ik==6) stmp2="16round";
    else if (ik==8) stmp2="20";
    else if (ik==9) stmp2="20round";
    ame_load(ame_o2scl,stmp2,false);
    
    for(size_t i=0;i<list.size();i++) {
      
      if (i%100==0) {
        cout << "Checking " << i << endl;
      }
      
      int Z=list[i].Z;
      int N=list[i].N;
      nucmass_ame::entry e=ame_o2scl.get_ZN(Z,N);
      
      if (e.NMZ!=list[i].NMZ) {
        cout << "Problem NMZ (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.NMZ << " "
             << list[i].NMZ << endl;
      }
      if (e.N!=list[i].N) {
        cout << "Problem N (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.N << " "
             << list[i].N << endl;
      }
      if (e.Z!=list[i].Z) {
        cout << "Problem Z (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.Z << " "
             << list[i].Z << endl;
      }
      if (e.A!=list[i].A) {
        cout << "Problem A (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.A << " "
             << list[i].A << endl;
      }
      if (((string)e.el)!=((string)list[i].el)) {
        cout << "Problem el (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.el << " "
             << list[i].el << endl;
      }
      if (((string)e.orig)!=((string)list[i].orig)) {
        cout << "Problem orig (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.orig << " "
             << list[i].orig << endl;
      }
      if (fabs(e.mass-list[i].mass)/fabs(e.mass)>1.0e-14) {
        cout << "Problem mass (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.mass << " "
             << list[i].mass << endl;
      }
      if (fabs(e.dmass-list[i].dmass)/fabs(e.dmass)>1.0e-14) {
        cout << "Problem dmass (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.dmass << " "
             << list[i].dmass << endl;
      }
      if (e.mass_acc!=list[i].mass_acc) {
        cout << "Problem mass_acc (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.mass_acc << " "
             << list[i].mass_acc << endl;
      }
      if (fabs(e.be-list[i].be)/fabs(e.be)>1.0e-14) {
        cout << "Problem be (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.be << " "
             << list[i].be << endl;
      }
      if (fabs(e.dbe-list[i].dbe)/fabs(e.dbe)>1.0e-14) {
        cout << "Problem dbe (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.dbe << " "
             << list[i].dbe << endl;
      }
      if (e.be_acc!=list[i].be_acc) {
        cout << "Problem be_acc (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.be_acc << " "
             << list[i].be_acc << endl;
      }
      if (fabs(e.beoa-list[i].beoa)/fabs(e.beoa)>1.0e-14) {
        cout << "Problem beoa (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.beoa << " "
             << list[i].beoa << endl;
      }
      if (fabs(e.dbeoa-list[i].dbeoa)/fabs(e.dbeoa)>1.0e-14) {
        cout << "Problem dbeoa (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.dbeoa << " "
             << list[i].dbeoa << endl;
      }
      if (e.beoa_acc!=list[i].beoa_acc) {
        cout << "Problem beoa_acc (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.beoa_acc << " "
             << list[i].beoa_acc << endl;
      }
      if (((string)e.bdmode)!=((string)list[i].bdmode)) {
        cout << "Problem bdmode (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.bdmode << " "
             << list[i].bdmode << endl;
      }
      if (fabs(e.bde-list[i].bde)/fabs(e.bde)>1.0e-14) {
        cout << "Problem bde (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.bde << " "
             << list[i].bde << endl;
      }
      if (fabs(e.dbde-list[i].dbde)/fabs(e.dbde)>1.0e-14) {
        cout << "Problem dbde (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.dbde << " "
             << list[i].dbde << endl;
      }
      if (e.bde_acc!=list[i].bde_acc) {
        cout << "Problem bde_acc (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.bde_acc << " "
             << list[i].bde_acc << endl;
      }
      if (fabs(e.A2-list[i].A2)/fabs(e.A2)>1.0e-14) {
        cout << "Problem A2 (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.A2 << " "
             << list[i].A2 << endl;
      }
      if (fabs(e.amass-list[i].amass)/fabs(e.amass)>1.0e-14) {
        cout << "Problem amass (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.amass << " "
             << list[i].amass << endl;
      }
      if (fabs(e.damass-list[i].damass)/fabs(e.damass)>1.0e-14) {
        cout << "Problem damass (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.damass << " "
             << list[i].damass << endl;
      }
      if (e.amass_acc!=list[i].amass_acc) {
        cout << "Problem amass_acc (i,Z,N): " << i << " " << Z << " "
             << N << endl;
        cout << "o2scl vs. current: " << e.amass_acc << " "
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
      
    const char *names[23]={"N-Z",
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

    // Hack to create an hdf_file object with 'write_access=true'
    hf.open_or_create("temp.o2");
    hf.close();
      
    //hf.open_or_create(outnames[ik]);
    hf.set_current_id(file);
    //hid_t file=hf.get_current_id();
    hf.seti("nrecords",list.size());
    hf.sets_fixed("comment",
                  ((string)"HDF5 version of Atomic Mass Evaluation. ")+
                  "data created for O2scl. "
                  "See https://neutronstars.utk.edu/code/o2scl for details.");

    herr_t status;
    
    status=H5TBmake_table(fnames[ik].c_str(),file,
                          outnames[ik].c_str(),
                          23,list.size(),sizeof(nucmass_ame::entry),
                          names,offset,field_type,100,0,1,&list[0]);

    hf.sets_fixed("orig_file",fnames[ik]);
    if (ik==0) {
      hf.sets_fixed("reference",
                    ((string)"G. Audi and A. H. Wapstra, ")+
                    "Nucl. Phys. A, 595 (1995) 409.");
    } else if (ik==1) {
      hf.sets_fixed("reference",
                    ((string)"G. Audi and A. H. Wapstra, ")+
                    "Nucl. Phys. A, 595 (1995) 409.");
    } else if (ik==2) {
      hf.sets_fixed("reference",
                    ((string)"G. Audi, A. H. Wapstra and C. Thibault, ")+
                    "Nucl. Phys. A, 729 (2003) 337.");
    } else if (ik==3) {
      hf.sets_fixed("reference",
                    ((string)"G. Audi, A. H. Wapstra and C. Thibault, ")+
                    "Nucl. Phys. A, 729 (2003) 337.");
    } else if (ik==4) {
      hf.sets_fixed
        ("reference",((string)"G. Audi, M. Wang, A. H. Wapstra, ")+
         "F. G. Kondev, M. MacCormick, X. Xu, and B. Pfeiffer, "+
         "Chin. Phys. C, 36 (2012) 1287; "+
         "M. Wang, G. Audi, A. H. Wapstra, "+
         "F. G. Kondev, M. MacCormick, X. Xu, and B. Pfeiffer, "+
         "Chin. Phys. C, 36 (2012) 1603.");
    } else if (ik>=5 && ik<=7) {
      hf.sets_fixed
        ("reference",((string)"W. J. Huang, G. Audi, M. Wang ")+
         "F. G. Kondev, S. Naimi, and X. Xu, "
         "Chin. Phys. C, 41 (2017) 030002; "+
         "M. Wang, G. Audi, F. G. Kondev, "+
         "W. J. Huang, , S. Naimi, and X. Xu, "
         "Chin. Phys. C, 41 (2017) 030003.");
    } else if (ik>=8 && ik<=10) {
      hf.sets_fixed
        ("reference",((string)"W. J. Huang, M. Wang ")+
         "F. G. Kondev, G. Audi, S. Naimi, X. Xu, "
         "Chin. Phys. C, 45 (2021) 030002; "+
         "M. Wang, W. J. Huang, F. G. Kondev, "+
         "G. Audi, and S. Naimi, "
         "Chin. Phys. C, 41 (2021) 030003.");
    }
      
    H5Tclose(string_type3);
    H5Tclose(string_type4);
    H5Tclose(string_type5);
    H5Fclose(file);
  }
    
  //}

  if (ik<7) {
    int sret2=system((((string)"h5diff ")+outnames[ik]+
                      " ../data/o2scl/nucmass/"+outnames[ik]).c_str());
  }
  
  return 0;
}
