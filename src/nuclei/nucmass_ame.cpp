/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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
#include <o2scl/nucmass.h>
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/hdf_io.h>
#include <o2scl/cloud_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

nucmass_ame::nucmass_ame() {
  reference="";
  last=0;
}

nucmass_ame::~nucmass_ame() {
}

double nucmass_ame::mass_excess(int Z, int N) {
  entry ret;
  ret=get_ZN(Z,N);
  if (ret.Z==0 && ret.N==0) return 0.0;
  return ret.mass/1.0e3;
}

bool nucmass_ame::is_included(int l_Z, int l_N) {

  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::is_included().",
		  exc_einval);
  }

  if (false) {

    int l_A=l_Z+l_N, mid=last, lo=0, hi=0;

    // Binary search for the correct A first
    if (mass[mid].Z+mass[mid].N!=l_A) {
      if (mass[mid].Z+mass[mid].N>l_A) {
	lo=0;
	hi=mid;
      } else {
	lo=mid;
	hi=n-1;
      }
      while (hi>lo+1) {
	int mp=(lo+hi)/2;
	if (mass[mid].Z+mass[mid].N<l_A) {
	  lo=mp;
	} else {
	  hi=mp;
	}
      }
      mid=lo;
      if (mass[mid].Z+mass[mid].N!=l_A) {
	mid=hi;
      }
      if (mass[mid].Z+mass[mid].N!=l_A) {
	return false;
      }
    }

    // The cached point was the right one, so we're done
    if (mass[mid].N==l_N) {
      if (fabs(mass[mid].mass)>1.0e-20 &&
	  fabs(mass[mid].mass)<1.0e90) {
	return true;
      } else {
	return false;
      }
    }

    // Now look for the right N among all the A's
    while (mass[mid].Z+mass[mid].N==l_A) {
      if (mass[mid].N==l_N) {
	if (fabs(mass[mid].mass)>1.0e-20 &&
	    fabs(mass[mid].mass)<1.0e90) {
	  return true;
	} else {
	  return false;
	}
      } else if (mass[mid].N>l_N) {
	if (mid==0) return false;
	mid--;
      } else {
	if (mid==((int)(n-1))) return false;
	mid++;
      }
    }

    return false;
  }

  for(size_t i=0;i<n;i++) {
    if (mass[i].Z==l_Z && mass[i].N==l_N) {
      return true;
    }
  }
  return false;
}

/*
bool nucmass_ame_exp::is_included(int l_Z, int l_N) {
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame_exp::is_included().",
	      exc_einval);
  }
  for(size_t i=0;i<n;i++) {
    if (mass[i].Z==l_Z && mass[i].N==l_N && mass[i].mass_acc==0) {
      return true;
    }
  }
  return false;
}
*/

nucmass_ame::entry nucmass_ame::get_ZN(int l_Z, int l_N) {
  nucmass_ame::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::get_ZN().",
	      exc_einval);
    return ret;
  }
  for(size_t i=0;i<n;i++) {
    if (mass[i].Z==l_Z && mass[i].N==l_N) {
      ret=mass[i];
    }
  }
  return ret;
}

nucmass_ame::entry nucmass_ame::get_ZA(int l_Z, int l_A) {
  nucmass_ame::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::get_ZA().",
	      exc_einval);
    return ret;
  }
  for(size_t i=0;i<n;i++) {
    if (mass[i].Z==l_Z && mass[i].A==l_A) {
      ret=mass[i];
    }
  }
  return ret;
}

nucmass_ame::entry nucmass_ame::get_elA(string l_el, int l_A) {
  nucmass_ame::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::get_elA().",
	      exc_einval);
    return ret;
  }
  for(size_t i=0;i<n;i++) {
    if (mass[i].el==l_el && mass[i].A==l_A) {
      ret=mass[i];
    }
  }
  return ret;
}

nucmass_ame::entry nucmass_ame::get(string nucleus) {
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame::get().",
	      exc_einval);
    nucmass_ame::entry ret;
    ret.Z=0;
    ret.A=0;
    ret.A2=0;
    ret.N=0;
    ret.NMZ=0;
    ret.el[0]='\0';
    ret.orig[0]='\0';
    ret.mass=0.0;
    ret.dmass=0.0;
    ret.amass=0.0;
    ret.damass=0.0;
    ret.bde=0.0;
    ret.dbde=0.0;
    ret.be=0.0;
    ret.dbe=0.0;
    ret.beoa=0.0;
    ret.dbeoa=0.0;
    ret.bde_acc=0.0;
    ret.mass_acc=0;
    ret.be_acc=0;
    ret.beoa_acc=0;
    ret.amass_acc=0;
    ret.bdmode[0]='\0';
    return ret;
  }
  if (isalpha(nucleus[1])) {
    string el=nucleus.substr(0,2);
    int A=o2scl::stoi(nucleus.substr(2,nucleus.length()-2));
    return get_elA(el,A);
  } 
  string el=nucleus.substr(0,1);
  int A=o2scl::stoi(nucleus.substr(1,nucleus.length()-1));
  return get_elA(el,A);
}

void nucmass_ame::load_ext(std::string name, std::string filename,
                            std::string nubase_file, bool exp_only,
                            int verbose) {

  nucmass_ame::entry ae;

  ifstream fin(filename.c_str());
  std::string line;

  if (verbose>0) {
    std::cout << "nucmass_ame()::load_ext(): "
              << "name,filename,nubase_file: " << name << " "
              << filename << " " << nubase_file << std::endl;
  }
  
  int count=0;

  mass.clear();

  if (name=="20") {
    for(size_t i=0;i<36;i++) getline(fin,line);
  } else if (name=="20round") {
    for(size_t i=0;i<34;i++) getline(fin,line);
  } else {
    for(size_t i=0;i<39;i++) getline(fin,line);
  }
  
  while (getline(fin,line)) {
    
    if (verbose>0 && count%100==0) {
      std::cout << "mass count: " << count << std::endl;
    }

    vector<string> entries;
    
    if (name=="20") {
      // AWS, 1/26/25: this is different from what's quoted in the
      // file, I think it may be a typo in the file.
      parse_fortran_format(line,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f14.6,f12.6,f13.5,f11.5,1x,a2,"+
                           "f13.5,f11.5,1x,i3,1x,f13.6,f12.6",entries);
    } else if (name=="20round") {
      parse_fortran_format(line,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f13.5,f11.5,f11.3,f9.3,1x,a2,"+
                           "f11.3,f9.3,1x,i3,1x,f12.5,f11.5",entries);
    } else if (name=="12" || name=="16round" || name=="16") {
      parse_fortran_format(line,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f13.5,f11.5,f11.3,f9.3,1x,a2,"+
                           "f11.3,f9.3,1x,i3,1x,f12.5,f11.5",entries);
    } else if (name=="03round" || name=="03") {
      parse_fortran_format(line,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f13.5,f11.5,f11.3,f9.3,1x,a2,"+
                           "f11.3,f9.3,1x,i3,1x,f12.5,f11.3,1x",entries);
    } else if (name=="95exp" || name=="95rmd") {
      parse_fortran_format(line,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f11.3,f9.3,f11.3,f9.3,4x,a2,"+
                           "f11.3,f9.3,2x,i3,1x,f10.3,f9.3",entries);
    } else {
      O2SCL_ERR("X.",o2scl::exc_efailed);
    }

    /*
      if (verbose>2) {
      cout << "Entries:" << endl;
      for(size_t k=0;k<entries.size();k++) {
      cout << k << " " << entries[k].size() << " '"
      << entries[k] << "'." << endl;
      }
      }
    */
      
    // Read N and Z first to compute NMZ and A
    ae.N=o2scl::stoi(entries[2]);
    ae.Z=o2scl::stoi(entries[3]);
    
    if (name=="20round" || name=="16round" || name=="03round") {
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
    
    if (name=="95rmd" || name=="95exp") {
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
    
    /*if (ae.N==126 && ae.Z==82) {
      cout << ae.N << " " << ae.Z << " x"
           << entries[10] << "x x" << entries[11] << "x x" 
           << entries[12] << "x" << endl;
      cout << ae.bdmode[0] << " " << ae.bdmode[1] << endl;
      }*/
    
    ae.A2=o2scl::stoi(entries[14]);
    parse(entries[15],entries[16],ae.amass,ae.damass,
          ae.amass_acc);
    
    mass.push_back(ae);
      
    count++;
    
  }

  n=count;

  last=mass.size()/2;
  
  if (name=="95exp" || name=="95rmd") {
    reference=((string)"G. Audi and A. H. Wapstra, ")+
      "Nucl. Phys. A, 595 (1995) 409.";
  } else if (name=="03" || name=="03round") {
    reference=((string)"G. Audi, A. H. Wapstra and C. Thibault, ")+
      "Nucl. Phys. A, 729 (2003) 337.";
  } else if (name=="12") {
    reference=((string)"G. Audi, M. Wang, A. H. Wapstra, ")+
      "F. G. Kondev, M. MacCormick, X. Xu, and B. Pfeiffer, "+
      "Chin. Phys. C, 36 (2012) 1287; "+
      "M. Wang, G. Audi, A. H. Wapstra, "+
      "F. G. Kondev, M. MacCormick, X. Xu, and B. Pfeiffer, "+
      "Chin. Phys. C, 36 (2012) 1603.";
  } else if (name=="16" || name=="16round") {
    reference=((string)"W. J. Huang, G. Audi, M. Wang ")+
      "F. G. Kondev, S. Naimi, and X. Xu, "
      "Chin. Phys. C, 41 (2017) 030002; "+
      "M. Wang, G. Audi, F. G. Kondev, "+
      "W. J. Huang, , S. Naimi, and X. Xu, "
      "Chin. Phys. C, 41 (2017) 030003.";
  } else if (name=="20" || name=="20round") {
    reference=((string)"W. J. Huang, M. Wang ")+
      "F. G. Kondev, G. Audi, S. Naimi, X. Xu, "
      "Chin. Phys. C, 45 (2021) 030002; "+
      "M. Wang, W. J. Huang, F. G. Kondev, "+
      "G. Audi, and S. Naimi, "
      "Chin. Phys. C, 41 (2021) 030003.";
  }

  fin.close();

  if (name=="20") {
    
    fin.open(nubase_file.c_str());
    
    count=0;

    if (name=="20") {
      for(size_t i=0;i<25;i++) getline(fin,line);
    }
    
    while (getline(fin,line)) {
      
      if (verbose>0 && count%100==0) {
        std::cout << "nubase count: " << count << std::endl;
      }

      vector<string> entries;

      parse_fortran_format(line,((string)"a3,1x,a4,3x,a5,a1,1x,f13.6,")+
                           "f11.6,f12.6,f11.6,a2,a1,a1,f9.4,"+
                           "a2,1x,a7,a14,a2,10x,a4,1x,a90",entries);

      if (verbose>2) {
        cout << "Entries:" << endl;
        for(size_t k=0;k<entries.size();k++) {
          cout << k << " '" << entries[k] << "'." << endl;
        }
      }

      int A=o2scl::stoi(entries[0]);
      int Z;
      if (entries[1].length()>3) {
        Z=o2scl::stoi(entries[1].substr(0,3));
      } else {
        Z=o2scl::stoi(entries[1]);
      }

      int ix=0;
      bool found=false;
      for(size_t i=0;i<n && found==false;i++) {
        if (mass[i].Z==Z && mass[i].A==A) {
          found=true;
          ix=i;
        }
      }
      if (found==false) {
        O2SCL_ERR("Couldn't find nubase nucleus.",o2scl::exc_efailed);
      }

      /*
      if (Z==3 && A==7) {
        cout << "Entries:" << endl;
        for(size_t k=0;k<entries.size();k++) {
          cout << k << " '" << entries[k] << "'." << endl;
        }
        //char ch;
        //cin >> ch;
      }
      */
      
      entry_nubase_20 enu20;
      if (entries[1].length()>3 && entries[1][3]!=' ') {
        enu20.Znote=entries[1][3];
      } else {
        enu20.Znote='\0';
      }
      if (verbose>1) cout << "Znote: " << ((int)enu20.Znote) << endl;
      string_to_char_array(entries[2],enu20.A_el,6);
      if (verbose>1) cout << "A_el: '" << enu20.A_el << "'" << endl;
      if (entries[3].length()>0 && entries[3][0]!=' ') {
        enu20.isomer=entries[3][0];
      } else {
        enu20.isomer='\0';
      }
      if (verbose>1) cout << "isomer: '" << enu20.isomer << "'" << endl;
      parse(entries[4],entries[5],enu20.mass,enu20.dmass,
            enu20.mass_acc);
      if (verbose>1) cout << "mass: " << enu20.mass << " " << enu20.dmass << " "
           << enu20.mass_acc << endl;
      remove_whitespace(entries[6]);
      if (entries[6]=="non-exist") {
        enu20.exc_energy=0.0;
        enu20.dexc_energy=0.0;
        enu20.exc_energy_acc=nucmass_ame::does_not_exist;
      } else {
        parse(entries[6],entries[7],enu20.exc_energy,enu20.dexc_energy,
              enu20.exc_energy_acc);
      }
      if (verbose>1) cout << "exc_energy: " << enu20.exc_energy << " "
                          << enu20.dexc_energy << " "
                          << enu20.exc_energy_acc << endl;
      if (entries.size()>8) {
        string_to_char_array(entries[8],enu20.origin,3);
      } else {
        enu20.origin[0]='\0';
      }
      if (verbose>1) cout << "origin: " << enu20.origin << endl;
      if (entries.size()>9 && entries[9].length()>0 && entries[9][0]!=' ') {
        enu20.isomer_unc=entries[9][0];
      } else {
        enu20.isomer_unc='\0';
      }
      if (verbose>1) cout << "isomer_unc: '"
                          << enu20.isomer_unc << "'" << endl;
      if (entries.size()>10 && entries[10].length()>0 &&
          entries[10][0]!=' ') {
        enu20.isomer_inv=entries[10][0];
      } else {
        enu20.isomer_inv='\0';
      }
      if (verbose>1) cout << "isomer_inv: '" << enu20.isomer_inv
                          << "'" << endl;
      if (entries.size()>11 && count_words(entries[11])>0) {
        remove_whitespace(entries[11]);
        if (entries[11]=="stbl") {
          enu20.hlife=0.0;
          enu20.hlife_acc=nucmass_ame::unstable;
        } else if (entries[11]=="p-unst") {
          enu20.hlife=0.0;
          enu20.hlife_acc=nucmass_ame::part_unstable;
        } else if (entries[11][0]=='>') {
          entries[11]=entries[11].substr
            (1,entries[11].length()-1);
          enu20.hlife=o2scl::stod(entries[11]);
          enu20.hlife_acc=nucmass_ame::lower_limit;
        } else if (entries[11][0]=='<') {
          entries[11]=entries[11].substr
            (1,entries[11].length()-1);
          enu20.hlife=o2scl::stod(entries[11]);
          enu20.hlife_acc=nucmass_ame::upper_limit;
        } else if (entries[11][0]=='~') {
          entries[11]=entries[11].substr
            (1,entries[11].length()-1);
          enu20.hlife=o2scl::stod(entries[11]);
          enu20.hlife_acc=nucmass_ame::approximate;
        } else {
          enu20.hlife=o2scl::stod(entries[11]);
          enu20.hlife_acc=0;
        }
      } else {
        enu20.hlife=0.0;
        enu20.hlife_acc=nucmass_ame::blank;
      }
      if (verbose>1) {
        cout << "hlife: " << enu20.hlife << " "
             << enu20.hlife_acc << endl;
      }
      if (entries.size()>12) {
        string_to_char_array(entries[12],enu20.hl_unit,3);
      } else {
        enu20.hl_unit[0]='\0';
      }
      if (verbose>1) cout << "hl_unit: " << enu20.hl_unit << endl;
      if (entries.size()>13) {
        string_to_char_array(entries[13],enu20.dhlife,8);
      } else {
        enu20.dhlife[0]='\0';
      }
      if (verbose>1) cout << "dhlife: " << enu20.dhlife << endl;
      if (entries.size()>14) {
        string_to_char_array(entries[14],enu20.spinp,15);
      } else {
        enu20.spinp[0]='\0';
      }
      if (verbose>1) cout << "spinp: " << enu20.spinp << endl;
      if (entries.size()>15 && count_words(entries[15])>0) {
        enu20.ENSDF_year=o2scl::stoi(entries[15]);
      } else {
        enu20.ENSDF_year=0;
      }
      if (verbose>1) {
        cout << "ENSDF_year: " << enu20.ENSDF_year << endl;
        cout << entries.size() << endl;
      }
      if (entries.size()>16 && count_words(entries[16])>0) {
        enu20.discovery=o2scl::stoi(entries[16]);
      } else {
        enu20.discovery=0;
      }
      if (verbose>1) cout << "discovery: " << enu20.discovery << endl;
      if (entries.size()>17) {
        if (entries[17].size()>90) {
          O2SCL_ERR("Decay intensity field too long for storage.",
                    o2scl::exc_einval);
        }
        string_to_char_array(entries[17],enu20.decay_intensity,91);
      } else {
        enu20.decay_intensity[0]='\0';
      }
      if (verbose>1) {
        cout << "decay_intensity: "
             << enu20.decay_intensity << endl;
      }

      mass[ix].props.push_back(enu20);

      count++;
    }

    fin.close();

    //for(size_t i=0;i<n;i++) {
    //cout << i << " " << mass[i].props.size() << endl;
    //}
    
  }
  
  return;
}

void nucmass_ame::load(std::string name, bool exp_only,
                        int verbose) {

  char *ed=getenv("O2SCL_EXT_DATA");
  std::string ext_data;
  if (ed) {
    ext_data=ed;
  } else {
    ext_data=".";
  }
    
  std::string prefix=o2scl::o2scl_settings.get_data_dir()+"/nucmass";

  cloud_file cf;
  cf.verbose=verbose;
    
  std::string filename, nubase_file;
  
  if (name=="20") {

    // These files don't need hashes because they're included in the repo
    filename=prefix+"/ame20/mass20.txt";
    nubase_file=prefix+"/ame20/nubase_4.mas20.txt";
    
  } else if (name=="20round") {
    
    // These files don't need hashes because they're included in the repo
    filename=prefix+"/ame20/mass20round.txt";
    nubase_file=prefix+"/ame20/nubase_4.mas20.txt";
    
  } else if (name=="16") {
    
    std::string sha=((std::string)"2167f57a2a98331e4649b2dd2b658a")+
      "9006ed4fba1975729ebfe52a42b4b9218a";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("mass16.txt",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame16/mass16.txt",sha,ext_data);

    sha=((std::string)"f3d08e4af75892ec4626805ca3465b")+
      "7925144d53e0bfee713f664c2abd4dd7c4";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("nubase2016.txt",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame16/nubase2016.txt",sha,ext_data);
    
    filename=ext_data+"/mass16.txt";
    nubase_file=ext_data+"/nubase2016.txt";
    
  } else if (name=="16round") {
    
    std::string sha=((std::string)"d34bc538daa65810aede5a7037e7712")+
      "6d5218b319db4b2ad08694b783bc4f249";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("mass16round.txt",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame16/mass16round.txt",sha,ext_data);

    sha=((std::string)"f3d08e4af75892ec4626805ca3465b")+
      "7925144d53e0bfee713f664c2abd4dd7c4";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("nubase2016.txt",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame16/nubase2016.txt",sha,ext_data);
    
    filename=ext_data+"/mass16round.txt";
    nubase_file=ext_data+"/nubase2016.txt";
    
  } else if (name=="12") {
    
    std::string sha=((std::string)"81e887c71c2c54c76caea36fd861b")+
      "195a7f3eeb77d04b520e05fa97e0eedd7f3";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("mass.mas12",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame12/mass.mas12",sha,ext_data);
    
    sha=((std::string)"d69cac4f34e01e5d92ac2c415492c9ee05de")+
      "2ca9b11e6cb3e71786ba66c8679c";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("nubase.mas12",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame12/nubase.mas12",sha,ext_data);
    
    filename=ext_data+"/mass.mas12";
    nubase_file=ext_data+"/nubase.mas12";
    
  } else if (name=="03") {
    
    std::string sha=((std::string)"33405560376f2adfb190beec44213")+
      "523ec79149804df94e436d608019a4c70d1";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("mass.mas03",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame03/mass.mas03",sha,ext_data);

    filename=ext_data+"/mass.mas03";
    nubase_file="";
    
  } else if (name=="03round") {
    
    std::string sha=((std::string)"1e951122a0c2531f14ca7f45343c")+
      "459f4b6d78353298af7c8d8b92fe58ecf403";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("mass.mas03round",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame03/mass.mas03round",sha,ext_data);

    filename=ext_data+"/mass.mas03round";
    nubase_file="";
    
  } else if (name=="95exp") {
    
    std::string sha=((std::string)"bf8f6fb685100467b2f522d8fb48")+
      "62089f4843b2354be3653e5b67488294eeb3";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("mass_exp.mas95",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame95/mass_exp.mas95",sha,ext_data);

    filename=ext_data+"/mass_exp.mas95";
    nubase_file="";
    
  } else if (name=="95rmd") {
    
    std::string sha=((std::string)"f05e9bf4041f2921f82a96186452b6")+
      "d5d23d57ae4c62f476e5dea40a41e60943";
    cf.hash_type=cloud_file::sha256;
    cf.get_file_hash
      ("mass_rmd.mas95",
       ((string)"https://awsteiner.org/")+
       "public_data/nucmass/ame95/mass_rmd.mas95",sha,ext_data);

    filename=ext_data+"/mass_rmd.mas95";
    nubase_file="";
    
  } else {
    
    std::string s=((std::string)"Invalid name '")+name+
      "' in nucmass_ame::load().";
    O2SCL_ERR(s.c_str(),exc_einval);
    
  }
  
  if (verbose>0) {
    std::cout << "nucmass_ame()::load(): "
              << "name,filename,nubase_file: " << name << " "
              << filename << " " << nubase_file << std::endl;
  }
  
  load_ext(name,filename,nubase_file,exp_only,verbose);
  
  return;
}
