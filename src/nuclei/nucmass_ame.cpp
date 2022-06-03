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
#include <o2scl/hdf_nucmass_io.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

nucmass_ame::nucmass_ame() {
  n=0;
  reference="";
  mass=0;
  last=0;
}

nucmass_ame::~nucmass_ame() {
  if (n>0) {
    delete[] mass;
  }
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
  bool found=false;
  for(size_t i=0;i<n && found==false;i++) {
    if (mass[i].Z==l_Z && mass[i].N==l_N) {
      ret=mass[i];
      found=true;
    }
  }
  if (found==false) {
    std::string err=((std::string)"Nucleus with (Z,N)=(")+
      o2scl::itos(l_Z)+","+o2scl::itos(l_N)+") not found in "+
      "nucmass_ame::get_ZN().";
    O2SCL_ERR(err.c_str(),exc_einval);
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
    ret.N=0;
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

nucmass_ame2::nucmass_ame2() {
  reference="";
  last=0;
}

nucmass_ame2::~nucmass_ame2() {
}

double nucmass_ame2::mass_excess(int Z, int N) {
  entry ret;
  ret=get_ZN(Z,N);
  if (ret.Z==0 && ret.N==0) return 0.0;
  return ret.mass/1.0e3;
}

bool nucmass_ame2::is_included(int l_Z, int l_N) {

  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame2::is_included().",
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
bool nucmass_ame2_exp::is_included(int l_Z, int l_N) {
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame2_exp::is_included().",
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

nucmass_ame2::entry nucmass_ame2::get_ZN(int l_Z, int l_N) {
  nucmass_ame2::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame2::get_ZN().",
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

nucmass_ame2::entry nucmass_ame2::get_ZA(int l_Z, int l_A) {
  nucmass_ame2::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame2::get_ZA().",
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

nucmass_ame2::entry nucmass_ame2::get_elA(string l_el, int l_A) {
  nucmass_ame2::entry ret;
  ret.Z=0;
  ret.A=0;
  ret.N=0;
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame2::get_elA().",
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

nucmass_ame2::entry nucmass_ame2::get(string nucleus) {
  if (n==0) {
    O2SCL_ERR("No masses loaded in nucmass_ame2::get().",
	      exc_einval);
    nucmass_ame2::entry ret;
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

void nucmass_ame2::load_ext(std::string model, std::string filename,
                            std::string nubase_file, bool exp_only) {

  nucmass_ame2::entry ae;

  ifstream fin(filename.c_str());
  std::string line;

  int count=0;

  int verbose=0;

  mass.clear();

  if (model=="20") {
    for(size_t i=0;i<36;i++) getline(fin,line);
  } else if (model=="20round") {
    for(size_t i=0;i<34;i++) getline(fin,line);
  }
  
  while (getline(fin,line)) {

    vector<string> entries;
    
    if (model=="20") {
      
      parse_fortran_format(line,((string)"a1,i3,i5,i5,i5,1x,a3,a4,")+
                           "1x,f14.6,f12.6,f15.5,f11.5,1x,a2,"+
                           "f13.5,f11.5,1x,i3,1x,f13.6,f12.6",entries);
      
      //} else if (model=="20round") {
      
    }

    /*
      if (verbose>2) {
      cout << "Entries:" << endl;
      for(size_t k=0;k<entries.size();k++) {
      cout << k << " " << entries[k].size() << " '" << entries[k] << "'." << endl;
      }
      }
    */
      
    // Read N and Z first to compute NMZ and A
    ae.N=o2scl::stoi(entries[2]);
    ae.Z=o2scl::stoi(entries[3]);
    if (model=="20round" || model=="16round" || model=="03round") {
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
    
    if (model=="95rmd" || model=="95exp") {
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
    
    mass.push_back(ae);
      
    count++;
    
  }

  cout << "count: " << count << endl;
  n=count;

  last=mass.size()/2;

  if (model=="20" || model=="20round") {
    
    reference=((string)"W. J. Huang, M. Wang ")+
      "F. G. Kondev, G. Audi, S. Naimi, X. Xu, "
      "Chin. Phys. C, 45 (2021) 030002; "+
      "M. Wang, W. J. Huang, F. G. Kondev, "+
      "G. Audi, and S. Naimi, "
      "Chin. Phys. C, 41 (2021) 030003.";
  }

  fin.close();

  if (model=="20") {
    
    fin.open(nubase_file.c_str());
    
    count=0;

    if (model=="20") {
      for(size_t i=0;i<25;i++) getline(fin,line);
    }
    
    while (getline(fin,line)) {

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

      entry_nubase_20 enu20;
      if (entries[1].length()>3 && entries[1][3]!=' ') {
        enu20.Znote=entries[1][3];
      } else {
        enu20.Znote='\0';
      }
      if (verbose>1) cout << "Znote: " << ((int)enu20.Znote) << endl;
      string_to_char_array(entries[2],enu20.A_el,6);
      if (verbose>1) cout << "A_el: '" << enu20.A_el << "'" << endl;
      if (entries[3].length()>0 && entries[4][0]!=' ') {
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
        enu20.exc_energy_acc=nucmass_ame2::does_not_exist;
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
      if (verbose>1) cout << "isomer_unc: '" << enu20.isomer_unc << "'" << endl;
      if (entries.size()>10 && entries[10].length()>0 && entries[10][0]!=' ') {
        enu20.isomer_inv=entries[10][0];
      } else {
        enu20.isomer_inv='\0';
      }
      if (verbose>1) cout << "isomer_inv: '" << enu20.isomer_inv << "'" << endl;
      if (entries.size()>11 && count_words(entries[11])>0) {
        remove_whitespace(entries[11]);
        if (entries[11]=="stbl") {
          enu20.hlife=0.0;
          enu20.hlife_acc=nucmass_ame2::unstable;
        } else if (entries[11]=="p-unst") {
          enu20.hlife=0.0;
          enu20.hlife_acc=nucmass_ame2::part_unstable;
        } else if (entries[11][0]=='>') {
          entries[11]=entries[11].substr
            (1,entries[11].length()-1);
          enu20.hlife=o2scl::stod(entries[11]);
          enu20.hlife_acc=nucmass_ame2::lower_limit;
        } else if (entries[11][0]=='<') {
          entries[11]=entries[11].substr
            (1,entries[11].length()-1);
          enu20.hlife=o2scl::stod(entries[11]);
          enu20.hlife_acc=nucmass_ame2::upper_limit;
        } else if (entries[11][0]=='~') {
          entries[11]=entries[11].substr
            (1,entries[11].length()-1);
          enu20.hlife=o2scl::stod(entries[11]);
          enu20.hlife_acc=nucmass_ame2::approximate;
        } else {
          enu20.hlife=o2scl::stod(entries[11]);
          enu20.hlife_acc=0;
        }
      } else {
        enu20.hlife=0.0;
        enu20.hlife_acc=nucmass_ame2::blank;
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
        string_to_char_array(entries[17],enu20.decay_intensity,91);
      } else {
        enu20.decay_intensity[0]='\0';
      }
      if (verbose>1) {
        cout << "decay_intensity: "
             << enu20.decay_intensity << endl;
      }

      mass[ix].props.push_back(enu20);
    }

    fin.close();

    //for(size_t i=0;i<n;i++) {
    //cout << i << " " << mass[i].props.size() << endl;
    //}
    
  }
  
  return;
}
