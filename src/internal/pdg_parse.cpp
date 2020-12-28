/*
  -------------------------------------------------------------------

  Copyright (C) 2012-2021, Andrew W. Steiner

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
#include <iostream>
#include <cctype>

#include <o2scl/table3d.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

bool has_digit(std::string s) {
  for(size_t i=0;i<s.length();i++) {
    if (std::isdigit(s[i])) return true;
  }
  return false;
}

typedef struct pdg_entry_s {
  int id;
  double mass;
  double mass_errp;
  double mass_errm;
  double width;
  double width_errp;
  double width_errm;
  string name;
  int charge;
  int isospin;
  int g_parity;
  int total_spin;
  int p_parity;
  int c_parity;
  int part_antipart_flag;
  int rank;
  char status;
  string quark_content;
} pdg_entry;

int main(int argc, char *argv[]) {

  cout.precision(10);

  vector<pdg_entry> db1, db2;
  
  string stemp;
  ifstream fin;
  fin.open("../../data/o2scl/pdg_mass_list.txt");

  // Find the first line which is not part of the header
  getline(fin,stemp);
  while (stemp[0]=='*') {
    getline(fin,stemp);
  }

  while (!fin.eof()) {

    // In the 2019 data file, one of the lines had a couple
    // extra spaces at the end.
    if (stemp.length()>=128) {
      if (false) {
	cout << stemp.substr(0,8) << "x";
	cout << stemp.substr(8,8) << "x";
	cout << stemp.substr(16,8) << "x";
	cout << stemp.substr(24,8) << "x";
	cout << stemp.substr(33,17) << "x";
	cout << stemp.substr(52,8) << "x";
	cout << stemp.substr(61,8) << "x";
	cout << stemp.substr(70,18) << "x";
	cout << stemp.substr(89,8) << "x";
	cout << stemp.substr(98,8) << "x";
	cout << stemp.substr(107,21) << endl;
      }

      pdg_entry pe;
      
      // For the IDs, if there are no numerical values, then its
      // blank so just set it equal to zero
      int id1=o2scl::stoi(stemp.substr(0,8));
      int id2;
      if (has_digit(stemp.substr(8,8))) {
	id2=o2scl::stoi(stemp.substr(8,8));
      } else {
	id2=0;
      }
      int id3;
      if (has_digit(stemp.substr(16,8))) {
	id3=o2scl::stoi(stemp.substr(16,8));
      } else {
	id3=0;
      }
      int id4;
      if (has_digit(stemp.substr(24,8))) {
	id4=o2scl::stoi(stemp.substr(24,8));
      } else {
	id4=0;
      }
      // if the mass is blank, then the errors are blank also
      double mass=-1.0, mass_errp=-1.0, mass_errm=+1.0;
      if (has_digit(stemp.substr(33,17))) {
	mass=o2scl::stod(stemp.substr(33,17));
	mass_errp=o2scl::stod(stemp.substr(52,8));
	mass_errm=o2scl::stod(stemp.substr(61,8));
      }
      // if the width is blank, then the errors are blank also
      double width=-1.0, width_errp=-1.0, width_errm=+1.0;
      if (has_digit(stemp.substr(70,18))) {
	width=o2scl::stod(stemp.substr(70,18));
	width_errp=o2scl::stod(stemp.substr(89,8));
	width_errm=o2scl::stod(stemp.substr(98,8));
      }
      // Split by spaces, then split the charges by commas
      string name_charge=stemp.substr(107,21);
      istringstream ins(name_charge);
      string name, charge;
      ins >> name >> charge;
      vector<string> charge_list;
      split_string_delim(charge,charge_list,',');
      for(size_t i=0;i<charge_list.size();i++) {
	cout << "{";
	if (i==0) {
	  cout << id1 << ",";
	} else if (i==1) {
	  cout << id2 << ",";
	} else if (i==2) {
	  cout << id3 << ",";
	} else {
	  cout << id4 << ",";
	}	  
	cout << mass << "," << mass_errp << "," << mass_errm << ",\n"
	     << width << "," << width_errp << "," << width_errm << ",\""
	     << name << "\",";
	if (charge_list[i]=="0") {
	  cout << "0}," << endl;
	} else if (charge_list[i]=="+") {
	  cout << "3}," << endl;
	} else if (charge_list[i]=="-") {
	  cout << "-3}," << endl;
	} else if (charge_list[i]=="++") {
	  cout << "6}," << endl;
	} else if (charge_list[i]=="-1/3") {
	  cout << "-1}," << endl;
	} else if (charge_list[i]=="+2/3") {
	  cout << "2}," << endl;
	} else {
	  O2SCL_ERR("Cannot interpret charge state.",o2scl::exc_efailed);
	}
      }
    }
    
    // Get the next line
    getline(fin,stemp);
  }

  fin.close();
  
  return 0;
}

