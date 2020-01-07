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
#include <iostream>

#include <o2scl/table3d.h>
#include <o2scl/constants.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_ext;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  string stemp;
  ifstream fin;
  fin.open("../../data/o2scl/pdg_mass_list.txt");
  
  for(size_t j=0;j<38;j++) {
    getline(fin,stemp);
  }

  while (!fin.eof()) {
    getline(fin,stemp);
    // For the IDs, if there are no numerical values, then its
    // blank so just set it equal to zero
    int id1=o2scl::stoi(stemp.substr(0,8));
    int id2=o2scl::stoi(stemp.substr(8,8));
    int id3=o2scl::stoi(stemp.substr(16,8));
    int id4=o2scl::stoi(stemp.substr(24,8));
    // if the mass is blank, then the errors are blank also
    double mass=o2scl::stod(stemp.substr(33,17));
    double mass_errp=o2scl::stod(stemp.substr(52,8));
    double mass_errm=o2scl::stod(stemp.substr(61,8));
    // if the width is blank, then the errors are blank also
    double width=o2scl::stod(stemp.substr(70,18));
    double width_errp=o2scl::stod(stemp.substr(89,8));
    double width_errm=o2scl::stod(stemp.substr(98,8));
    // Split by spaces, then split the charges by commas
    string name_charge=o2scl::stod(stemp.substr(107,21));
  }

  return 0;
}

