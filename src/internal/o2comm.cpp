/*
  -------------------------------------------------------------------

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

  -------------------------------------------------------------------
*/
#include <iostream>
#include <fstream>
#include <string>

#include <o2scl/string_conv.h>
#include <o2scl/misc.h>
#include <o2scl/hdf_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  hdf_file hf;
  hf.open(argv[1]);

  if (argc>2) {
    
    string line, final;
    ifstream fin;
    fin.open(argv[2]);
    while (getline(fin,line)) {
      if (final.length()!=0) final+='\n';
      final+=line;
    }
    fin.close();
  
    hf.sets("comment",final);
    
  }

  string olds;
  hf.gets("comment",olds);
  cout << olds << endl;

  hf.close();

  return 0;
}
