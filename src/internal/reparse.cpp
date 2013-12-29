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
  Reparse the EOS data files from the old to the new format
*/

#include <iostream>
#include <string>

// for exit()
#include <cstdlib>

#include <o2scl/hdf_file.h>
#include <o2scl/misc.h>

using namespace std;
using namespace o2scl;

int main(int argc, char *argv[]) {

  if (argc<3) {
    cout << "Usage: reparse <input filename> <output filename>" << endl;
    exit(-1);
  }

  cout << argv[1] << " " << argv[2] << endl;
  string fname=argv[1];
  string outname=argv[2];
  ifstream fin(fname.c_str());
  hdf_file hdf;
  hdf.open(outname.c_str());

  string type, name;
  while (fin >> type) {
    if (type[0]=='#') {
      getline(fin,type);
    } else if (type=="double") {
      double d;
      fin >> name;
      fin >> d;
      cout << "double " << name << " " << d << endl;
      hdf.setd(name,d);
    } else if (type=="bool") {
      bool b;
      fin >> name;
      string tmp;
      fin >> tmp;
      cout << "bool " << name << " " << tmp << endl;
      if (tmp=="true" || tmp=="1") {
	hdf.seti(name,1);
      } else {
	hdf.seti(name,0);
      }
    } else if (type=="int") {
      int i;
      fin >> name;
      fin >> i;
      cout << "int " << name << " " << i << endl;
      hdf.seti(name,i);
    } else if (type=="string") {
      string s;
      fin >> name;
      getline(fin,s);
      if (count_words(s)==0) {
	getline(fin,s);
      }
      cout << "string " << name << " " << s << endl;
      hdf.sets(name,s);
    } else {
      cout << "Unknown type: " << type << endl;
      exit(-1);
    }
  }

  fin.close();
  hdf.close();

  return 0;
}
