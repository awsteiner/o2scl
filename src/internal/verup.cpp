/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
#include <sstream>
#include <string>

// for exit()
#include <cstdlib>

using namespace std;

int stoi(string s) {
  int ret;
  istringstream ins(s);
  if (ins >> ret) {
    return ret;
  }
  return 0;
}

int main(int argv, char *argc[]) {

  string line, word;
  
  if (argv<4) {
    cout << "Usage: " << endl;
    cout << " verup 1 version Doxyfile outfile" << endl;
    cout << " verup 2 version header.tex outfile" << endl;
    exit(-1);
  }

  ifstream fin(argc[3]);
  ofstream fout(argc[4]);
  if (stoi(argc[1])==1) {

    // Modify Doxyfile
    
    while (getline(fin,line)) {
      istringstream *ins=new istringstream(line);
      (*ins) >> word;
      if (line.length()>0 && word=="PROJECT_NUMBER") {
	fout << "PROJECT_NUMBER = " << argc[2] << endl;
      } else {
	fout << line << endl;
      }
      delete ins;
    }

  } else {

    // Modify header.tex

    while (getline(fin,line)) {
      istringstream *ins=new istringstream(line);
      (*ins) >> word;
      if (line.length()>0 && word=="Version") {
	fout << "Version\n" << argc[2] << endl;
	delete ins;
	getline(fin,line);
      } else {
	fout << line << endl;
	delete ins;
      }
    }

  }
  fin.close();
  fout.close();
  return 0;
}
