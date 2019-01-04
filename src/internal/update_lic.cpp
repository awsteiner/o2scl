/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
  AWS: This is just some quick code to verify that the license information
  at the top of the source files is correct.
 */
#include <iostream>
#include <fstream>

using namespace std;

int llines=22;

std::string base[22]={
  "/*",
  "  -------------------------------------------------------------------",
  "  ",
  "  Copyright (C) 2006-2018, Andrew W. Steiner",
  "  ",
  "  This file is part of O2scl.",
  "  ",
  "  O2scl is free software; you can redistribute it and/or modify",
  "  it under the terms of the GNU General Public License as published by",
  "  the Free Software Foundation; either version 3 of the License, or",
  "  (at your option) any later version.",
  "  ",
  "  O2scl is distributed in the hope that it will be useful,",
  "  but WITHOUT ANY WARRANTY; without even the implied warranty of",
  "  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the",
  "  GNU General Public License for more details.",
  "  ",
  "  You should have received a copy of the GNU General Public License",
  "  along with O2scl. If not, see <http://www.gnu.org/licenses/>.",
  "",
  "  -------------------------------------------------------------------",
  "*/"};

int main(int argc, char *argv[]) {
  std::string fn=argv[1];
  ifstream fin(fn.c_str());
  std::string flist[llines];
  bool match=true;
  cout << "Checking: " << fn << endl;
  for(size_t i=0;i<llines;i++) {
    getline(fin,flist[i]);
    if (i!=3 && i!=2 && i!=4 && i!=6 && i!=11 && i!=16 && i!=19 && 
	flist[i]!=base[i]) {
      match=false;
      cout << "Failed on line " << i << ": " << endl;
      cout << "Base: " << base[i] << endl;
      cout << "File: " << flist[i] << endl;
    }
  }
  fin.close();
  if (match==true) {
    cout << "Matched." << endl;
  }
  /*
  ofstream fout("temp");
  for(size_t i=0;i<llines-1;i++) {
    fout << llist[i] << endl;
  }
  std::string stmp;
  while (getline(fin,stmp)) {
    fout << stmp << endl;
  }
  fin.close();
  fout.close();
  
  std::string cmd;
  cmd="mv temp ";
  cmd+=fn;
  system(cmd.c_str());
  */

  return 0;
}
