/*
  -------------------------------------------------------------------
  
  Copyright (C) 2009-2017, Andrew W. Steiner
  
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
  Just a short program to analyze and summarize O2scl test results.
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
// Need cstdlib for exit()
#include <cstdlib>
#include <sstream>

using namespace std;

int main(int argv, char *argc[]) {

  ifstream fin("testlist");
  if (!fin) {
    cerr << "Could not open list of files to test." << endl;
    exit(-1);
  }

  vector<string> fns;

  if (argv<2 || ((string)argc[1])!="summary") {
    
    while(!fin.eof()) {
      string tmp;
      fin >> tmp;
      if (tmp.length()>0) {
	fns.push_back(tmp);
      }
    }

  } else {

    while(!fin.eof()) {
      string tmp, dir;
      getline(fin,tmp);
      istringstream ins(tmp);
      if (ins >> dir) {
	string file;
	while (ins >> file) {
	  fns.push_back(dir+"/"+file);
	}
      }
    }

  }
  
  cout << "Summarizing test results from files: " << endl;
  for(size_t i=0;i<fns.size();i++) {
    cout << fns[i] << endl;
  }
  cout << endl;
  
  bool failed=false;

  int ntests=0;

  const size_t nch=300;
  vector<string>::iterator it;
  for(it=fns.begin();it!=fns.end();it++) {
    FILE *fp=fopen(it->c_str(),"r");
    if (fp==0) {
      cerr << "Failed to open test file: " << (*it) << endl;
      failed=true;
    } else {
      fseek(fp,-nch,SEEK_END);
      vector<string> lines;
      while (!feof(fp)) {
	char ctmp[nch];
	char *cret=fgets(ctmp,nch,fp);
	// Function fgets returns &ctmp[0] on success and 0 on failure
	if (cret==&ctmp[0]) {
	  lines.push_back(ctmp);
	}
      }
      fclose(fp);
      if (lines.size()>=2) {
	if (lines[lines.size()-1]!="All tests passed.\n") {
	  cerr << "Some tests failed in file: " << (*it) << endl;
	  string cmd=((string)"cat ")+*it;
	  int sret=system(cmd.c_str());
	  failed=true;
	} else {
	  istringstream ins(lines[lines.size()-2]);
	  int itmp;
	  ins >> itmp;
	  ntests+=itmp;
	}
      } else {
	cerr << "Some tests failed in file: " << (*it) << endl;
	string cmd=((string)"cat ")+*it;
	int sret=system(cmd.c_str());
	failed=true;
      }
    }
  }

  if (failed==true) {
    cout << "One or more O2scl tests failed." << endl;
    exit(-1);
  }

  cout << ntests << " total tests." << endl;
  cout << "All O2scl tests passed." << endl;

  return 0;
}
