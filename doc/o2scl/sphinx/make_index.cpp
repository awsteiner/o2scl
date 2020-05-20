/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020, Andrew W. Steiner
  
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
#include <vector>
#include <set>
#include <map>

using namespace std;

int main(void) {

  // kk=0 for classes and kk=1 for functions
  for (size_t kk=0;kk<2;kk++) {

    // -------------------------------------------------------
    // Parse the list file

    if (kk==0) {
      cout << "Parsing class_list:" << endl;
    } else {
      cout << "Parsing function_list:" << endl;
    }
    
    // Open the file of functions or classes
    string fname_in;
    if (kk==0) fname_in="class_list";
    else fname_in="function_list";
    ifstream fin(fname_in);

    // List of items, name as "first" and namespace as "second"
    std::map<std::string,std::string> list;
    // list of duplicate items (names only)
    std::set<std::string> list_dup;

    while (!fin.eof()) {
      
      string s;
      std::getline(fin,s);

      // There is sometimes a line with zero length at the end,
      // so we skip it
      if (s.length()>0) {

	// Exit early if there is a parsing error
	if (s.find("<name>")==std::string::npos) {
	  cout << "No <name>" << endl;
	  cout << s << endl;
	  exit(-1);
	}
	if (s.find("</name>")==std::string::npos) {
	  cout << "No </name>" << endl;
	  cout << s << endl;
	  exit(-1);
	}

	// Extract the name
	size_t loc1=s.find("<name>");
	size_t loc2=s.find("</name>");
	size_t loct=loc1+6;
	if (loc2<loct) {
	  cout << "Missing name" << endl;
	  cout << s << endl;
	  exit(-1);
	}
	size_t len=loc2-loc1-6;
	s=s.substr(loc1+6,len);

	// Extract the namespace
	string ns;
	if (kk==0 && s.find("::")!=std::string::npos) {
	  size_t loc=s.find("::");
	  ns=s.substr(0,loc);
	  s=s.substr(loc+2,s.length()-loc-2);
	}
	
	if (kk==0) {
	  cout << "Namespace: " << ns << " class: " << s << endl;
	} else {
	  cout << "Namespace: " << ns << " function: " << s << endl;
	}

	// Handle multiple entry or add to the main list
	if (list.find(s)!=list.end()) {
	  cout << "Multiple entry " << s << endl;
	  list_dup.insert(s);
	  if (kk==0) {
	    cerr << "Duplicate classes not allowed." << endl;
	    exit(-1);
	  }
	} else if (ns!=((string)"boost")) {
	  list.insert(std::make_pair(s,ns));
	}

	// End of 'if (s.length()>0)'
      }
      
      // End of 'while (!fin.eof())'
    }

    fin.close();

    cout << endl;

    // -------------------------------------------------------
    // Now create the rst files for the non-duplicate entries
    
    cout << "---------------------------------------------------" << endl;
    cout << "Creating rst files:" << endl;

    // Proceed through the list
    for (std::map<std::string,std::string>::iterator it=list.begin();
	 it!=list.end();it++) {

      string s=it->first, ns=it->second;

      // Make sure it's not a duplicate
      if (list_dup.find(s)==list_dup.end()) {

	// Open the rst file
	string fname_out="class/";
	if (kk==1) fname_out="function/";
	fname_out+=s+".rst";
	ofstream fout(fname_out);

	// Header
	if (kk==0) {
	  fout << ":ref:`O2scl <o2scl>` : :ref:`Class List`\n" << endl;
	  fout << ".. _" << s << ":\n" << endl;
	} else {
	  fout << ":ref:`O2scl <o2scl>` : :ref:`Function List`\n" << endl;
	}
	
	string head="Class "+s;
	if (kk==1) head="Function "+s;
	fout << head << endl;
	for(size_t i=0;i<head.length();i++) {
	  fout << "=";
	}
	fout << endl;
	fout << endl;

	// Output the class or function directive
	if (kk==0) {
	  fout << ".. doxygenclass:: " << ns << "::" << s << endl;
	} else {
	  fout << ".. doxygenfunction:: " << ns << "::" << s << endl;
	}

	// Close the file
	fout.close();

      } else {
	// Skip items in the list of duplicates
	cout << "Skipping " << it->first << endl;
      }
    }

    cout << endl;
    cout << "---------------------------------------------------" << endl;

    // -------------------------------------------------------
    // For functions, create the rst files for duplicate entries
    
    if (kk==1) {

      // Iterate over the list
      for (std::set<std::string>::iterator it=list_dup.begin();
	   it!=list_dup.end();it++) {

	// The function name
	string s=*it;

	// Open the rst file
	string fname_out="function/";
	fname_out+=s+".rst";
	ofstream fout(fname_out);

	// Header
	fout << ":ref:`O2scl <o2scl>` : :ref:`Function List`\n" << endl;
	
	string head="Functions "+s;
	fout << head << endl;
	for(size_t i=0;i<head.length();i++) {
	  fout << "=";
	}
	fout << endl;
	fout << endl;

	// Now open the namespace xml file and read it to find the various
	// argument lists
	ifstream fin("../xml/namespaceo2scl.xml");
	
	while (!fin.eof()) {
	  
	  string s2;
	  std::getline(fin,s2);

	  // Proceed through the file until we find the right
	  // definition
	  if (s2.find("<definition>")!=std::string::npos &&
	      s2.find(s)!=std::string::npos) {

	    // Output the name part of the sphinx directive
	    fout << ".. doxygenfunction:: " << s << "(";

	    size_t arg_count=0;
	    // Skip the argstring
	    std::getline(fin,s2);
	    // Skip the name
	    std::getline(fin,s2);
	    // Read the param
	    std::getline(fin,s2);
	    
	    cout << s << "(";
	    
	    while (s2.find("<param>")!=std::string::npos) {
	      string type, declname, end;
	      std::getline(fin,type);
	      if (type.length()>23) {
		type=type.substr(16,type.length()-23);
	      }
	      if (type.find(" &amp;")!=std::string::npos) {
		type.replace(type.find(" &amp;"),6,"&");
	      }
	      std::getline(fin,declname);
	      if (declname.length()>31) {
		declname=declname.substr(20,declname.length()-31);
	      }
	      std::getline(fin,end);
	      std::getline(fin,s2);
	      
	      if (arg_count!=0) {
		cout << ", " << type;
	      } else {
		cout << type;
	      }

	      if (arg_count!=0) {
		fout << ", " << type;
	      } else {
		fout << type;
	      }
	      arg_count++;
	      //char ch;
	      //cin >> ch;
	    }
	    
	    cout << ")" << endl;
	    fout << ")\n" << endl;
	  }
	    
	}
	fin.close();
	fout.close();
      }
    }
    
  }
  
  return 0;
}
