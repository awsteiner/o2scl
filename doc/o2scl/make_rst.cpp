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

/** \brief Inside string \c s, extract the element inside tag \c tag
    and store the element in \c result
*/
int xml_get_tag_element(string s, string tag, string &result) {
  // Remove enclosing <>
  if (tag[0]=='<') {
    if (tag[tag.length()-1]!='>') {
      return 1;
    }
    tag=tag.substr(1,tag.length()-2);
  }
  size_t tag_len=tag.length();
  string tag_start=((string)"<")+tag+">";
  string tag_end=((string)"</")+tag+">";

  // Exit early if there is a parsing error
  size_t istart=s.find(tag_start);
  if (istart==std::string::npos) {
    return 2;
  }
  size_t iend=s.find(tag_end);
  if (iend==std::string::npos) {
    return 3;
  }
  if (istart+tag_len+2>iend) return 4;
  
  // Extract the name
  size_t loc1=s.find("<name>");
  size_t loc2=s.find("</name>");
  size_t loct=loc1+6;
  if (loc2<loct) {
    return 5;
  }
  size_t len=iend-istart-tag_len-2;
  result=s.substr(istart+tag_len+2,len);
  return 0;
}

/*
  This code parses 'function_list' and 'class list' to create .rst
  files for each class and function (or list of overloaded functions).
*/
int main(int argc, char *argv[]) {

  if (argc<2) {
    cerr << "Requires context argument." << endl;
    exit(-1);
  }

  string context=argv[1];
  cout << "Using context " << context << endl;

  size_t kk_max=2;
  if (context==((string)"eos")) kk_max=1;
  
  // kk=0 for classes and kk=1 for functions
  for (size_t kk=0;kk<kk_max;kk++) {

    // -------------------------------------------------------
    // Parse the list file

    cout << "---------------------------------------------------" << endl;
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

	if (xml_get_tag_element(s,"name",s)!=0) {
	  cerr << "Failed to find name in " << s << endl;
	  exit(-1);
	}

	// Extract the namespace
	string ns;
	if (kk==0 && s.find("::")!=std::string::npos) {
	  size_t loc=s.find("::");
	  ns=s.substr(0,loc);
	  s=s.substr(loc+2,s.length()-loc-2);
	}
	
	// Replace &lt; with &
	while (s.find("&lt;")!=std::string::npos) {
	  s.replace(s.find("&lt;"),4,"<");
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
    cout << "Creating rst files for non-duplicate functions:" << endl;

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
	if (context==((string)"main")) {
	  if (kk==0) {
	    fout << ":ref:`O2scl <o2scl>` : :ref:`Class List`\n" << endl;
	    fout << ".. _" << s << ":\n" << endl;
	  } else {
	    fout << ":ref:`O2scl <o2scl>` : :ref:`Function List`\n" << endl;
	  }
	} else if (context==((string)"part")) {
	  if (kk==0) {
	    fout << ":ref:`O2scl_part <o2sclp>` : :ref:`Class List`\n"
		 << endl;
	    fout << ".. _" << s << ":\n" << endl;
	  } else {
	    fout << ":ref:`O2scl_part <o2sclp>` : "
		 << ":ref:`Function List`\n" << endl;
	  }
	} else if (context==((string)"eos")) {
	  if (kk==0) {
	    fout << ":ref:`O2scl_eos <o2scle>` : :ref:`Class List`\n"
		 << endl;
	    fout << ".. _" << s << ":\n" << endl;
	  } else {
	    fout << ":ref:`O2scl_eos <o2scle>` : "
		 << ":ref:`Function List`\n" << endl;
	  }
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
    cout << "Creating rst files for duplicate functions:" << endl;
    
    // -------------------------------------------------------
    // For functions, create the rst files for duplicate entries
    
    if (kk==1) {

      // Open the namespace xml file and read it into memory
      vector<string> ns_file;
      if (context==((string)"main")) {
	ifstream fin("../xml/namespaceo2scl.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
	fin.open("../xml/namespaceo2scl__acol.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
	fin.open("../xml/namespaceo2scl__cblas.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
	fin.open("../xml/namespaceo2scl__hdf.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
	fin.open("../xml/namespaceo2scl__linalg.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
      } else if (context==((string)"part")) {
	ifstream fin("../../xml/namespaceo2scl.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
	fin.open("../../xml/namespaceo2scl__hdf.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
      } else if (context==((string)"eos")) {
	ifstream fin("../../xml/namespaceo2scl.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
	fin.open("../../xml/namespaceo2scl__hdf.xml");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  ns_file.push_back(s2);
	}
	fin.close();
      }

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
	if (context==((string)"main")) {
	  fout << ":ref:`O2scl <o2scl>` : :ref:`Function List`\n" << endl;
	} else if (context==((string)"part")) {
	  fout << ":ref:`O2scl_part <o2sclp>` : "
	       << ":ref:`Function List`\n" << endl;
	} else if (context==((string)"eos")) {
	  fout << ":ref:`O2scl_eos <o2scle>` : "
	       << ":ref:`Function List`\n" << endl;
	}
	
	string head="Functions "+s;
	fout << head << endl;
	for(size_t i=0;i<head.length();i++) {
	  fout << "=";
	}
	fout << endl;
	fout << endl;

	// Iterate through the namespace file
	for(size_t i=0;i<ns_file.size();i++) {
	  
	  string s2=ns_file[i];
	  
	  // Proceed through the file until we find the right
	  // definition
	  if (s2.find("<definition>")!=std::string::npos &&
	      s2.find(s)!=std::string::npos) {

	    size_t arg_count=0;

	    string argsstring;
	    // Read the argstring
	    i++;
	    argsstring=ns_file[i];
	    if (xml_get_tag_element(argsstring,"argsstring",argsstring)!=0) {
	      cerr << "Failed to find argsstring in " << argsstring << endl;
	      exit(-1);
	    }
	    // Replace &amp; with &
	    while (argsstring.find("&amp;")!=std::string::npos) {
	      argsstring.replace(argsstring.find("&amp;"),5,"&");
	    }
	    // Replace "&lt; " with <
	    while (argsstring.find("&lt; ")!=std::string::npos) {
	      argsstring.replace(argsstring.find("&lt; "),5,"<");
	    }
	    // Replace "&lt;" with <
	    while (argsstring.find("&lt;")!=std::string::npos) {
	      argsstring.replace(argsstring.find("&lt;"),4,"<");
	    }
	    // Replace " &gt;" with >
	    while (argsstring.find(" &gt;")!=std::string::npos) {
	      argsstring.replace(argsstring.find(" &gt;"),5,">");
	    }
	    // Replace "&gt;" with >
	    while (argsstring.find("&gt;")!=std::string::npos) {
	      argsstring.replace(argsstring.find("&gt;"),4,">");
	    }
	    // Replace &quot; with "
	    while (argsstring.find("&quot;")!=std::string::npos) {
	      argsstring.replace(argsstring.find("&quot;"),6,"\"");
	    }
	    // Replace all = with ' = '
	    vector<size_t> equal_list;
	    for(size_t j=0;j<argsstring.length();j++) {
	      if (argsstring[j]=='=') {
		argsstring.replace(j,1," = ");
		j++;
	      }
	    }

	    // Read the name line
	    string name;
	    i++;
	    name=ns_file[i];
	    
	    if (xml_get_tag_element(name,"name",name)!=0) {
	      cerr << "Failed to find name in " << name << endl;
	      exit(-1);
	    }

	    // Only proceed if the name and s match
	    if (name==s) {
	      
	      // Output the sphinx directive
	      fout << ".. doxygenfunction:: " << s << argsstring
		   << "\n" << endl;
	      cout << s << argsstring << endl;
	      
	      // Read the line <param> (if present)
	      i++;
	      s2=ns_file[i];
	      
	      //cout << s << "(";

	      // Read all of the parameters
	      bool params_done=false;
	      while (params_done==false) {

		string type, declname, end;
		i++;
		type=ns_file[i];
		if (xml_get_tag_element(type,"type",type)!=0) {
		  cerr << "Failed to find type in " << type << endl;
		  exit(-1);
		}
		
		//cout << "Type: " << type << endl;
		if (type.find(" &amp;")!=std::string::npos) {
		  type.replace(type.find(" &amp;"),6,"&");
		}
		
		i++;
		declname=ns_file[i];
		if (xml_get_tag_element(declname,"declname",declname)!=0) {
		  cerr << "Failed to find declname in " << declname << endl;
		  exit(-1);
		}
		//cout << "Declname: " << declname << endl;
		
		// Read <defval> (if present) and </param> lines
		i++;
		end=ns_file[i];
		//cout << "End: " << end << endl;
		if (end.find("<defval>")!=std::string::npos) {
		  i++;
		  end=ns_file[i];
		}
		
		// Read the next line
		i++;
		s2=ns_file[i];
		//cout << "s2: " << s2 << endl;

		/*
		// Update file and screen output
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
		*/
	      
		arg_count++;

		if (s2.find("<param>")==std::string::npos) {
		  params_done=true;
		  //cout << "X " << s2 << endl;
		}
		
		//char ch;
		//cin >> ch;
	      }

	      // Output right parenthesis
	      //cout << ")" << endl;
	      //fout << ")\n" << endl;
	      
	    } else {
	      //cout << "No match " << name << " " << s << endl;
	    }
	  }
	    
	}

	// Close this rst file
	fout.close();
      }
      
    }
    
  }
  
  return 0;
}
