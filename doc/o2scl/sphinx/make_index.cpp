#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>

using namespace std;

int main(void) {

  // kk=0 for classes and kk=1 for functions
  for (size_t kk=0;kk<2;kk++) {

    if (kk==0) {
      cout << "Parsing class_list:" << endl;
    } else {
      cout << "Parsing function_list:" << endl;
    }
    
    // Read list
    string fname_in;
    if (kk==0) fname_in="class_list";
    else fname_in="function_list";
    ifstream fin(fname_in);
    
    std::map<std::string,std::string> list;
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
	} else if (ns!=((string)"boost")) {
	  list.insert(std::make_pair(s,ns));
	}

	// End of 'if (s.length()>0)'
      }
      
      // End of 'while (!fin.eof())'
    }

    fin.close();

    cout << endl;
    cout << "Creating files:" << endl;

    for (std::map<std::string,std::string>::iterator it=list.begin();
	 it!=list.end();it++) {

      string s=it->first, ns=it->second;
      
      if (list_dup.find(s)==list_dup.end()) {
	string fname_out="class/";
	if (kk==1) fname_out="function/";
	fname_out+=s+".rst";
	ofstream fout(fname_out);

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
	if (kk==0) {
	  fout << ".. doxygenclass:: " << ns << "::" << s << endl;
	} else {
	  fout << ".. doxygenfunction:: " << ns << "::" << s << endl;
	}
	fout.close();

      } else {
	cout << "Skipping " << it->first << endl;
      }
    }

    cout << endl;
    cout << "---------------------------------------------------" << endl;

    // For functions, we need to handle duplicates separately
    if (kk==1) {

      // For items with duplicates
      for (std::set<std::string>::iterator it=list_dup.begin();
	   it!=list_dup.end();it++) {

	string s=*it;
	string cmd="grep -1 \"";
	cmd+=s+"</definition>\" ../xml/namespaceo2scl.xml | "+
	  "grep -i argsstring > /tmp/mi.txt";
	int ret=system(cmd.c_str());

	string fname_out="function/";
	fname_out+=s+".rst";
	ofstream fout(fname_out);

	fout << ":ref:`O2scl <o2scl>` : :ref:`Function List`\n" << endl;
	
	string head="Function "+s;
	fout << head << endl;
	for(size_t i=0;i<head.length();i++) {
	  fout << "=";
	}
	fout << endl;
	fout << endl;

	ifstream fin("/tmp/mi.txt");
	while (!fin.eof()) {
	  string s2;
	  std::getline(fin,s2);
	  if (s2.length()>=27) {
	    s2=s2.substr(20,s2.length()-33);
	    cout << s << " " << s2 << endl;
	    fout << ".. doxygenfunction:: " << s
		 << s2 << "\n" << endl;
	  }
	}
	fin.close();
	fout.close();
      }
    }
    
  }
  
  return 0;
}
