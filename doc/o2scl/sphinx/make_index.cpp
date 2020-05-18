#include <iostream>
#include <fstream>
#include <vector>
#include <set>

using namespace std;

int main(void) {

  for (size_t kk=0;kk<2;kk++) {

    // Read list
    string fname_in;
    if (kk==0) fname_in="class_list";
    else fname_in="function_list";
    ifstream fin(fname_in);
    
    std::set<std::string> list;
    std::set<std::string> list_dup;
    
    while (!fin.eof()) {
      
      string s;
      std::getline(fin,s);
      
      if (s.length()>0) {
	
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
	string ns;
	if (s.find("::")!=std::string::npos) {
	  size_t loc=s.find("::");
	  ns=s.substr(0,loc);
	  s=s.substr(loc+2,s.length()-loc-2);
	}
	if (kk==0) {
	  cout << "Namespace: " << ns << " class: " << s << endl;
	} else {
	  cout << "Namespace: " << ns << " function: " << s << endl;
	}
	if (list.find(s)!=list.end()) {
	  cout << "Multiple entry " << s << endl;
	  list_dup.insert(s);
	} else {
	  list.insert(s);
	}

	// End of 'if (s.length()>0)'
      }
      
      // End of 'while (!fin.eof())'
    }

    fin.close();

    for (std::set<std::string>::itertor it=list.begin();
	 it!=list.end();it++) {

      if (list_dup.find(*it)==list_dup.end()) {
	string fname_out="class/";
	if (kk==1) fname_out="function/";
	fname_out+=s+".rst";
	ofstream fout(fname_out);
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
	cout << "Skipping " << *it << endl;
      }
    }
    
  }
  
  return 0;
}
