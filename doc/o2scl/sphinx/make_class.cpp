#include <iostream>
#include <fstream>
#include <vector>
#include <set>

using namespace std;

int main(void) {

  string s;
  ifstream fin("class_list");
  std::set<std::string> list;
  while (fin >> s) {
    string ns;
    if (s.find("::")!=std::string::npos) {
      size_t loc=s.find("::");
      ns=s.substr(0,loc);
      s=s.substr(loc+2,s.length()-loc-2);
    }
    cout << "Namespace: " << ns << " class: " << s << endl;
    if (list.find(s)!=list.end()) {
      cout << "Multiple entry " << s << endl;
      return 1;
    }
    list.insert(s);
    string fname="class/";
    fname+=s+".rst";
    ofstream fout(fname);
    fout << s << endl;
    for(size_t i=0;i<s.length();i++) {
      fout << "=";
    }
    fout << endl;
    fout << endl;
    fout << ".. doxygenclass:: " << ns << "::" << s << endl;
    fout.close();
  }
  
  return 0;
}
