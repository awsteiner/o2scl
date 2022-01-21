/*
  -------------------------------------------------------------------
  
  Copyright (C) 2020-2022, Andrew W. Steiner
  
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

#ifdef O2SCL_PUGIXML
#include "pugixml.hpp"
#endif

using namespace std;

/** \brief Convert <,>, ,: to _
 */
std::string underscoreify(std::string s) {
  std::string s2=s;
  for(size_t i=0;i<s2.length();i++) {
    if (s2[i]==' ' || s2[i]=='<' || s2[i]==':' || s2[i]=='>') {
      s2[i]='_';
    }
  }
  while (s2.find("__")!=std::string::npos) {
    s2.replace(s2.find("__"),2,"_");
  }
  return s2;
}

/** \brief Inside string \c s, extract the element inside tag \c tag
    and store the element in \c result
*/
int xml_get_tag_element(string s, string tag, string &result) {
  // Remove enclosing <> in tag
  if (tag[0]=='<') {
    if (tag[tag.length()-1]!='>') {
      return 1;
    }
    tag=tag.substr(1,tag.length()-2);
  }
  size_t tag_len=tag.length();
  string tag_start=((string)"<")+tag+">";
  string tag_end=((string)"</")+tag+">";

  // Exit early if the tag is not found
  size_t istart=s.find(tag_start);
  if (istart==std::string::npos) {
    return 2;
  }
  size_t iend=s.find(tag_end);
  if (iend==std::string::npos) {
    return 3;
  }
  if (istart+tag_len+2>iend) return 4;
  
  // Make sure the <name> tag is also present in the string somewhere,
  // return early if it is not
  size_t loc1=s.find("<name>");
  size_t loc2=s.find("</name>");
  size_t loct=loc1+6;
  if (loc2<loct) {
    return 5;
  }
  
  // Finally, return the result, adjusting by 2 for the < and >
  size_t len=iend-istart-tag_len-2;
  result=s.substr(istart+tag_len+2,len);
  return 0;
}

int main(int argc, char *argv[]) {

#ifdef O2SCL_PUGIXML
  if (argc<3) {
    cerr << "Requires input file and context argument, "
         << "either \"eos\" or \"main\"."
         << endl;
    exit(-1);
  }

  std::string in_file=argv[1];
  cout << "Using input file " << in_file << endl;
  string context=argv[2];
  cout << "Using context " << context << endl;

  // -----------------------------------------------------------------
  // Step 1: Parse the doxygen XML to find all the classes and
  // functions
  
  pugi::xml_document doc;

  vector<std::string> class_list;
  vector<std::string> function_list[2];

  pugi::xml_parse_result result=doc.load_file(in_file.c_str());
  if (!result) {
    cout << result.description() << endl;
    cout << "Failed to read file." << endl;
    exit(-1);
  }

  int xml_verbose=0;
  
  pugi::xml_node dindex=doc.first_child();
  
  if (xml_verbose>0) {
    std::cout << "1: " << dindex.name() << endl;
  }
  int i=0;
  
  for (pugi::xml_node_iterator it=dindex.begin();it!=dindex.end();++it) {
    if (xml_verbose>0) {
      cout << i << " " << it->name() << " " << it->child_value("name") << endl;
      cout << "  kind: " << it->attribute("kind").value() << std::endl;
    }
    if (((string)it->attribute("kind").value())==((string)"class")) {
      class_list.push_back(it->child_value("name"));
    }
    if (((string)it->attribute("kind").value())==((string)"namespace")) {
      int j=0;
      for (pugi::xml_node_iterator it2=it->begin();
           it2!=it->end();++it2) {
        if (xml_verbose>1) {
          cout << j << " " << it2->name() << " "
               << it2->child_value("name") << endl;
          cout << "  kind: " << it2->attribute("kind").value() << std::endl;
        }
        if (((string)it2->attribute("kind").value())==((string)"function")) {
          function_list[0].push_back(it->child_value("name"));
          function_list[1].push_back(it2->child_value("name"));
        }
        j++;
      }
    }
    i++;
  }

  if (xml_verbose>1) {
    for(size_t i=0;i<class_list.size();i++) {
      cout << i << " " << class_list[i] << endl;
    }
    for(size_t i=0;i<function_list[0].size();i++) {
      cout << i << " " << function_list[0][i] << " "
           << function_list[1][i] << endl;
    }
  }
  exit(-1);
  
  size_t kk_max=2;
  if (context==((string)"eos")) kk_max=1;
  //if (context==((string)"main")) kk_max=3;
  if (context==((string)"main")) kk_max=2;
  
  // kk=0 for classes and kk=1 for functions
  for (size_t kk=0;kk<kk_max;kk++) {

    // -------------------------------------------------------
    // Parse the list file

    cout << "---------------------------------------------------" << endl;
    if (kk==0) {
      cout << "Parsing class_list:" << endl;
    } else if (kk==1) {
      cout << "Parsing function_list:" << endl;
    } else {
      cout << "Parsing file_list:" << endl;
    }
    cout << endl;
    
    // Open the file of functions or classes
    string fname_in;
    if (kk==0) fname_in="class_list";
    else if (kk==1) fname_in="function_list";
    else fname_in="file_list";
    ifstream fin(fname_in);

    // List of items, name as "first" and namespace as "second"
    std::map<std::string,std::string> list;
    // list of duplicate items (names only)
    std::set<std::string> list_dup;

    while (!fin.eof()) {

      // Read each line of the file
      string s;
      std::getline(fin,s);

      // There is sometimes a line with zero length at the end,
      // so we skip it
      if (s.length()>0) {

	string ns;
        
	// Extract namespace for a function
	if (kk==1) {
	  size_t loc=s.find("namespace")+14;
	  size_t cnt=0;
	  for(size_t j=loc;j<s.length();j++) {
	    if (s[j]=='_') cnt++;
	    if (cnt==3) {
	      size_t loc2=j;
	      string stemp=s.substr(loc,loc2-loc);
	      if (stemp==((string)"__acol")) {
		ns="o2scl_acol";
	      } else if (stemp==((string)"__cblas")) {
		ns="o2scl_cblas";
	      } else if (stemp==((string)"__hdf")) {
		ns="o2scl_hdf";
	      } else if (stemp==((string)"__linalg")) {
                ns="o2scl_linalg";
	      } else {
		ns="o2scl";
	      }
	      j=s.length();
	    }
	  }
	}
	
	if (xml_get_tag_element(s,"name",s)!=0) {
	  cerr << "Failed to find name in " << s << endl;
	  exit(-1);
	}

        // Fix parsing for lu_decomp_array_2d()
        if (s==((string)"LU_decomp_array_2d")) {
          ns="o2scl_linalg_bracket";
        }
        
	// Extract the namespace for a class
	if (kk==0 && s.find("::")!=std::string::npos) {
	  size_t loc=s.find("::");
	  ns=s.substr(0,loc);
          if (ns=="o2scl" || ns=="o2scl_auto_format" ||
              ns=="o2scl_linalg" || ns=="o2scl_hdf" ||
              ns=="o2scl_const" || ns=="o2scl_acol" ||
              ns=="o2scl_cgs" || ns=="o2scl_cgsm" ||
              ns=="o2scl_mks") {
            s=s.substr(loc+2,s.length()-loc-2);
          } else {
            ns="";
          }
	}
	
	// Replace &lt; with <
	while (s.find("&lt;")!=std::string::npos) {
	  s.replace(s.find("&lt;"),4,"<");
	}
	// Replace &gt; with >
	while (s.find("&gt;")!=std::string::npos) {
	  s.replace(s.find("&gt;"),4,">");
	}
      
        if (s==((string)"operator<<")) {
          if (ns=="o2scl") {
            ns="o2scl_auto_format";
          } else {
            ns="o2scl";
          }
        }
        
	if (kk==0) {
	  cout << "Namespace: " << ns << " class: " << s << endl;
	} else if (kk==1) {
	  cout << "Namespace: " << ns << " function: " << s << endl;
	} else { 
	  cout << "File: " << s << endl;
	}

        //cout << "Here4: " << ns << " " << s << endl;

        // Polytrope solve is already documented as part of the
        // nstar_rot class, so we don't need to duplicate the
        // documentation.
        if (s!="nstar_rot::polytrope_solve") {
          
          // Handle multiple entry or add to the main list
          if (list.find(s)!=list.end()) {
            cout << "Multiple entry " << s << endl;
            list_dup.insert(s);
            if (kk==0) {
              cerr << "Duplicate classes not allowed." << endl;
              exit(-1);
            }
            if (kk==2) {
              cerr << "Duplicate files not allowed." << endl;
              exit(-1);
            }
          } else if (ns!=((string)"boost")) {
            list.insert(std::make_pair(s,ns));
          }
          
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
    if (kk==0) {
      cout << "Creating rst files for classes:" << endl;
    } else if (kk==1) {
      cout << "Creating rst files for non-duplicate functions:" << endl;
      cout << endl;
    } else {
      cout << "Creating rst files for files:" << endl;
    }

    // Proceed through the list
    for (std::map<std::string,std::string>::iterator it=list.begin();
	 it!=list.end();it++) {

      string s=it->first, ns=it->second;

      // Make sure it's not a duplicate
      if (list_dup.find(s)==list_dup.end()) {

	// Open the rst file
	string fname_out="class/";
	if (kk==1) fname_out="function/";
	else if (kk==2) fname_out="file/";
	fname_out+=s+".rst";
        if (s.length()<4 || s[0]!='o' || s[1]!='p' ||
            s[2]!='e' || s[3]!='r') {
          fname_out=underscoreify(fname_out);
        }
	ofstream fout(fname_out);

	// Title
	string head="Class "+s;
	// If we're in "class mode" and the namespace is not empty,
	// then add the namespace to the header
	if (kk==0 && ns.length()>0) {
	  head+=" ("+ns+")";
	} else if (kk==1) {
	  // If we're in function mode, then switch the to a
	  // function header
	  head="Function "+s;
	  if (ns.length()>0) {
	    head+=" ("+ns+")";
	  }
	} else if (kk==2) {
	  head="File "+s;
	}
	fout << head << endl;
	for(size_t i=0;i<head.length();i++) {
	  fout << "=";
	}
	fout << endl;
	fout << endl;

	// Links to top-level 
	if (context==((string)"main")) {
	  if (kk==0) {
	    fout << ":ref:`O2scl <o2scl>` : :ref:`Class List`\n" << endl;
	    fout << ".. _" << s << ":\n" << endl;
	  } else if (kk==1) {
	    fout << ":ref:`O2scl <o2scl>` : :ref:`Function List`\n" << endl;
	  } else {
	    fout << ":ref:`O2scl <o2scl>` : :ref:`File List`\n" << endl;
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
	
	// Output the class or function directive
	if (kk==0) {
	  if (ns.length()>0) {
	    fout << ".. doxygenclass:: " << ns << "::" << s << endl;
	  } else {
	    fout << ".. doxygenclass:: " << s << endl;
	  }
	} else if (kk==1) {
	  if (context==((string)"part")) {
	    fout << ".. doxygenfunction:: " << s << endl;
	  } else {
	    fout << ".. doxygenfunction:: " << ns << "::" << s << endl;
	  }
	} else {
	  fout << ".. doxygenfile:: " << s << endl;
	}

	// Close the file
	fout.close();

      } else {
	// Skip items in the list of duplicates
	cout << "Skipping " << it->first << " because it is a duplicate."
             << endl;
      }
    }

    // -------------------------------------------------------
    // For functions, create the rst files for duplicate entries
    
    if (kk==1) {

      cout << endl;
      cout << "---------------------------------------------------" << endl;
      cout << "Creating rst files for duplicate functions:" << endl;
      cout << endl;
    
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
	fin.open("../xml/namespaceo2scl__auto__format.xml");
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

        if (s!="operator<<") {
        
          // Open the rst file
          string fname_out="function/";
          fname_out+=s+".rst";
          if (s.length()<4 || s[0]!='o' || s[1]!='p' ||
              s[2]!='e' || s[3]!='r') {
            fname_out=underscoreify(fname_out);
          }
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
                    cout << "ns_file[i]: " << ns_file[i] << endl;
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
    
  }
  
#endif
  
  return 0;
}
